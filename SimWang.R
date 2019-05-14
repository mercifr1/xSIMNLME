#####################################################
#   SteinFojo Stretched simulate
#   Q2-2019
#   Francois Mercier 
#####################################################
library(magrittr) #' useful for set_colnames
library(tidyverse)
library(mvtnorm)
set.seed(922)
#' ####### Fit with SAEMIX ########
library(saemix)
indf<-readRDS("./indf.rds")
#Retain 80 patients only
set.seed(1005)
idtrain<-sample(indf$UID, 80, replace=F)
trndf<-indf %>% filter(UID %in% idtrain)
din<-saemixData(name.data=trndf, header=T, name.group=c("UID"),
                name.predictors=c("TIME"), name.response=c("SLD2"),
                units=list(x="Days", y="mm"), name.X="TIME") 
op<-list(directory="./saemix", algorithm=c(1, 1), nbiter.saemix=c(3000, 1000), 
         nb.chains=4, seed=0510207, save=F, save.graphs=F)
m04saem<-function(psi, id, xidep) {
  tim<-xidep[, 1]
  lBSLD<-psi[id, 1]
  lKG<-psi[id, 2]
  lKS<-psi[id, 3]
  KGpow<-psi[id, 4]
  GRW=(exp(lKG)^KGpow*tim)
  SHR=exp(-exp(lKS)*tim)
  ypred<-exp(lBSLD)*SHR+GRW
  return(ypred)
}
m04.1.mod<-saemixModel(model=m04saem, description="Wang modified",
                       psi0=matrix(c(4.5, -7, -5, 1), ncol=4, byrow=TRUE, dimnames=list(NULL, c("lBSLD", "lKG", "lKS", "KGpow"))),
                       transform.par=c(0, 0, 0, 0), fixed.estim=c(1, 1, 1, 1),
                       covariance.model=matrix(c(1, 0, 0, 0,
                                                 0, 1, 0, 0,
                                                 0, 0, 1, 0,
                                                 0, 0, 0, 1), ncol=4, byrow=TRUE),
                       omega.init=matrix(c(1, 0, 0, 0,
                                           0, 1, 0, 0,
                                           0, 0, 1, 0,
                                           0, 0, 0, 1), ncol=4, byrow=TRUE),
                       error.model="constant")
saem.m04.1<-saemix(m04.1.mod, din, op)
#' ####### Fit with SAEMIX ########
#' ===============================================
#' 
#' STEP I.
#' Sampling from uncertainty matrix (aka priors in Bayesian terminology)
#' 
#' ===============================================
#' Typical/population-level values for thetas (means)
th.lbsld<-4.5
th.lkg<-(-3.2)
th.lks<-(-5.1)
th.KGpow<-0.8
#' Typical/population-level values for omegas (IIV)
om.lbsld<-sqrt(0.6)
om.lkg<-sqrt(0.6)
om.lks<-sqrt(0.5)
om.KGpow<-sqrt(0.1)
#' Uncertainty (SE) for thetas
un.th.lbsld<-0.09*0.09        #SE^2 (SD^2=var)  
un.th.lkg<-0.75*0.75    
un.th.lks<-0.1*0.1    
un.th.KGpow<-0.18*0.18    
#' Uncertainty (SE) for omegas
un.om.lbsld<-0.1*0.1    
un.om.lkg<-1.3*1.3
un.om.lks<-0.1*0.1
un.om.KGpow<-0.1*0.1
#' Correlations between thetas or omegas
um11<-sqrt(un.th.lbsld); um22<-sqrt(un.th.lkg); um33<-sqrt(un.th.lks); um44<-sqrt(un.th.KGpow); 
uv11<-sqrt(un.om.lbsld); uv22<-sqrt(un.om.lkg); uv33<-sqrt(un.om.lks); uv44<-sqrt(un.om.KGpow); 
#' if we had correlation, e.g. 0.1 between TVbsl and TVkg, this would be implemented like that:
#' ucorm11m22<-0.1
#' ucovm12<-ucorm11m22*sqrt(um11)*sqrt(um22)
Umat<-matrix(c(um11, 0, 0, 0, 0, 0, 0, 0,
               0, um22, 0, 0, 0, 0, 0, 0,
               0, 0, um33, 0, 0, 0, 0, 0,
               0, 0, 0, um44, 0, 0, 0, 0,
               0, 0, 0, 0, uv11, 0, 0, 0,
               0, 0, 0, 0, 0, uv22, 0, 0,
               0, 0, 0, 0, 0, 0, uv33, 0,
               0, 0, 0, 0, 0, 0, 0, uv44), nrow=8)
#' Uncertainty vector
#' sigma is the vcov matrix
#' method="svd" for single value decomposition
Params<-rmvnorm(5, 
                mean=c(th.lbsld, th.lkg, th.lks, th.KGpow, om.lbsld, om.lkg, om.lks, om.KGpow), 
                sigma=Umat, 
                method="svd")
#' set_colnames is a convenient function from 'magrittr'
dfParams<-data.frame(Rep=1:5, Params) %>%
  set_colnames(., c("Replicate", "th.lbsld", "th.lkg", "th.lks", "th.KGpow", 
                    "om.lbsld", "om.lkg", "om.lks", "om.KGpow")) %>%
  mutate(futime=round(runif(5, min=1, max=500), 0))
#' ===============================================
#' 
#' STEP II.
#' Generate regression lines at the population level i.e. w/o IIV
#' 
#' ===============================================
nlm01_pred<-function(mlbsld, mlkg, mlks, mKGpow, t){
  GRW=(exp(mlkg)^mKGpow*t)
  SHR=exp(-exp(mlks)*t)
  y=exp(mlbsld)*SHR+GRW
  return(y)
}
znlm01_preds<-dfParams %>% 
  split(.$Rep) %>% 
  map_dfr(., ~nlm01_pred(.$th.lbsld, .$th.lkg, .$th.lks, .$th.KGpow, t=-2:500))
#str(zlm01_preds)
znlm01p<-gather(znlm01_preds) %>% 
  rename(repID=key, ppred=value) %>%
  mutate(repID=as.numeric(repID), t=rep(-2:500, 5))
#str(zlm01p)
gwang01<-ggplot(znlm01p, aes(t/7, ppred))+
  geom_line(aes(group=repID))+
  scale_x_continuous("Weeks", breaks=8*0:16)+
  scale_y_continuous("Response", limits=c(0, 120))+
  theme_minimal()
ggsave("./Figs/Wang01.png", plot=gwang01, width=5, height=5/sqrt(2))

#' ===============================================
#' 
#' STEP III.
#' Generate regression lines at the population level i.e. w/o IIV
#' with different follow up duration across replicates
#' 
#' ===============================================
nlm02_pred<-function(mlbsld, mlkg, mlks, mKGpow, mtend){
  t=seq(from=-1, to=mtend, by=1)
  GRW=(exp(mlkg)^mKGpow*t)
  SHR=exp(-exp(mlks)*t)
  ppred=exp(mlbsld)*SHR+GRW
  return(data.frame(ppred, t))
}
znlm02_preds<-dfParams %>% 
  split(.$Rep) %>% 
  map(., ~nlm02_pred(.$th.lbsld, .$th.lkg, .$th.lks, .$th.KGpow, .$futime))
#str(znlm02_preds)
znlm02p<-bind_rows(znlm02_preds, .id="repID")
ggplot(znlm02p, aes(t/7, ppred))+
  geom_line(aes(group=repID))+
  theme_minimal()