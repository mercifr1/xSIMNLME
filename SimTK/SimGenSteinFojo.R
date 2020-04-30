#####################################################
#   Generalized SteinFojo simulate
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
m05saem<-function(psi, id, xidep) {
  tim<-xidep[, 1]
  lBSLD<-psi[id, 1]
  lKG<-psi[id, 2]
  lKS<-psi[id, 3]
  zf<-psi[id, 4]
  GRW=exp((exp(lKG)/1000)*tim)
  SHR=exp(-(exp(lKS)/1000)*tim)
  ypred<-exp(lBSLD)*((1-zf)*GRW+zf*SHR)
  return(ypred)
}
m05.1.mod<-saemixModel(model=m05saem, description="Generalized SteinFojo",
                       psi0=matrix(c(4.6, 1.6, 0.9, 0.8), ncol=4, byrow=TRUE, dimnames=list(NULL, c("lBSLD", "lKG", "lKS", "zf"))),
                       transform.par=c(0, 0, 0, 3), fixed.estim=c(1, 1, 1, 1),
                       covariance.model=matrix(c(1, 0, 0, 0,
                                                 0, 1, 0, 0,
                                                 0, 0, 1, 0,
                                                 0, 0, 0, 1), ncol=4, byrow=TRUE),
                       error.model="constant")
saem.m05.1<-saemix(m05.1.mod, din, op)
#' ####### Fit with SAEMIX ########
#' ===============================================
#' 
#' STEP I.
#' Sampling from uncertainty matrix (aka priors in Bayesian terminology)
#' 
#' ===============================================
#' Typical/population-level values for thetas (means)
th.lbsld<-4.47
th.lkg<-0.79
th.lks<-1.72
th.zf<-0.81
#' Typical/population-level values for omegas (IIV)
om.lbsld<-sqrt(0.61)
om.lkg<-sqrt(0.29)
om.lks<-sqrt(0.74)
om.zf<-sqrt(0.81)
#' Uncertainty (SE) for thetas
un.th.lbsld<-0.093*0.093        #SE^2 (SD^2=var)  
un.th.lkg<-0.199*0.199    
un.th.lks<-0.141*0.141    
un.th.zf<-0    
#' Uncertainty (SE) for omegas
un.om.lbsld<-0.1*0.1    
un.om.lkg<-0.16*0.16
un.om.lks<-0.17*0.17
un.om.zf<-0
#' Correlations between thetas or omegas
um11<-sqrt(un.th.lbsld); um22<-sqrt(un.th.lkg); um33<-sqrt(un.th.lks); um44<-sqrt(un.th.zf); 
uv11<-sqrt(un.om.lbsld); uv22<-sqrt(un.om.lkg); uv33<-sqrt(un.om.lks); uv44<-sqrt(un.om.zf); 
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
                mean=c(th.lbsld, th.lkg, th.lks, th.zf, om.lbsld, om.lkg, om.lks, om.zf), 
                sigma=Umat, 
                method="svd")
#' set_colnames is a convenient function from 'magrittr'
dfParams<-data.frame(Rep=1:5, Params) %>%
  set_colnames(., c("Replicate", "th.lbsld", "th.lkg", "th.lks", "th.zf", "om.lbsld", "om.lkg", "om.lks", "om.zf")) %>%
  mutate(futime=round(runif(5, min=1, max=500), 0))
#' ===============================================
#' 
#' STEP II.
#' Generate regression lines at the population level i.e. w/o IIV
#' 
#' ===============================================
nlm01_pred<-function(mlbsld, mlkg, mlks, mzf, t){
  y=exp(mlbsld)*((1-mzf)*exp((exp(mlkg)/1000)*t)+mzf*exp(-(exp(mlks)/1000)*t))
  return(y)
}
znlm01_preds<-dfParams %>% 
  split(.$Rep) %>% 
  map_dfr(., ~nlm01_pred(.$th.lbsld, .$th.lkg, .$th.lks, .$th.zf, t=-2:500))
#str(zlm01_preds)
znlm01p<-gather(znlm01_preds) %>% 
  rename(repID=key, ppred=value) %>%
  mutate(repID=as.numeric(repID), t=rep(-2:500, 5))
#str(zlm01p)
ggenfojo01<-ggplot(znlm01p, aes(t/7, ppred))+
  geom_line(aes(group=repID))+
  scale_x_continuous("Weeks", breaks=8*0:16)+
  scale_y_continuous("Response", limits=c(0, 200))+
  theme_minimal()
ggsave("./Figs/GenFojo01.png", plot=ggenfojo01, width=5, height=5/sqrt(2))