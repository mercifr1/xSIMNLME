#####################################################
#   SteinFojo simulate
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
m02saem<-function(psi, id, xidep) {
  tim<-xidep[, 1]
  lBSLD<-psi[id, 1]
  lKG<-psi[id, 2]
  lKS<-psi[id, 3]
  GRW=exp((exp(lKG)/1000)*tim)
  SHR=exp(-(exp(lKS)/1000)*tim)
  ypred<-exp(lBSLD)*(GRW+SHR-1)
  return(ypred)
}
m02.1.mod<-saemixModel(model=m02saem, description="SteinFojo",
                       psi0=matrix(c(4.6, 1.1, 1.1), ncol=3, byrow=TRUE, dimnames=list(NULL, c("lBSLD", "lKG", "lKS"))),
                       transform.par=c(0, 0, 0), fixed.estim=c(1, 1, 1),
                       covariance.model=matrix(c(1, 0, 0,
                                                 0, 1, 0,
                                                 0, 0, 1), ncol=3, byrow=TRUE),
                       omega.init=matrix(c(1, 0, 0,
                                           0, 1, 0,
                                           0, 0, 1), ncol=3, byrow=TRUE),
                       error.model="constant")
saem.m02.1<-saemix(m02.1.mod, din, op)
#' ####### Fit with SAEMIX ########
#' ===============================================
#' 
#' STEP I.
#' Sampling from uncertainty matrix (aka priors in Bayesian terminology)
#' 
#' ===============================================
#' Typical/population-level values for thetas (means)
th.lbsld<-4.7
th.lkg<-(-0.093)
th.lks<-1.6
#' Typical/population-level values for omegas (IIV)
om.lbsld<-sqrt(0.6)
om.lkg<-sqrt(0.35)
om.lks<-sqrt(0.44)
#' Uncertainty (SE) for thetas
un.th.lbsld<-0.09*0.09        #SE^2 (SD^2=var)  
un.th.lkg<-0.1*0.1    
un.th.lks<-0.1*0.1    
#' Uncertainty (SE) for omegas
un.om.lbsld<-0.1*0.1    
un.om.lkg<-0.09*0.09
un.om.lks<-0.1*0.1
#' Correlations between thetas or omegas
um11<-sqrt(un.th.lbsld); um22<-sqrt(un.th.lkg); um33<-sqrt(un.th.lks); 
uv11<-sqrt(un.om.lbsld); uv22<-sqrt(un.om.lkg); uv33<-sqrt(un.om.lks); 
#' if we had correlation, e.g. 0.1 between TVbsl and TVkg, this would be implemented like that:
#' ucorm11m22<-0.1
#' ucovm12<-ucorm11m22*sqrt(um11)*sqrt(um22)
Umat<-matrix(c(um11, 0, 0, 0, 0, 0,
               0, um22, 0, 0, 0, 0,
               0, 0, um33, 0, 0, 0,
               0, 0, 0, uv11, 0, 0,
               0, 0, 0, 0, uv22, 0,
               0, 0, 0, 0, 0, uv33), nrow=6)
#' Uncertainty vector
#' sigma is the vcov matrix
#' method="svd" for single value decomposition
Params<-rmvnorm(5, 
                mean=c(th.lbsld, th.lkg, th.lks, om.lbsld, om.lkg, om.lks), 
                sigma=Umat, 
                method="svd")
#' set_colnames is a convenient function from 'magrittr'
dfParams<-data.frame(Rep=1:5, Params) %>%
  set_colnames(., c("Replicate", "th.lbsld", "th.lkg", "th.lks", "om.lbsld", "om.lkg", "om.lks")) %>%
  mutate(futime=round(runif(5, min=1, max=500), 0))
#' ===============================================
#' 
#' STEP II.
#' Generate regression lines at the population level i.e. w/o IIV
#' 
#' ===============================================
nlm01_pred<-function(mlbsld, mlkg, mlks, t){
  y=exp(mlbsld)*(exp((exp(mlkg)/1000)*t)+exp(-(exp(mlks)/1000)*t)-1)
  return(y)
}
znlm01_preds<-dfParams %>% 
  split(.$Rep) %>% 
  map_dfr(., ~nlm01_pred(.$th.lbsld, .$th.lkg, .$th.lks, t=-2:500))
#str(zlm01_preds)
znlm01p<-gather(znlm01_preds) %>% 
  rename(repID=key, ppred=value) %>%
  mutate(repID=as.numeric(repID), t=rep(-2:500, 5))
#str(zlm01p)
gfojo01<-ggplot(znlm01p, aes(t/7, ppred))+
  geom_line(aes(group=repID))+
  scale_x_continuous("Weeks", breaks=8*0:16)+
  scale_y_continuous("Response", limits=c(0, 200))+
  theme_minimal()
ggsave("./Figs/Fojo01.png", plot=gfojo01, width=5, height=5/sqrt(2))

#' ===============================================
#' 
#' STEP III.
#' Generate regression lines at the population level i.e. w/o IIV
#' with different follow up duration across replicates
#' 
#' ===============================================
nlm02_pred<-function(mlbsld, mlkg, mlks, mtend){
  t=seq(from=-1, to=mtend, by=1)
  ppred=exp(mlbsld)*(exp((exp(mlkg)/1000)*t)+exp(-(exp(mlks)/1000)*t)-1)
  return(data.frame(ppred, t))
}
znlm02_preds<-dfParams %>% 
  split(.$Rep) %>% 
  map(., ~nlm02_pred(.$th.lbsld, .$th.lkg, .$th.lks, .$futime))
#str(znlm02_preds)
znlm02p<-bind_rows(znlm02_preds, .id="repID")
ggplot(znlm02p, aes(t/7, ppred))+
  geom_line(aes(group=repID))+
  theme_minimal()
#' ===============================================
#' 
#' STEP IV.
#' Generate regression lines at the individual level i.e. with IIV
#' with different follow up duration across replicates
#' without residual error
#' 
#' ===============================================
Nind=40
Tmax=500
etas<-function(omlbsld, omlkg, omlks, Nind){
  etalbsld<-rnorm(Nind, 0, omlbsld)
  etalkg<-rnorm(Nind, 0, omlkg)
  etalks<-rnorm(Nind, 0, omlks)
  return(data.frame(etalbsld, etalkg, etalks))
}
#etas(om.lbsld, om.lks, 4)
desmat<-data.frame(ID=1:Nind, itend=runif(Nind, 1, Tmax)) %>% 
  split(.$ID)%>%
  map(., ~expand.grid(t=-1:.$itend))
desmatp<-bind_rows(desmat, .id="ID") %>% mutate(ID=as.numeric(as.character(ID)))
nlme01_Level2<-function(mlbsld, mlkg, mlks, omlbsld, omlkg, omlks, Nind, desmatrix){
  myetas<-etas(omlbsld, omlkg, omlks, Nind)
  myparms<-myetas %>% mutate(mlbsld, mlkg, mlks, ID=1:Nind)
  myalls<-left_join(desmatrix, myparms, by="ID")
  return(myalls)
}
znlme01_Level12<-dfParams %>% 
  split(.$Rep) %>% 
  map(., ~nlme01_Level2(mlbsld=.$th.lbsld, mlkg=.$th.lkg, mlks=.$th.lks, 
                        omlbsld=.$om.lbsld, omlkg=.$om.lkg, omlks=.$om.lks, Nind=Nind, desmatrix=desmatp))
znlm02_Level12p<-bind_rows(znlme01_Level12, .id="repID")
odf<-znlm02_Level12p %>% 
  mutate(ipred=exp(mlbsld+etalbsld)*(exp((exp(mlkg+etalkg)/1000)*t)+exp(-(exp(mlks+etalks)/1000)*t)-1))
ggplot(odf, aes(t/7, ipred))+
  geom_line(aes(group=ID, colour=as.factor(ID)))+
  facet_wrap(~repID)+
  scale_colour_viridis_d(guide=F)+
  theme_minimal()
