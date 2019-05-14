#####################################################
#   Claret-Bruno simulate
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
m03saem<-function(psi, id, xidep) {
  tim<-xidep[, 1]
  lBSLD<-psi[id, 1]
  lKG<-psi[id, 2]
  lKS<-psi[id, 3]
  lLAM<-psi[id, 4]
  GRW=(exp(lKG)*tim)
  SHR=(exp(lKS)/exp(lLAM))*(1-exp(-exp(lLAM)*tim))
  ypred<-exp(lBSLD)*exp(GRW-SHR)
  return(ypred)
}
m03.1.mod<-saemixModel(model=m03saem, description="ClaretBruno",
                       psi0=matrix(c(4.5, -7, -5, -5), ncol=4, byrow=TRUE, dimnames=list(NULL, c("lBSLD", "lKG", "lKS", "lLAM"))),
                       transform.par=c(0, 0, 0, 0), fixed.estim=c(1, 1, 1, 1),
                       covariance.model=matrix(c(1, 1, 1, 1,
                                                 1, 1, 1, 1,
                                                 1, 1, 1, 1,
                                                 1, 1, 1, 1), ncol=4, byrow=TRUE),
                       omega.init=matrix(c(1, 0, 0, 0,
                                           0, 1, 0, 0,
                                           0, 0, 1, 0,
                                           0, 0, 0, 1), ncol=4, byrow=TRUE),
                       error.model="constant")
saem.m03.1<-saemix(m03.1.mod, din, op)
#' ####### Fit with SAEMIX ########
#' ===============================================
#' 
#' STEP I.
#' Sampling from uncertainty matrix (aka priors in Bayesian terminology)
#' 
#' ===============================================
#' Typical/population-level values for thetas (means)
th.lbsld<-4.5
th.lkg<-(-7.2)
th.lks<-(-5.1)
th.llam<-(-4.9)
#' Typical/population-level values for omegas (IIV)
om.lbsld<-sqrt(0.6)
om.lkg<-sqrt(0.7)
om.lks<-sqrt(0.37)
om.llam<-sqrt(0.81)
#' Uncertainty (SE) for thetas
un.th.lbsld<-0.09*0.09        #SE^2 (SD^2=var)  
un.th.lkg<-0.335*0.335    
un.th.lks<-0.1*0.1    
un.th.llam<-0.208*0.208    
#' Uncertainty (SE) for omegas
un.om.lbsld<-0.1*0.1    
un.om.lkg<-0.51*0.51
un.om.lks<-0.1*0.1
un.om.llam<-0.21*0.21
#' Correlations between thetas or omegas
um11<-sqrt(un.th.lbsld); um22<-sqrt(un.th.lkg); um33<-sqrt(un.th.lks); um44<-sqrt(un.th.llam); 
uv11<-sqrt(un.om.lbsld); uv22<-sqrt(un.om.lkg); uv33<-sqrt(un.om.lks); uv44<-sqrt(un.om.llam); 
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
                mean=c(th.lbsld, th.lkg, th.lks, th.llam, om.lbsld, om.lkg, om.lks, om.llam), 
                sigma=Umat, 
                method="svd")
#' set_colnames is a convenient function from 'magrittr'
dfParams<-data.frame(Rep=1:5, Params) %>%
  set_colnames(., c("Replicate", "th.lbsld", "th.lkg", "th.lks", "th.llam", 
                    "om.lbsld", "om.lkg", "om.lks", "om.llam")) %>%
  mutate(futime=round(runif(5, min=1, max=500), 0))
#' ===============================================
#' 
#' STEP II.
#' Generate regression lines at the population level i.e. w/o IIV
#' 
#' ===============================================
nlm01_pred<-function(mlbsld, mlkg, mlks, mllam, t){
  GRW=(exp(mlkg)*t)
  SHR=(exp(mlks)/exp(mllam))*(1-exp(-exp(mllam)*t))
  y=exp(mlbsld)*exp(GRW-SHR)
  return(y)
}
znlm01_preds<-dfParams %>% 
  split(.$Rep) %>% 
  map_dfr(., ~nlm01_pred(.$th.lbsld, .$th.lkg, .$th.lks, .$th.llam, t=-2:500))
#str(zlm01_preds)
znlm01p<-gather(znlm01_preds) %>% 
  rename(repID=key, ppred=value) %>%
  mutate(repID=as.numeric(repID), t=rep(-2:500, 5))
#str(zlm01p)
gbruno01<-ggplot(znlm01p, aes(t/7, ppred))+
  geom_line(aes(group=repID))+
  scale_x_continuous("Weeks", breaks=8*0:16)+
  scale_y_continuous("Response", limits=c(0, 120))+
  theme_minimal()
ggsave("./Figs/Bruno01.png", plot=gbruno01, width=5, height=5/sqrt(2))

#' ===============================================
#' 
#' STEP III.
#' Generate regression lines at the population level i.e. w/o IIV
#' with different follow up duration across replicates
#' 
#' ===============================================
nlm02_pred<-function(mlbsld, mlkg, mlks, mllam, mtend){
  t=seq(from=-1, to=mtend, by=1)
  GRW=(exp(mlkg)*t)
  SHR=(exp(mlks)/exp(mllam))*(1-exp(-exp(mllam)*t))
  ppred=exp(mlbsld)*exp(GRW-SHR)
  return(data.frame(ppred, t))
}
znlm02_preds<-dfParams %>% 
  split(.$Rep) %>% 
  map(., ~nlm02_pred(.$th.lbsld, .$th.lkg, .$th.lks, .$th.llam, .$futime))
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
etas<-function(omlbsld, omlkg, omlks, omllam, Nind){
  etalbsld<-rnorm(Nind, 0, omlbsld)
  etalkg<-rnorm(Nind, 0, omlkg)
  etalks<-rnorm(Nind, 0, omlks)
  etallam<-rnorm(Nind, 0, omllam)
  return(data.frame(etalbsld, etalkg, etalks, etallam))
}
#etas(om.lbsld, om.lks, 4)
desmat<-data.frame(ID=1:Nind, itend=runif(Nind, 1, Tmax)) %>% 
  split(.$ID)%>%
  map(., ~expand.grid(t=-1:.$itend))
desmatp<-bind_rows(desmat, .id="ID") %>% mutate(ID=as.numeric(as.character(ID)))
nlme01_Level2<-function(mlbsld, mlkg, mlks, mllam, omlbsld, omlkg, omlks, omllam, Nind, desmatrix){
  myetas<-etas(omlbsld, omlkg, omlks, omllam, Nind)
  myparms<-myetas %>% mutate(mlbsld, mlkg, mlks, mllam, ID=1:Nind)
  myalls<-left_join(desmatrix, myparms, by="ID")
  return(myalls)
}
znlme01_Level12<-dfParams %>% 
  split(.$Rep) %>% 
  map(., ~nlme01_Level2(mlbsld=.$th.lbsld, mlkg=.$th.lkg, mlks=.$th.lks, mllam=.$th.lks, 
                        omlbsld=.$om.lbsld, omlkg=.$om.lkg, omlks=.$om.lks, omllam=.$om.llam, 
                        Nind=Nind, desmatrix=desmatp))
znlm02_Level12p<-bind_rows(znlme01_Level12, .id="repID")
odf<-znlm02_Level12p %>% 
  mutate(ipred=exp(mlbsld+etalbsld)*exp(
    (exp(mlkg+etalkg)*t)+
      (exp(mlks+etalks)/exp(mllam+etallam))*(1-exp(-exp(mllam+etallam)*t))
  ))
odfo<-odf %>% filter(ipred<1000)
ggplot(odfo, aes(t/7, ipred))+
  geom_line(aes(group=ID, colour=as.factor(ID)))+
  facet_wrap(~repID)+
  scale_colour_viridis_d(guide=F)+
  theme_minimal()
