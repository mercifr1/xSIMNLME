




library(magrittr) #' useful for set_colnames
library(tidyverse)
library(mvtnorm)

set.seed(922)

#' ===============================================
#' 
#' STEP I.
#' Sampling from uncertainty matrix (aka priors in Bayesian terminology)
#' 
#' ===============================================

#' Typical/population-level values for thetas (means)
th.int<-10
th.slo<-2

#' Typical/population-level values for omegas (IIV)
om.int<-3
om.slo<-4

#' Uncertainty (SE) for thetas
un.th.int<-1.1*1.1      #SE^2 (SD^2=var) th.int
un.th.slo<-0.9*0.9    #SE^2 (SD^2=var) th.slo

#' Uncertainty (SE) for omegas
un.om.int<-0.05*0.05    #SE^2 (SD^2=var) om.int
un.om.slo<-0.07*0.07    #SE^2 (SD^2=var) om.slo

#' Correlations between thetas or omegas
corth<-0.1
corom<-0.1

covth<-corth*sqrt(un.th.int)*sqrt(un.th.slo)
covom<-corom*sqrt(un.om.int)*sqrt(un.om.slo)
Umat<-matrix(c(un.th.int, covth, 0, 0,
               covth, un.th.slo, 0, 0,
               0, 0, un.om.int, covom,
               0, 0, covom, un.om.slo), nrow=4)

#' Uncertainty vector
#' sigma is the vcov matrix
#' method="svd" for single value decomposition
Params<-rmvnorm(5, 
  mean=c(th.int, th.slo, om.int, om.slo), 
  sigma=Umat, 
  method="svd")

#' set_colnames is a convenient function from 'magrittr'
dfParams<-data.frame(Rep=1:5, Params) %>%
  set_colnames(., c("Replicate", "th.int", "th.slo", "om.int", "om.slo")) %>%
  mutate(futime=round(runif(5, min=1, max=20), 0))


#' ===============================================
#' 
#' STEP II.
#' Generate regression lines at the population level i.e. w/o IIV
#' 
#' ===============================================

lm01_pred<-function(mint, mslo, t){
  y=mint+mslo*t
  return(y)
}

zlm01_preds<-dfParams %>% 
  split(.$Rep) %>% 
  map_dfr(., ~lm01_pred(.$th.int, .$th.slo, t=-2:10))
#str(zlm01_preds)
zlm01p<-gather(zlm01_preds) %>% 
  rename(repID=key, ppred=value) %>%
  mutate(repID=as.numeric(repID), t=rep(-2:10, 5))
#str(zlm01p)

ggplot(zlm01p, aes(t, ppred))+
  geom_point()+
  geom_line(aes(group=repID))+
  theme_minimal()


#' ===============================================
#' 
#' STEP III.
#' Generate regression lines at the population level i.e. w/o IIV
#' with different follow up duration across replicates
#' 
#' ===============================================

lm02_pred<-function(mint, mslo, mtend){
  t=seq(from=-1, to=mtend, by=1)
  ppred=mint+mslo*t
  return(data.frame(ppred, t))
}

zlm02_preds<-dfParams %>% 
  split(.$Rep) %>% 
  map(., ~lm02_pred(.$th.int, .$th.slo, .$futime))
#str(zlm02_preds)

zlm02p<-bind_rows(zlm02_preds, .id="repID")

ggplot(zlm02p, aes(t, ppred))+
  geom_point()+
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

Nind=4
Tmax=10

etas<-function(omint, omslo, Nind){
  etaint<-rnorm(Nind, 0, omint)
  etaslo<-rnorm(Nind, 0, omslo)
  return(data.frame(etaint, etaslo))
}
#etas(om.int, om.slo, 4)

desmat<-data.frame(ID=1:Nind, itend=runif(Nind, 1, Tmax)) %>% 
  split(.$ID)%>%
  map(., ~expand.grid(t=-1:.$itend))
desmatp<-bind_rows(desmat, .id="ID") %>% mutate(ID=as.numeric(as.character(ID)))

lme01_Level2<-function(mint, mslo, omint, omslo, Nind, desmatrix){
  myetas<-etas(omint, omslo, Nind)
  myparms<-myetas %>% mutate(mint, mslo, ID=1:Nind)
  myalls<-left_join(desmatrix, myparms, by="ID")
  return(myalls)
}
#zlme01_Level2<-lme01_Level2(mint=12, mslo=2, omint=3, omslo=4, Nind=4, desmatrix=desmatp)

zlme01_Level12<-dfParams %>% 
  split(.$Rep) %>% 
  map(., ~lme01_Level2(mint=.$th.int, mslo=.$th.slo, omint=.$om.int, omslo=.$om.slo, Nind=4, desmatrix=desmatp))
zlm02_Level12p<-bind_rows(zlme01_Level12, .id="repID")

odf<-zlm02_Level12p %>% mutate(ipred=(mint+etaint)+(mslo+etaslo)*t)

ggplot(odf, aes(t, ipred))+
  geom_point(aes(colour=as.factor(ID)))+
  geom_line(aes(group=ID, colour=as.factor(ID)))+
  facet_wrap(~repID)+
  scale_colour_viridis_d(guide=F)+
  theme_minimal()


