#'####################################################
#'
#'   BP42169
#'   Simulation study 
#'   Francois Mercier 
#'   Q1-2020
#'   
#'####################################################



library(tidyverse)
library(mrgsolve)
library(patchwork)


#'####################################################
#' PART I: Using fixed effect model
#'####################################################

#' Reference model is CEA-TCB
#'----------------------------------------------------

code <- '
$PARAM TVCL=5, WT=70, DOSE=0, T=0, 
TVLBS=4.2, TVLKG=-2.3, TVLKS=1.5, TVLPHI=0.7

$OMEGA 1 1.14 0.77 0.5 0.83

$SIGMA 10 

$PRED
double CL = TVCL*pow(WT/70,0.75)*exp(ETA(1));
double KG = exp(TVLKG)*exp(ETA(2));
double KS = exp(TVLKS+AUC)*exp(ETA(3));
double PHI = exp(TVLPHI*(1+ETA(4))+(AUC*10))/(1+exp(TVLPHI*(1+ETA(4))+(AUC*10)));
double BS = exp(TVLBS)*exp(ETA(5));

capture AUC = DOSE/CL;
capture Y = (BS*((1-PHI)*exp(KG*T)+PHI*exp(-KS*T)))+EPS(1);
capture BAS = BS;
'

tk1<-mcode_cache("tk1", code)

indf<-expand_grid(ID=1:30, DOSE=c(0, 0.3, 1, 3, 10), T=0.1*(0:10)) %>%
  mutate(WT=exp(rnorm(n(), log(80), 1)),
         ID_DOSE=(100*ID+DOSE))


#head(indf)
set.seed(1527)
sndf<-mrgsim_d(tk1, indf, carry.out="T,DOSE,ID_DOSE") %>%
  as.data.frame() %>%
  filter(BAS<300)

ggplot(sndf, aes(DOSE, AUC))+
  geom_point()+
  theme_minimal()
ggplot(sndf, aes(T, Y))+
  geom_line(aes(group=ID_DOSE), alpha=0.3)+
  geom_point(alpha=0.3)+
  facet_wrap(~DOSE, nrow=1)+
  theme_minimal()


ggsave("../aTYRP1TCB/Pgm/Frac.ran.jpg", frac.ran, width=7, height=4, dpi=300)
