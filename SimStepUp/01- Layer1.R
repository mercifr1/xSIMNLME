#'####################################################
#'
#'   Simulation study to show the impact of w/i patient
#'   dose adjustment (e.g. during a dose-escalation study)
#'   on signal/noise ratio for SLD values
#'   
#'   Layer1: Individual level profile
#'   
#'   Francois Mercier 
#'   Q1-2020
#'   
#'####################################################



library(tidyverse)
library(mrgsolve)
library(patchwork)


#' Empirical PKPD model
#'----------------------------------------------------

code <- '
$PARAM TVCL=5, WT=70, DOSE=0, T=0, 
TVLBS=4.2, TVLKG=0.5, TVLKS=0.5, TVLPHI=0.7

$OMEGA 1 1.14 0.77 0.5 0.83

$SIGMA 7 

$PRED
double CL = TVCL*pow(WT/70,0.75)*exp(ETA(1));
double KG = exp(TVLKG)*exp(ETA(2));
double KS = exp(TVLKS+(AUC/10))*exp(ETA(3));
double PHI = exp(TVLPHI*(1+ETA(4))+(AUC*10))/(1+exp(TVLPHI*(1+ETA(4))+(AUC*10)));
double BS = exp(TVLBS)*exp(ETA(5));

capture AUC = DOSE/CL;
capture Y = (BS*((1-PHI)*exp(KG*T)+PHI*exp(-KS*T)))+EPS(1);
capture BAS = BS;
'

tk1<-mcode_cache("tk1", code)


#' Time and design matrices
#'----------------------------------------------------

#' Need 20 patients per dose level
#' To each ID corresponds a weight; CL is fixed to 5 for all patients
#' 

func_days<-function(nTAm, skel){
  rTAdys<-jitter(skel, amount=10)         #' Jitter TA days by 'amount' days
  TAyr<-round(rTAdys/365, 3)              #' Round to 3 digit years
  return(data.frame(TAnum=0:nTA, TAyr))
} 

nCoh<-5
nInd<-20
nTA<-7
dsnmat<-expand_grid(Coh=1:nCoh, ID=1:nInd) %>%
  mutate(UID=paste0("C", Coh, "ID", ID))

skel_eo6w<-(0:nTA)*(7*6)                    #' Basic TA skeleton: TA every 6 weeks
skel_exp<-7*c(0, 6, 12, 20, 28, 36, 44, 52) #' every 6 weeks, then 8 then 12.

indf0<-dsnmat %>% 
  split(.$UID) %>% 
  map_dfr(., ~func_days(nTA=nTA, skel=skel_eo6w), .id="UID")

indf1<-dsnmat %>% 
  split(.$UID) %>% 
  map_dfr(., ~func_days(nTA=nTA, skel=skel_exp), .id="UID")


#' Flat dosing
#' -----------
doses<-c(0, 0.3, 1, 3, 10)
flt.dose<-rep(doses, each=(nInd+nCoh)*(nTA+1))

flt<-left_join(indf0, dsnmat, by="UID") %>%
  mutate(T=TAyr, 
         ID=as.numeric(Coh*100+ID),
         DOSE=case_when(Coh==1 ~0,
                        Coh==2 ~0.3,
                        Coh==3 ~1,
                        Coh==4 ~3,
                        Coh==5 ~10,
                        TRUE~9999),
         WT=rep(exp(rnorm(nInd*nCoh, log(75), .1)), each=nTA+1))


#' StepUp dosing
#' -----------
dose0<-c(0, 0.3, 1, 3, 10)
dose1<-c(0.3, 1, 3, 10, 30)

tswitch<-tibble(TSWITCH=sample(seq(from=3, to=7), nCoh*nInd, replace=T),
                    UID=as.character(unique(dsnmat$UID)))

stepup0<-left_join(indf0, dsnmat, by="UID") %>%
  left_join(., tswitch, by="UID") %>%
  mutate(
    T=TAyr, 
    ID=as.numeric(Coh*100+ID),
    DOSE0=rep(dose0, each=nInd*(nTA+1)),
      DOSE1=rep(dose1, each=nInd*(nTA+1)),
    WT=rep(exp(rnorm(nInd*nCoh, log(75), .1)), each=nTA+1)) %>%
  group_by(ID) %>%
  mutate(DOSE=ifelse(DOSE0!=0 & TAnum>=TSWITCH, DOSE1, DOSE0)) %>%
  ungroup()

#' StepUp dosing and alternative TA sampling scheme
#' -----------
stepup1<-left_join(indf1, dsnmat, by="UID") %>%
  left_join(., tswitch, by="UID") %>%
  mutate(
    T=TAyr, 
    ID=as.numeric(Coh*100+ID),
    DOSE0=rep(dose0, each=nInd*(nTA+1)),
         DOSE1=rep(dose1, each=nInd*(nTA+1)),
    WT=rep(exp(rnorm(nInd*nCoh, log(75), .1)), each=nTA+1)) %>%
  group_by(ID) %>%
  mutate(DOSE=ifelse(DOSE0!=0 & TAnum>=TSWITCH, DOSE1, DOSE0)) %>%
  ungroup()


#' Simulate
#'----------------------------------------------------

set.seed(1789)
snd.flt<-mrgsim_d(tk1, as.data.frame(flt), carry.out="T,DOSE,ID") %>%
  as.data.frame() %>%
  filter(BAS<300) %>% mutate(DESN=" Flat dosing") %>%
  left_join(., flt, by=c("ID", "T", "DOSE"))

#set.seed(1789)
#snd.stepup0<-mrgsim_d(tk1, as.data.frame(stepup0), carry.out="T,DOSE,ID") %>%
#  as.data.frame() %>%
#  filter(BAS<300) %>% 
#  mutate(DESN="Dose Increase (0)") %>%
#  left_join(., stepup0, by=c("ID", "T", "DOSE")) %>%
#  select(-DOSE0, -DOSE1, -TSWITCH)

set.seed(1789)
snd.stepup1<-mrgsim_d(tk1, as.data.frame(stepup1), carry.out="T,DOSE,ID") %>%
  as.data.frame() %>%
  filter(BAS<300) %>% 
  mutate(DESN="Dose Increase") %>%
  left_join(., stepup1, by=c("ID", "T", "DOSE")) %>%
  select(-DOSE0, -DOSE1, -TSWITCH) 

#snd.mix<-rbind(snd.flt, snd.stepup0, snd.stepup1) %>%
snd.mix<-rbind(snd.flt, snd.stepup1) %>%
  filter(Y<600) %>%
  mutate(YC=round(ifelse(Y<5, 2.5, Y), 1))


#' Graphical display
#'----------------------------------------------------

#' ggplot(snd.mix, aes(DOSE, AUC))+
#'   geom_point()+
#'   theme_minimal()

#' Fig 01: Dose-dependent trend
ggplot(snd.mix, aes(T, YC))+
  geom_line(aes(group=ID), lwd=.8, alpha=0.3)+
  #  geom_point(alpha=0.3)+
  facet_grid(DESN~Coh)+
  scale_x_continuous("Time (year)")+
  scale_y_continuous("SLD (mm)")+
  theme_minimal()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank())


#' Fig 02: Dose-dependent trend
#'         + sampling scheme
ggplot(snd.mix, aes(T, YC))+
  #  geom_line(aes(group=ID), colour="grey", alpha=0.3)+
  geom_point(size=1, alpha=0.3)+
  facet_grid(DESN~Coh)+
  scale_x_continuous("Time (year)")+
  scale_y_continuous("SLD (mm)")+
  theme_minimal()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank())


#' Shaving off:
#' - 20% progression from nadir
#' - 20% PD due to new lesion (assuming exponential)
#'----------------------------------------------------

tnewles<-tibble(TNEWL=2*rexp(nInd*nCoh, rate=1.5),
                UID=as.character(unique(dsnmat$UID)))
#hist(tnewles$TNEWL)

shv<-left_join(snd.mix, tnewles, by=c("UID")) %>%
  mutate(SHAV=ifelse(T>TNEWL, 1, 0))


#' Graphical display retained in the end
#'----------------------------------------------------

subtxt<-paste0("20 patients per cohort; with doses 0, 0.3, 1, 3, 10 mg and including layers L0 to L4.")
captxt<-paste0("R packages: tidyverse, mrgsolve. mercief3/20200504")
g01_Blur<-shv %>%
  filter(SHAV==0) %>%
  ggplot(., aes(T, YC))+
  geom_line(aes(group=ID), colour="grey70", lwd=.8, alpha=0.2)+
  geom_point(size=2, colour="darkblue", alpha=0.3)+
  #  geom_point(alpha=0.3)+
  facet_grid(~DESN)+
#  facet_grid(DESN~Coh)+
  scale_x_continuous("Time (year)", breaks=c(0, 0.5))+
  scale_y_continuous("SLD (mm)", trans="log", breaks=c(1, 3, 10, 30, 100, 300))+
  labs(title="100 patients on both sides",
       subtitle=subtxt,
       caption=captxt)+
  theme_minimal()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank())
g01_Blur
ggsave("./SimStepUp/Fig1.png", g01_Blur, width=7, height=4, dpi=300)

subtxt<-paste0("20 patients per cohort; including layers L0 to L4.")
captxt<-paste0("R packages: tidyverse, mrgsolve. mercief3/20200504")
g02_Blur<-shv %>%
  filter(SHAV==0) %>%
  ggplot(., aes(T, YC))+
  geom_line(aes(group=ID), colour="grey70", lwd=.8, alpha=0.2)+
  geom_point(size=2, colour="darkblue", alpha=0.3)+
  #  geom_point(alpha=0.3)+
  facet_grid(DESN~Coh, 
             labeller=labeller(Coh=c('1'="0 mg", '2'="0.3 mg", '3'="1 mg", '4'="3 mg", '5'="10 mg")))+
  scale_x_continuous("Time (year)", breaks=c(0, 0.5))+
  scale_y_continuous("SLD (mm)", trans="log", breaks=c(1, 3, 10, 30, 100, 300))+
  labs(title="Split per dose cohort",
       subtitle=subtxt,
       caption=captxt)+
  theme_minimal()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank())
g02_Blur
ggsave("./SimStepUp/Fig2.png", g02_Blur, width=7, height=4, dpi=300)

subtxt<-paste0("20 patients per cohort; including layers L0 to L4.")
captxt<-paste0("R packages: tidyverse, mrgsolve. mercief3/20200504")
g03_Blur<-shv %>%
  filter(SHAV==0) %>%
  ggplot(., aes(T, YC))+
  geom_line(aes(group=ID), colour="grey70", lwd=.8, alpha=0.2)+
  geom_point(size=2, colour="darkblue", alpha=0.3)+
  geom_line(data=snd.fixed, aes(T, Y), lwd=1.7, colour="orange2")+
  #  geom_point(alpha=0.3)+
  facet_grid(DESN~Coh, 
             labeller=labeller(Coh=c('1'="0 mg", '2'="0.3 mg", '3'="1 mg", '4'="3 mg", '5'="10 mg")))+
  scale_x_continuous("Time (year)", breaks=c(0, 0.5))+
  scale_y_continuous("SLD (mm)", trans="log", breaks=c(1, 3, 10, 30, 100, 300))+
  labs(title="Split per dose cohort; Typical trend materialized in orange",
       subtitle=subtxt,
       caption=captxt)+
  theme_minimal()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank())
g03_Blur
ggsave("./SimStepUp/Fig3.png", g03_Blur, width=7, height=4, dpi=300)


