#'####################################################
#'
#'   Simulation study to show the impact of w/i patient
#'   dose adjustment (e.g. during a dose-escalation study)
#'   on signal/noise ratio for SLD values
#'   
#'   Layer0: Population level profile
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

$PRED
double CL = TVCL*pow(WT/70, 0.75);
double KG = exp(TVLKG);
double KS = exp(TVLKS+(AUC/10));
double PHI = exp(TVLPHI+(AUC*10))/(1+exp(TVLPHI+(AUC*10)));
double BS = exp(TVLBS);

capture AUC = DOSE/CL;
capture Y = (BS*((1-PHI)*exp(KG*T)+PHI*exp(-KS*T)));
capture BAS = BS;
'

tk0<-mcode_cache("tk0", code)


#' Time matrix
#'----------------------------------------------------

illus_AUC<-data.frame(DOSE=c(0, 0.3, 1, 3, 10)) %>%
  mutate(AUC=DOSE/5) %>%
  mutate(KS=exp(0.5+(AUC/10)),
          PHI=exp(0.7+(AUC*10))/(1+exp(0.7+(AUC*10))))

g01_ybreaks<-c(0, 0.3, 1, 3, 10)
g01_ylabs<-as.character(c(0, "", 1, 3, 10))
g01<-ggplot(illus_AUC, aes(DOSE, AUC))+
  geom_point()+
  geom_line()+
  scale_y_continuous("AUC (unit)", breaks=unique(illus_AUC$AUC))+
  scale_x_continuous("Dose (mg Q3W iv)", breaks=g01_ybreaks, labels=g01_ylabs)+
  theme_minimal()+
  theme(panel.grid.minor=element_blank())

g02_ybreaks<-round(unique(illus_AUC$KS), 2)
g02_ylabs<-as.character(c(g02_ybreaks[1], "", "", g02_ybreaks[4], g02_ybreaks[5]))
g02<-ggplot(illus_AUC, aes(KS, AUC))+
  geom_point()+
  geom_line()+
  scale_y_continuous("", breaks=unique(illus_AUC$AUC))+
  scale_x_continuous("KS (1/time)", breaks=g02_ybreaks, labels=g02_ylabs)+
  theme_minimal()+
  theme(panel.grid.minor=element_blank())

g03_ybreaks<-round(unique(illus_AUC$PHI), 2)
g03_ylabs<-as.character(g03_ybreaks)
captxt<-paste0("R packages: tidyverse. mercief3/20200512")
g03<-ggplot(illus_AUC, aes(PHI, AUC))+
  geom_point()+
  geom_line()+
  scale_y_continuous("", breaks=unique(illus_AUC$AUC))+
  scale_x_continuous("PHI", breaks=g03_ybreaks, labels=g03_ylabs)+
  labs(caption=captxt)+
  theme_minimal()+
  theme(panel.grid.minor=element_blank())
  
gillus<-g01+g02+g03
ggsave("./SimStepUp/Fig0.png", gillus, width=10, height=4, dpi=300)



#' Time matrix
#'----------------------------------------------------

func_days<-function(nTAm, skel){
  rTAdys<-jitter(skel, amount=10)         #' Jitter TA days by 'amount' days
  TAyr<-round(rTAdys/365, 3)              #' Round to 3 digit years
  return(data.frame(TAnum=0:nTA, TAyr))
} 

#' **Nb of Ind must corresponds to nb of doses in this program**
nInd<-5
nTA<-7
skel_eo6w<-(0:nTA)*(7*6)                    #' Basic TA skeleton: TA every 6 weeks
skel_exp<-7*c(0, 6, 12, 20, 28, 36, 44, 52) #' every 6 weeks, then 8 then 12.

indf0<-tibble(ID=1:nInd) %>% 
  split(.$ID) %>% 
  map_dfr(., ~func_days(nTA=nTA, skel=skel_eo6w), .id="ID")

indf1<-tibble(ID=1:nInd) %>% 
  split(.$ID) %>% 
  map_dfr(., ~func_days(nTA=nTA, skel=skel_exp), .id="ID")


#' Design matrix: Flat dose
#'----------------------------------------------------

#' To each ID corresponds a different dose e.g. 5 ID = 5 doses
#' and to each ID corresponds a weight
doses<-c(0, 0.3, 1, 3, 10)
flt<-indf0 %>%
  mutate(T=TAyr, ID=as.numeric(as.character(ID)),
         DOSE=rep(doses, each=nTA+1),
         WT=rep(exp(rnorm(nInd, log(75), .1)), each=nTA+1))


#' Design matrix: Within-pat dose increase
#' ---------------------------------

#' To each ID corresponds a different dose e.g. 5 ID = 5 doses
#' and to each ID corresponds a weight
dose0<-c(0, 0.3, 1, 3, 10)
dose1<-c(0.3, 1, 3, 10, 30)
tswitch<-c(999, 4, 4, 3, 3)
stepup0<-indf0 %>%
  mutate(T=TAyr, ID=as.numeric(as.character(ID)),
         TSWITCH=rep(tswitch, each=nTA+1),
         DOSE0=rep(dose0, each=nTA+1),
         DOSE1=rep(dose1, each=nTA+1),
         WT=rep(exp(rnorm(nInd, log(75), .1)), each=nTA+1)) %>%
  group_by(ID) %>%
    mutate(DOSE=ifelse(row_number()>TSWITCH, DOSE1, DOSE0)) %>%
  ungroup()


#' Design matrix: Within-pat dose increase
#'      with alternative TA sampling scheme
#' ---------------------------------

#' To each ID corresponds a different dose e.g. 5 ID = 5 doses
#' and to each ID corresponds a weight
stepup1<-indf1 %>%
  mutate(T=TAyr, ID=as.numeric(as.character(ID)),
         TSWITCH=rep(tswitch, each=nTA+1),
         DOSE0=rep(dose0, each=nTA+1),
         DOSE1=rep(dose1, each=nTA+1),
         WT=rep(exp(rnorm(nInd, log(75), .1)), each=nTA+1)) %>%
  group_by(ID) %>%
  mutate(DOSE=ifelse(row_number()>TSWITCH, DOSE1, DOSE0)) %>%
  ungroup()

 
#' Simulate
#'----------------------------------------------------

set.seed(1789)
snd.flt<-mrgsim_d(tk0, as.data.frame(flt), carry.out="T,DOSE") %>%
  as.data.frame() %>%
  filter(BAS<300) %>% mutate(DESN=" Flat dosing")

#set.seed(1789)
#snd.stepup0<-mrgsim_d(tk0, as.data.frame(stepup0), carry.out="T,DOSE") %>%
#  as.data.frame() %>%
#  filter(BAS<300) %>% 
#  mutate(DESN="Dose Increase (0)")

set.seed(1789)
snd.stepup1<-mrgsim_d(tk0, as.data.frame(stepup1), carry.out="T,DOSE") %>%
  as.data.frame() %>%
  filter(BAS<300) %>% 
  mutate(DESN="Dose Increase")

#snd<-rbind(snd.flt, snd.stepup0, snd.stepup1)
snd.fixed<-rbind(snd.flt, snd.stepup1) %>%
  mutate(Coh=ID) %>%
  filter(Y<1000)


#' Graphical display
#'----------------------------------------------------

#' ggplot(snd.flt, aes(DOSE, AUC))+
#'   geom_point()+
#'   theme_minimal()

#' Fig 01: Dose-dependent trend
ggplot(snd.fixed, aes(T, Y))+
  geom_line(aes(group=Coh), lwd=1)+
#  geom_point(alpha=0.3)+
  facet_grid(DESN~Coh)+
  scale_x_continuous("Time (year)")+
  scale_y_continuous("SLD (mm)")+
  theme_minimal()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank())


#' Fig 02: Dose-dependent trend
#'         + sampling scheme
ggplot(snd, aes(T, Y))+
#  geom_line(aes(group=ID), colour="grey", alpha=0.3)+
  geom_point(size=2)+
  facet_grid(DESN~ID)+
  scale_x_continuous("Time (year)")+
  scale_y_continuous("SLD (mm)")+
  theme_minimal()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank())


