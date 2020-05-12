

library(tidyverse)
library(hrbrthemes)

lBS<-4.5
lKG<-0.2
lKS<-1.5
zf<-0.7

BS<-exp(lBS)
KG<-exp(lKG)
KS<-exp(lKS)

TKdf<-tibble(tim=seq(0, 1, length=11)) %>%
  mutate(ypred=BS*((1-zf)*exp(KG*tim)+zf*exp(-KS*tim)))

ggplot(TKdf, aes(tim, ypred))+
  geom_point()+
  geom_line()+
  theme_ipsum()
