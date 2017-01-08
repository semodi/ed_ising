dat <- read.csv("~/Documents/Physics/Würzburg/2016/ED_ising/outs/E0_N.dat",head = F)
library(dplyr)
library(ggplot2)
dat <- dat %>% rename(N = V2) %>% rename(E0 = V1)
dat <- dat %>% mutate(set = 1)
dat$E0 <- dat$E0 -100
linear_model <- lm(E0 ~ N,dat);linear_model

dat2 <- read.csv("~/Documents/Physics/Würzburg/2016/ED_ising/outs/E0_N_2.dat",head = F)
dat2 <- dat2 %>% rename(N = V2) %>% rename(E0 = V1)
dat2 <- dat2 %>% mutate(set= 2)
dat2$E0 <- dat2$E0-40
dat <- rbind(dat,dat2)

dat$set <- as.factor(dat$set)

dat %>% filter(set == 1) %>% ggplot(aes(N,E0))+geom_point()+
  ggtitle(expression(atop("Ground state energy w.r.t system size", atop(italic("slope = -1.27"), ""))))

dat %>% filter(set == 1) %>% ggplot(aes(N,E0/N))+geom_point()+
  ggtitle("Ground state energy density w.r.t system size")
