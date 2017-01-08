readdat <- read.csv("~/Documents/Physics/Würzburg/2016/ED_ising/out.dat",head = F)
library(dplyr)
library(ggplot2)

readdat <- readdat%>% rename(eigenv = V1) %>% rename(momentum = V2) %>% rename(parity = V3)
readdat$parity <- as.factor(readdat$parity)
glimpse(readdat)

L <- 16

leading <- readdat %>% group_by(momentum,parity) %>% filter(min_rank(eigenv) <= 5)

#offset <- 1/24
#scaling <- L/(2*pi)
offset <- min(leading$eigenv)
scaling <- 1/8*1/(sort(leading$eigenv)[2]-sort(leading$eigenv)[1])
leading <- leading %>% mutate(eigenv0 = scaling*(eigenv-offset))


qplot(leading$momentum,leading$eigenv0,color = leading$parity,alpha = .2)


