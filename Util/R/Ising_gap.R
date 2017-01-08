dat <- read.csv("~/Documents/Physics/Würzburg/2016/ED_ising/outs/diff.dat",head = F)
library(dplyr)
library(ggplot2)

dat <- dat %>% rename(E = V1) %>% rename(N = V2)
Nlist <- unique(dat$N)

dE <- c()

for(i in Nlist){
  min <- dat %>% filter(N == i) %>% select(E) %>% min()
  max <- dat %>% filter(N == i) %>% select(E) %>% max()
  dE <- rbind(dE,max-min)
}
diff <- data.frame(N = Nlist,dE)
diff %>% ggplot(aes(N,dE/N))+geom_point()+geom_smooth(method = "lm", formula = y~I(1/x^2))+
  labs(title="Gap size w.r.t system size")
  
fit <- lm(log(dE) ~ log(N),data = diff)
summary(fit)
diff <- diff %>% mutate(linear = dE*N)

diff %>% ggplot(aes(log(N),log(dE)))+geom_point()+ geom_smooth()+
  labs(title="Gap size * Length w.r.t. system size") +
  xlab("log(N)") + ylab("log(dE)") +
  geom_abline(color = "red",slope = 0, intercept = pi/2,show.legend = T)



