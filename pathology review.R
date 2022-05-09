
rm(list = ls())
### load libraries
library(lme4)
library(lattice)
library(tidyverse)
library(magrittr)
### load dataset
df = read.csv("/Users/FBINSATTER/Desktop/ChRCC- paper 2/Pathology_review2.1.csv") ## file path
colnames(df)
df %<>% select(1:8)
## view data    
head(df)
### develop model
lmm = lme4::lmer(OG.DX ~ H.E + 
                   (1 | Scorer), data = df)
lmm
stargazer::stargazer(lmm, type = "text",
                     digits = 3,
                     star.cutoffs = c(0.05, 0.01, 0.001),
                     digit.separator = "")


lmm2 = lme4::lmer(OG.DX ~ H.E + 
                    S100A1.positive + 
                    S100A1.Staining + 
                    (1 | Scorer), data = df)

lmm3 = lme4::lmer(OG.DX ~  
                    S100A1.positive + 
                    S100A1.Staining + 
                    (1 | Scorer), data = df)
lmm3
stargazer::stargazer(lmm3, type = "text",
                     digits = 3,
                     star.cutoffs = c(0.05, 0.01, 0.001),
                     digit.separator = "")
library(sjplot)
sjp.lmer(lmm2)

### chi square test
df$Scorer = as.factor(df$Scorer)
chisq.test(df$OG.DX, df[ , c(3,6,8)])

df %>% group_by(Scorer) %>% select(2, 3, 8)
df %>% group_by(Scorer) %>% 
  select(2,3,8) %>% table()

chisq.test(df$H.E[df$Scorer == 1], df$H.E[df$Scorer == 2])

chisq.test(df$OG.DX[df$Scorer == 1], df$H.E[df$Scorer == 3])
library(gplots)
balloonplot(df[ , 2] , df[ , c(3, 6:7)], df[ , 8])

t.test(df[ , 2], df[ , 3])
t.test(df$H.E[df$Scorer == 1], df$H.E[df$Scorer == 2])

t.test(df$OG.DX[df$Scorer == 1], df$H.E[df$Scorer == 3])
t.test(df$H.E[df$Scorer == 1], df$H.E[df$Scorer == 3])

t.test(df$S100A1.positive[df$Scorer == 1], df$S100A1.positive[df$Scorer == 3])

df$S100A1.Staining = df$S100A1.Staining %>% replace_na(0)
df %>% group_by(Scorer) %>% 
  select(S100A1.Staining) %>% table()
