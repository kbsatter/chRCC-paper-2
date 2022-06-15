
rm(list = ls())

library(tidyverse)
library(caret)
library(magrittr)
# library()

pheno = read.csv("/Users/FBINSATTER/Desktop/chRCC/PhenoData.csv")
pheno %<>% filter(Histology2 != "N")
X_test = read.csv("./RCCs/Pre-BatchCorrectiondata.csv", row.names = 1)
X_test %<>% t %>% as.data.frame() %>% select(1:30) %>% 
  scale(center = T, scale = T)
# COGS = read_csv("/Users/khaled/Downloads/chRCC/GS30_allinformation.csv")

df = read.csv("./RCCs/COGS-Discovery Data.csv", row.names = 1)
identical(pheno$geo_accession, rownames(df))
Histology = ifelse(pheno$Histology2 == "chRCC", 1, 2) %>% as.factor()
df %<>% scale(center = T, scale = T)
df = cbind.data.frame(Histology, df)
index = createDataPartition(df$Histology, p = 0.8, list = F, times = 1)

X_train = df[ index, ]
X_val = df[ -index, ]

Histology2 = c(rep(1, 16), rep(2, 14)) %>% as.factor()
X_test = cbind.data.frame(Histology2, X_test)
## change DNAI3 for WDR63
colnames(X_test)[6] = "WDR63"
### caret random forest
control = trainControl(method = "repeatedcv",
                       number = 10,
                       repeats = 3)

metric = "Accuracy"
set.seed(142)
mtry = sqrt(ncol(X_train))
tunegrid = expand.grid(.mtry = mtry)
# 
### check with UMAP
# library(umap)
# u1 = umap(X_train[ , -1])
# u1$layout %>% as.data.frame() %>%
#   ggplot(aes(V1, V2, color = Histology)) +
#   geom_point() + theme_bw() + xlab("") + ylab("") +
#   ggtitle("chRCC-RO Discovery Data (Scaled)")
# # 
# u2 = umap(X_test[ , -1])
# u2$layout %>% as.data.frame() %>%
#   ggplot(aes(V1, V2, color = Histology2)) +
#   geom_point() + theme_bw() + xlab("") + ylab("") +
#   ggtitle("chRCC-RO NanoString Data")


### GLMNET
set.seed(142)
my_glmnet_model = train(
  Histology ~ .,
  method = "glmnet",
  data = X_train,
  tuneGrid = expand.grid(alpha = 0, lambda = 10^seq(-3,3,length = 100)),
  trControl = trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 5,
                           summaryFunction = defaultSummary))

plot(my_glmnet_model$finalModel, xvar = "dev")
coefficients <- coef(my_glmnet_model$finalModel, my_glmnet_model$bestTune$lambda)

prediction = my_glmnet_model %>% predict(X_test)
confusionMatrix(Histology2, prediction)
summary(my_glmnet_model)
glm.res = X_test %>% select(1) %>% add_column(prediction)
