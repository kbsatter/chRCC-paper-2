
rm(list = ls())

library(tidyverse)
library(caret)
library(magrittr)
# library()

pheno = read.csv("/Users/khaled/Downloads/chRCC/PhenoData.csv")
pheno %<>% filter(Histology2 != "N")
X_test = read.csv("./RCCs/Pre-BatchCorrectiondata.csv", row.names = 1)
X_test %<>% t %>% as.data.frame() %>% select(1:30) %>% 
  scale(center = T, scale = T)
# COGS = read_csv("/Users/khaled/Downloads/chRCC/GS30_allinformation.csv")

X_train = read.csv("./RCCs/COGS-Discovery Data.csv", row.names = 1)
identical(pheno$geo_accession, rownames(X_train))
## change DNAI3 for WDR63
### caret random forest
control = trainControl(method = "repeatedcv",
                       number = 10,
                       repeats = 3)

metric = "Accuracy"
set.seed(142)

mtry = sqrt(ncol(X_train))
tunegrid = expand.grid(.mtry = mtry)

### create data.frame
X_train %<>% scale(center = T, scale = T)
Histology = ifelse(pheno$Histology2 == "chRCC", 1,2) %>% as.factor()
X_train = cbind.data.frame(Histology, X_train)

Histology2 = c(rep(1, 16), rep(2, 14)) %>% as.factor()
X_test = cbind.data.frame(Histology2, X_test)

colnames(X_test)[7] = "WDR63"
# 
# ### check with UMAP
# library(umap)
# u1 = umap(X_train[ , -1])
# u1$layout %>% as.data.frame() %>% 
#   ggplot(aes(V1, V2, color = Histology)) +
#   geom_point() + theme_bw() + xlab("") + ylab("") +
#   ggtitle("chRCC-RO Discovery Data (Scaled)")
# 
# u2 = umap(X_test[ , -1])
# u2$layout %>% as.data.frame() %>% 
#   ggplot(aes(V1, V2, color = Histology2)) +
#   geom_point() + theme_bw() + xlab("") + ylab("") +
#   ggtitle("chRCC-RO NanoString Data")


### GLMNET
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

barplot(coefficients)
prediction = my_glmnet_model %>% predict(X_test)
confusionMatrix(Histology2, prediction)
summary(my_glmnet_model)
glm.res = X_test %>% select(1) %>% add_column(prediction)
write.csv(glm.res, "GLM-Result.csv")

final_score = cbind(rownames(X_test), X_test$Histology2, prediction)
write.csv(final_score, "2ndScore.csv")
