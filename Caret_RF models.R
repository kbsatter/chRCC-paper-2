
rm(list = ls())
library(tidyverse)
library(caret)
library(magrittr)
# library()
set.seed(142)
pheno = read.csv("/Users/FBINSATTER/Desktop/chRCC/PhenoData.csv")
pheno %<>% filter(Histology2 != "N")
X_test = read.csv("./RCCs/Pre-BatchCorrectiondata.csv", row.names = 1)
X_test %<>% t %>% as.data.frame() %>% select(1:30) %>% log2 %>% 
  scale(center = T, scale = T)

X_train = read.csv("trainData.csv", row.names = 1)
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
X_train %<>% select(-1) %>% scale(center = T, scale = T)
Histology = ifelse(pheno$Histology2 == "chRCC", 1,2) %>% as.factor()
X_train = cbind.data.frame(Histology, X_train)

Histology2 = c(rep(1, 16), rep(2, 14)) %>% as.factor()
X_test = cbind.data.frame(Histology2, X_test)

colnames(X_test)[7] = "WDR63"

### check with UMAP
library(umap)
u1 = umap(X_train[ , -1])
u1$layout %>% as.data.frame() %>% 
  ggplot(aes(V1, V2, color = pheno$Histology2)) +
  geom_point() + theme_bw() + xlab("") + ylab("") +
  ggtitle("chRCC-RO Discovery Data (Scaled)")

u2 = umap(X_test[ , -1],
          n_neighbors = 5,
          mis_dist = 0.001)
u2$layout %>% as.data.frame() %>% 
  ggplot(aes(V1, V2, color = Histology2)) +
  geom_point() + theme_bw() + xlab("") + ylab("") +
  ggtitle("chRCC-RO NanoString Data") + labs(color = "Histology")

### rf
rf.defaults = train(Histology ~ .,
                    data = X_train,
                    method = "rf",
                    metric = "Accuracy",
                    tunegrid = tunegrid,
                    trcontrol = control)
plot(rf.defaults$finalModel)
plot(rf.defaults)
rf.predict = rf.defaults %>% predict(X_test)
confusionMatrix(Histology2, rf.predict)



### svm
svm.default = train(Histology ~ .,
                    data = X_train,
                    method = "svmLinear",
                    metric = "Accuracy",
                    tunegrid = tunegrid,
                    trcontrol = control)

svm.default
svm.predict = svm.default %>% predict(X_test)
confusionMatrix(Histology2, svm.predict)

test_score = cbind(rownames(X_test), Histology2, rf.predict, svm.predict)
write.csv(test_score, file = "Supervised_Models.csv")



####

library(lattice)
library(AppliedPredictiveModeling)
transparentTheme(trans = .9)

featurePlot(x = X_train[ , c("KRT7", "S100A1", "AQP6",
                            "HOOK2")], 
            y = X_train$Histology,
            plot = "density", 
            scales = list(x = list(relation="free"), 
                          y = list(relation="free")), 
            adjust = 1.5, 
            pch = "|", 
            layout = c(4, 1), 
            auto.key = list(columns = 3))

featurePlot(x = X_test[ , c("KRT7", "S100A1", "AQP6",
                             "HOOK2")], 
            y = X_test$Histology2,
            plot = "density", 
            scales = list(x = list(relation="free"), 
                          y = list(relation="free")), 
            adjust = 1.5, 
            pch = "|", 
            layout = c(4, 1), 
            auto.key = list(columns = 3))

### base SVM
library(tidyverse)    # data manipulation and visualization
library(kernlab)      # SVM methodology
library(e1071)        # SVM methodology
library(RColorBrewer) # customized coloring of plots

svmfit = svm(X_train[ , 1] ~ ., data = X_train[, -1], kernel = "linear", scale = F)
plot(svmfit, X_train, BSPRY~KRT7)
