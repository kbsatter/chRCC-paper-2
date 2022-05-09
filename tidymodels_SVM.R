

rm(list = ls())
## Development of SVM, RF and GLMNET models
library(tidymodels)
library(tidyverse)
library(magrittr)

## rm(list = ls())
## load data
discovery_data = read_csv("trainData.csv")
glimpse(COGS)

GeneID = colnames(discovery_data)[-c(1,2)]
colnames(discovery_data)[1] = "sampleID"
pheno = read.csv("/Users/FBINSATTER/Desktop/chRCC/PhenoData.csv")
pheno %<>% filter(Histology2 != "N")
identical(discovery_data$sampleID, pheno$geo_accession)

Histology = pheno$Histology2 %>% as.factor()
exp_data = as.matrix(discovery_data[ , 3:30])
exp_data = scale(exp_data, scale = T, center = T)
X_train = cbind.data.frame(Histology, exp_data)

## test data
X_test = read.csv("RCCs/Pre-BatchCorrectiondata.csv", row.names = 1)
X_test %<>% t %>% as.data.frame() %>% select(1:30) %>% scale(center = T, scale = T)

Histology2 = c(rep("chRCC", 16), rep("Oncocytoma", 14)) %>% as.factor()
X_test = cbind.data.frame(Histology2, X_test)
### svm models
svm_linear_spec = svm_poly(degree = 1) %>% 
  set_mode("classification") %>% 
  set_engine("kernlab")

svm_linear_fit = svm_linear_spec %>% 
  set_args(cost = .2) %>% 
  fit(Histology ~ ., data = X_train)

svm_linear_fit

library(kernlab)
svm_linear_fit %>% ### multidimension can not plot
  extract_fit_engine() %>% 
  plot()

svm_linear_wf = workflow() %>%
  add_model(svm_linear_spec %>% set_args(cost = tune())) %>%
  add_formula(Histology ~ .)
set.seed(142)
sim_data_fold <- vfold_cv(X_train, strata = Histology)
param_grid <- grid_regular(cost(), levels = 10)

tune_res <- tune_grid(
  svm_linear_wf, 
  resamples = sim_data_fold, 
  grid = param_grid
)

autoplot(tune_res)

best_cost <- select_best(tune_res, metric = "accuracy")

svm_linear_final <- finalize_workflow(svm_linear_wf, best_cost)

svm_linear_fit <- svm_linear_final %>% fit(X_train)

augment(svm_linear_fit, new_data = X_test) %>%
  conf_mat(truth = Histology2, estimate = .pred_class)
