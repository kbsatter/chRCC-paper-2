
## Development of SVM, RF and GLMNET models
library(tidymodels)
library(tidyverse)
library(magrittr)

## rm(list = ls())
## load data
setwd("./RCCs")
discovery_data = read_csv("/Users/khaled/Downloads/chRCC/TestData1.csv")
COGS = read_csv("/Users/khaled/Downloads/chRCC/GS30_allinformation.csv")
glimpse(COGS)

discovery_data = discovery_data[ discovery_data$...1 %in% COGS$Gene, ]
GeneID = discovery_data$GENEID2
exp_data = as.matrix(discovery_data[ , 3:108])

pheno = read.csv("/Users/khaled/Downloads/chRCC/PhenoData.csv")
pheno %<>% filter(Histology2 != "N")
exp_data = exp_data[ , match(pheno$geo_accession, colnames(exp_data))]
exp_data %<>% t %>% as.data.frame()
colnames(exp_data) = GeneID


library(limma)
design = model.matrix(~ pheno$Histology2)
fit = lmFit(t(exp_data), design)
fit = eBayes(fit)
topgene = topTable(fit, n = 30)

topgene %>% mutate(FC = abs(logFC)) %>% select(FC)# %>% summary

SD.data = (apply(exp_data, 2, sd))
summary(abs(topgene$logFC)/SD.data)
### power analysis
library(pwr)
1.25/.6998653
pwr.t.test(n = NULL, d = 1.9, sig.level = 0.05, power = .80, type = "two.sample")
pwr.chisq.test(w = 1, N = NULL, sig.level = 0.05, power = 0.80, df = 30)

EnhancedVolcano:: EnhancedVolcano(topgene, rownames(topgene), x = "logFC", y = "adj.P.Val")
