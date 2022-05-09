


rm(list = ls())
library(tidyverse)
library(umap)
library(magrittr)
library(ComplexHeatmap)
library(ggalluvial)
library(caret)
library(e1071)
library(caTools)
library(tidymodels)
library(embed)
df = read.csv("./RCCs/Pre-BatchCorrectiondata.csv", row.names = 1)

df2 = t(df[ 1:30, ])

set.seed(142)

u1 = df2 %>% scale(center = T, scale = T) %>%
  umap(n_neighbors = 10,
          min_dist = 0.001)

colnames(df2)
rownames(df2)[1:16]
Group = c(rep(1, 16), rep(2, 14)) %>% as.factor()

u1$layout %>% as_tibble() %>% 
  ggplot(aes(V1,V2)) +
  geom_point(aes(color = Group)) + theme_bw() +
  theme(text = element_text(size = 12))



auc = colAUC(df2, Group, plotROC = F, alg = "ROC")
gene_auc = t(auc)
gene_auc = data.frame(rownames(gene_auc), gene_auc[ , 1])
colnames(gene_auc) = c("Gene", "AUC")  
gene_auc2 = gene_auc %>% top_n(5) %>% arrange(-AUC)

auc2 = df2 %>% as.data.frame() %>% 
  select(all_of(gene_auc2$Gene)) %>% 
  colAUC(Group, plotROC = T, alg = "ROC")
gene_auc %>% 
  ggplot(aes(x = Gene, y = AUC)) +
  geom_bar(stat = "identity", show.legend = F) + theme_classic() +
  theme(text = element_text(angle = 90)) + xlab("") 

Histology = pheno$Histology2
u3 = umap(testdata2, n_neighbors = 50,
          min_dist = 0.01)
u3$layout %>% as.data.frame() %>% 
  ggplot(aes(V1,V2, color = Histology)) +
  geom_point() +
  theme(text = element_text(size = 12))
 
### tidymodels
library(tidymodels)
library(embed)
## get original data
testdata = read.csv("/Users/FBINSATTER/Desktop/chRCC/TestData1.csv", row.names = 1)
COGS = read.csv("/Users/FBINSATTER/Desktop/chRCC/GS30_allinformation.csv")

testdata2 = testdata[ match(COGS$Gene, rownames(testdata)),  ]
pheno = read.csv("/Users/FBINSATTER/Desktop/chRCC/PhenoData.csv")
pheno %<>% filter(Histology2 != "N")
testdata2 = testdata2[ , match(pheno$geo_accession, colnames(testdata2))]
testdata2 %<>% t %>% as.data.frame()
match(rownames(testdata2), pheno$geo_accession)
geneID = COGS[match(colnames(testdata2), COGS$Gene), "SYMBOL"]
colnames(testdata2) = geneID
Histology2 = pheno$Histology2
train_data = cbind.data.frame(Histology2, testdata2)

supervised <- 
  recipe(Histology2 ~ ., data = train_data) %>%
  step_center(all_predictors()) %>% 
  step_scale(all_predictors()) %>% 
  step_umap(all_predictors(), outcome = vars(Histology2), num_comp = 2) %>% 
  prep(training = train_data)

colnames(df2)[6] = "WDR63"
theme_set(theme_bw())
Group = ifelse(Group == 1, "chRCC", "RO") %>% as.factor()
bake(supervised, new_data = df2, starts_with("umap")) %>% 
  ggplot(aes(x = UMAP1, y = UMAP2, col = Group)) + 
  geom_point(alpha = 1, size = 2) +
  theme(text = element_text(size = 12))


rm(testdata)
# save.image("Supervisedmodel.Rdata")

gene_auc %>% arrange(AUC)
gene_auc %>% arrange(-AUC)


Group = c(rep("chRCC", 16), rep("RO", 14))
