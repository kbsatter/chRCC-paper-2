
rm(list = ls())
library(tidyverse)
library(magrittr)
library(caTools)
df = read.csv("./RCCs/Pre-BatchCorrectiondata.csv", row.names = 1)
df2 = t(df[ 1:30, ])
colnames(df2)
rownames(df2)[1:16]
Group = c(rep("chRCC", 16), rep("RO", 14)) %>% as.factor()

auc = colAUC(df2, Group, plotROC = F, alg = "ROC")
gene_auc = t(auc)
gene_auc = data.frame(rownames(gene_auc), gene_auc[ , 1])
colnames(gene_auc) = c("Gene", "AUC")  
gene_auc2 = gene_auc %>% top_n(5) %>% arrange(-AUC)

auc2 = df2 %>% as.data.frame() %>% 
  select(all_of(gene_auc2$Gene)) %>% 
  colAUC(Group, plotROC = T, alg = "ROC")
gene_auc2 %>% 
  ggplot(aes(x = Gene, y = AUC)) +
  geom_bar(stat = "identity", show.legend = F) + theme_classic() +
  theme(text = element_text(angle = 90)) + xlab("") 

df2 = log2(df2) %>% as.data.frame()
df2 %>%
  select(AQP6, BSPRY, ESRP1, HOOK2, SPINT2) %>% 
  add_column(Group) %>% 
  ggplot(aes(x = Group, y = AQP6, fill = Group)) +
  geom_boxplot()
df2 %>% add_column(Group) %>% 
  dplyr::select(Group, KRT7, S100A1, AQP6, BSPRY, HOOK2, LIMS1) %>%
  gather(Measure, Value, -Group) %>%
  ggplot(aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot() +
  facet_wrap(~Measure
             , scales = "free_y") +
  theme_classic() + theme(text = element_text(size = 15)) +
  scale_fill_manual(labels = c("chRCC", "RO"),
                     values = c("limegreen", "violet"))
df2 %>% add_column(Group) %>% 
  dplyr::select(Group, AQP6, BSPRY, ESRP1) %>%
  gather(Measure, Value, -Group) %>%
  ggplot(aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot() +
  facet_wrap(~Measure
             , scales = "free_y") +
  theme_classic() +
  scale_fill_manual(labels = c("chRCC", "RO"),
                    values = c("limegreen", "violet"))
###
pheno = read.csv("/Users/khaled/Downloads/chRCC/PhenoData.csv")
pheno %<>% filter(Histology2 != "N")
X_train = read.csv("./RCCs/COGS-Discovery Data.csv", row.names = 1)
identical(pheno$geo_accession, rownames(X_train))

X_train %>% select(AQP6) %>% 
  add_column(pheno$Histology2) %>% 
  ggplot(aes(x = pheno$Histology2, y = AQP6)) +
  geom_boxplot(aes(fill = pheno$Histology2))

Histology = pheno$Histology2
Histology[Histology == "Oncocytoma"] = "RO"
X_train = cbind.data.frame(Histology, X_train)
X_train %>% 
  dplyr::select(Histology, KRT7, S100A1, AQP6, BSPRY, HOOK2, LIMS1) %>%
  gather(Measure, Value, -Histology) %>%
  ggplot(aes(x = Histology, y = Value, fill = Histology)) +
  geom_boxplot() +
  facet_wrap(~Measure
             , scales = "free_y") +
  theme_classic() + theme(text = element_text(size = 15)) +
  scale_fill_manual(labels = c("chRCC", "RO"),
                     values = c("limegreen", "violet"))

df2 %>% add_column(Group) %>% 
  dplyr::select(Group, AP1M2:SUCLA2) %>%
  gather(Measure, Value, -Group) %>%
  ggplot(aes(x = Group, y = Value, color = Group)) +
  geom_violin() +
  facet_wrap(~Measure
             , scales = "free_y") +
  theme_classic() +
  scale_color_manual(labels = c(1, 2),
                     values = c("limegreen", "violet"))

