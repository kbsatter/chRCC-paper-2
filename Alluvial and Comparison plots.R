

## alluvial for unsupevised data
rm(list = ls())
library(ggalluvial)
library(tidyverse)
library(magrittr)
fscore = read.csv("/Users/Khaled/Downloads/Thesis/Final Score.csv")


fscore %>% 
  ggplot( aes(axis1 = Histology, axis3 = k_df.cluster,
              axis4 = DBU.clusters, axis5 = nmf_score, axis2 = cc_kmeans,
              axis6 = best_class)) +
  scale_x_discrete(limits = c("Histology", "CC.Kmeans", "Kmeans",
                              "DBU", "NMF", "COnsensus.HC"))+
  xlab("") +
  geom_alluvium(aes(fill = Histology)) +
  scale_fill_manual(values = c("limegreen","violet")) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal() + theme(text = element_text(size = 12)) +
  theme(text = element_text(color = "black", size = 12)) +
  ggtitle("chRCC-Onc",
          "Unsupervised Model Comparison")


## alluvial for supevised data
s.score = readxl::read_xlsx("/Users/khaled/Downloads/Thesis/Final Clusters.xlsx",
                            sheet = 2)
s.score %<>% as.data.frame()
colnames(s.score)[2:6] = c("Histology", "SVM", "RF", "GLMnet", "Supervised.UMAP")
s.score %>% 
  ggplot( aes(axis1 = Histology, axis2 = Supervised.UMAP,
              axis3 = SVM, axis4 = RF, axis5 = GLMnet)) +
  scale_x_discrete(limits = c("Histology", "Supervised.UMAP", "SVM",
                              "RF", "GLMnet"))+
  xlab("") +
  geom_alluvium(aes(fill = Histology)) +
  scale_fill_manual(values = c("limegreen","violet")) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal() + theme(text = element_text(size = 12)) +
  theme(text = element_text(color = "black", size = 12)) +
  ggtitle("chRCC-Onc",
          "Supervised Model Comparison")


### plot for sensitivity, specificity

snplot = readxl::read_xlsx("/Users/khaled/Downloads/Thesis/Score Table.xlsx",
                           sheet = 3, col_names = T)
snplot %<>% as.data.frame()
rownames(snplot) = snplot$Metric
snplot %<>% select(-1) %>% t %>% as.data.frame()

Accuracy = snplot %>% rownames_to_column("Model") %>% 
  ggplot(aes(x = Model, y = Accuracy, fill = Model)) +
  geom_bar(stat = "identity", show.legend = F) + theme_bw() +
  xlab("")

Specificity = snplot %>% rownames_to_column("Model") %>% 
  ggplot(aes(x = Model, y = Specificity, fill = Model)) +
  geom_bar(stat = "identity", show.legend = F) + theme_bw() +
  xlab("")

Sensitivity = snplot %>% rownames_to_column("Model") %>% 
  ggplot(aes(x = Model, y = Sensitivity, fill = Model)) +
  geom_bar(stat = "identity", show.legend = F) + theme_bw() +
  xlab("")

Neg.Pred.Value = snplot %>% rownames_to_column("Model") %>% 
  ggplot(aes(x = Model, y = `Neg. Pred value`, fill = Model)) +
  geom_bar(stat = "identity", show.legend = F) + theme_bw() +
  xlab("")
ggpubr::ggarrange(Accuracy, Sensitivity, Specificity, Neg.Pred.Value)

