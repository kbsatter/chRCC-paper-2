
rm(list=ls(all=TRUE))
library(pheatmap)
library(tidyverse)
library(limma)
library(gplots)
library(magrittr)
library(umap)

#geometric mean function
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

#read in and reshape dataframe
setwd("~/Downloads/ChRCC- paper 2/ChRCC-2/RCCs")
ncounter_data <- read.csv("All_sample2.csv", row.names = 1)
colnames(ncounter_data)
ncounter_data1 = as.matrix(ncounter_data[ , 3:32])

samp_ids <- colnames(ncounter_data1)
length(samp_ids)

samp_ids_split<-str_split(samp_ids,"_", simplify = T)

### Step 0: Visualize raw data to detect batch effects and outliers  #############
#make heatmap
my_palette <- colorRampPalette(c("red", "white", "blue"))(n = 30)
Cartridge = c(rep('Run 1', 11), rep("Run 2", 11), rep("Run 3", 8))
Histology = c(rep("chRCC", 16), rep("RO", 14))
mydf <- data.frame(row.names = samp_ids_split[ , 3], Cartridge = as.factor(Cartridge),
                   Histology = Histology)
colnames(ncounter_data1) = samp_ids_split[ , 3]

pheatmap(log2(ncounter_data1), color=rev(my_palette),cluster_cols = F, cluster_rows = F, 
         annotation_col = mydf,main = "Nanostring raw data",
         fontsize_row = 5,fontsize_col = 5)

pheatmap(log2(ncounter_data1[ 1:30, ]), color=rev(my_palette),cluster_cols = T, cluster_rows = T, 
         annotation_col = mydf,main = "Nanostring raw data",
         fontsize_row = 5,fontsize_col = 5)

"NanoString Raw data does not show the difference between the histological subtype"

rownames(ncounter_data)
COGS = rownames(ncounter_data)[1:30]
Housekeeping = rownames(ncounter_data)[31:33]
pos_con = rownames(ncounter_data)[42:47]
neg_con = rownames(ncounter_data)[34:41]

########      Step 1: Remove samples and genes with low counts  ################
myexpdat<-as.matrix(ncounter_data1[COGS,])
background_val <- max(ncounter_data1[grep("NEG_",rownames(ncounter_data1)),])

########      Step 2: background thresholding  ################

ncounter_data1 = ncounter_data1 - background_val
ncounter_data1 %<>% as.matrix()
ncounter_data1[which(ncounter_data1 <= 0)] = 1

pheatmap(log2(ncounter_data1), color=rev(my_palette),cluster_cols = F, cluster_rows = F, 
         annotation_col = mydf,main = "Nanostring raw data low count removed",
         fontsize_row = 8,fontsize_col = 8)

### sample COGS004,22, 23 did not have enough signal
"
Batch effect was not highlighted by background threshold.
"

########  Step 3: Normalize samples using Positive Control and Housekeeping Genes  ##########
##pos control norm
#visual inspection of pos control counts
pos_data = ncounter_data1[pos_con, ]
# pos_data = cbind(rownames(pos_data), pos_data)
# colnames(pos_data)[1] = "Probes"
data.frame(t(pos_data)) %>% 
  gather(Probe, Count, c("POS_A","POS_B","POS_C","POS_D","POS_E","POS_F")) -> 
  posdat_resh
posdat_resh %>% 
  ggplot(aes(x=Probe,y=log2(Count))) +
  geom_jitter() +
  coord_cartesian(ylim=c(0,20))
pos_data %>% log2 %>% t %>% boxplot()

# some qc work showed POS_E & POS_F wasn't detect for all samples, so will rm from normalization
myposdat<-as.matrix(ncounter_data1[pos_con[1:4],])
geomeans_persamp_prenorm = apply(myposdat,2,gm_mean)
NF_PC <- geomeans_persamp_prenorm/gm_mean(geomeans_persamp_prenorm)
posdat_norm <- sweep(myposdat, 2, NF_PC, "/")

data.frame(t(posdat_norm)) %>% 
  gather(Probe, Count, c("POS_A","POS_B","POS_C","POS_D")) -> 
  posdat_norm_resh
ggplot(posdat_norm_resh,aes(x=Probe,y=log2(Count))) +
  geom_jitter() +
  coord_cartesian(ylim=c(0,20)) + theme_classic()

posdat_norm %>% log2 %>% t %>% boxplot()
#Normalize data
ncounter_data1_pos_norm <- sweep(ncounter_data1, 2, NF_PC, "/")

## HK norm
myHKdat<-as.matrix(ncounter_data1_pos_norm[as.vector(Housekeeping),])

data.frame(t(myHKdat)) %>% gather(Probe, Count, as.vector(Housekeeping)) -> HKdat_resh
ggplot(HKdat_resh,aes(x=Probe,y=log2(Count)))+
  geom_jitter() +
  coord_cartesian(ylim=c(0,20)) + theme_classic()

#calculate normalization vectors
myHKdat_geomeans_persamp_prenorm<-apply(myHKdat,2,gm_mean) #myHKdat1 or myHKdat
NF_HK <- myHKdat_geomeans_persamp_prenorm/gm_mean(myHKdat_geomeans_persamp_prenorm)
myHKdat_norm <- sweep(myHKdat, 2, NF_HK, "/") #myHKdat1 or myHKdat


data.frame(t(myHKdat_norm)) %>% 
  gather(Probe, Count, as.vector(Housekeeping)) -> 
  myHKdat_norm_resh
ggplot(myHKdat_norm_resh,aes(x=Probe,y=log2(Count)))+
  geom_jitter() +
  coord_cartesian(ylim=c(0,20)) + theme_classic()

#HK norm data
ncounter_data2 <- sweep(ncounter_data1_pos_norm, 2, NF_HK, "/")

#check distribution of normalization factors
summary(NF_PC)
summary(NF_HK)

#Data after low count rm genes and samples, background thresholding, pos control norm, and HK gene norm
pheatmap(log2(ncounter_data2), color=rev(my_palette),
         cluster_cols = T, cluster_rows = T, 
         annotation_col = mydf, 
         main = "Nanostring norm data low count removed",
         fontsize_row = 8,fontsize_col = 8)
pheatmap(log2(ncounter_data2[ 1:30, ]), color=rev(my_palette),
         cluster_cols = F, cluster_rows = F, 
         annotation_col = mydf, 
         main = "Nanostring norm data low count removed",
         fontsize_row = 8,fontsize_col = 8)
pheatmap(log2(ncounter_data2[ 1:30, ]), color=rev(my_palette),
         cluster_cols = T, cluster_rows = T, 
         annotation_col = mydf, 
         main = "Nanostring norm POST Normalization",
         fontsize_row = 8,fontsize_col = 8)

"
Most of the steps listed above are standard and recommended by nanostring
"
umap_NSdata_before_BC = umap(t(ncounter_data2[1:30, ]))
umap_NSdata_before_BC$layout %>% 
  as.data.frame() %>% 
  ggplot(aes(V1,V2, color = Histology)) +
  geom_point() +theme_bw()
#################  Step 4: Visualize ####################
mydata<-log2(ncounter_data2[COGS, ]+1)
mydata_withHK<-log2(ncounter_data2[c(COGS,Housekeeping),]+1)

#PCA
#before
pca1 = prcomp(t(mydata), scale. = TRUE)
scores1 = as.data.frame(pca1$x)
ggplot(data = scores1, aes(x = PC1, y = PC2, label = rownames(scores1))) +
  geom_point(aes(color = mydf$Cartridge, size = 0.6)) +
  labs(color = "Cartridge") +theme_bw() +
  ggtitle("PCA plot of Nanostring Processed Before Batch effect removal")


#write dataframes for downstream analyses
date = gsub("-","",Sys.Date())
# write.csv(mydata,paste0("UEC_IGS_Nanostring_Data_PostQC_",date,".csv"))
# write.csv(mydata_withHK,paste0("UEC_IGS_Nanostring_Data_PostQC_",date,"_withHK.csv"))


#################  Step 5: Batch Effect Removal ####################

library(sva)
Histology = as.factor(Histology)
mydata_withHK_batch <- ComBat(dat = (mydata_withHK), batch = Cartridge, mod = Histology)
rownames(mydata_withHK_batch)

mydata_batch <- mydata_withHK_batch[rownames(mydata) %in% Housekeeping == F, ]
mydata_batch %<>% round(2)

u1 = umap::umap(t(mydata_withHK_batch[1:30, ]))
u1$layout %>% as_tibble() %>% 
  ggplot(aes(V1,V2, color = Histology)) +
  geom_point() + theme_classic()


#PCA
#after
pca2 = prcomp(t(mydata_batch), scale. = TRUE)
scores2 = as.data.frame(pca2$x)
ggplot(data = scores2, aes(x = PC1, y = PC2, label = rownames(scores2))) +
  geom_point(aes(color = as.factor(Cartridge))) + theme_bw() +
  ggtitle("PCA plot of Nanostring Processed After Batch effect removal")

#Data after low count rm genes and samples, background thresholding, pos control norm, and HK gene norm, along with batch effect removal
pheatmap(mydata_batch, color=rev(my_palette),
         cluster_cols = T, clustering_distance_cols = "euclidean", clustering_method = "ward.D2",
         cluster_rows = T, 
         annotation_col = mydf,main = "Nanostring norm data batch effect removed",
         fontsize_row = 8,fontsize_col = 8)

#scaled by gene
pheatmap(mydata_withHK_batch, color=rev(my_palette),
         cluster_cols = T, clustering_distance_cols = "euclidean", clustering_method = "ward.D2",
         cluster_rows = F, scale = "row",
         annotation_col = mydf,
         main = "Nanostring norm data batch effect removed",
         fontsize_row = 8,fontsize_col = 8)

pheatmap(mydata_withHK_batch[ 1:30, ], color=rev(my_palette),
         cluster_cols = T, clustering_distance_cols = "euclidean", clustering_method = "ward.D2",
         cluster_rows = T, scale = "row",
         annotation_col = mydf[ , 2],
         main = "Nanostring norm data batch effect removed",
         fontsize_row = 8,fontsize_col = 8)


library(ComplexHeatmap)
TA = HeatmapAnnotation(Histology = mydf$Histology)
mydata_withHK_batch %>% as.data.frame() %>% select(1:30) %>% as.matrix() %>% 
  scale(center = T, scale = T) %>%
  Heatmap(cluster_rows = T, cluster_columns = T,
          row_names_side = "left",
          col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
          row_names_gp = gpar(cex = 1),
          column_names_gp = gpar(cex = 0.5),
          row_dend_width = unit(2, "cm"),
          clustering_distance_rows = "maximum",
          clustering_method_rows = "ward.D",
          column_km = 2, top_annotation = TA, km = 3
          )
ncounter_data2
pheatmap(log2(ncounter_data2[1:33, ]), color=rev(my_palette),
         cluster_cols = T, clustering_distance_cols = "euclidean", clustering_method = "ward.D2",
         cluster_rows = F, scale = "row",
         annotation_col = mydf,
         main = "Nanostring norm data batch effect removed",
         fontsize_row = 8,fontsize_col = 8)

#write dataframes for downstream analyses
date<-gsub("-","",Sys.Date())
write.csv(mydata_batch, "ProcessedData.csv")


