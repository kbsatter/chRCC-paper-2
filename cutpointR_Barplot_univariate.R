
rm(list = ls())
### cutploint R

library(tidyverse)
library(magrittr)
library(cutpointr)

df = read.csv("./RCCs/Pre-BatchCorrectiondata.csv", row.names = 1)
df2 = df[1:30, ] %>% t %>% log2 %>% round(2) %>% as.data.frame()

Group = c(rep(1, 16), rep(2, 14))
df2 = cbind.data.frame(Group, df2)

genes = colnames(df2[-1])

cp = c()
# i = 3
for (i in seq_along(genes)) {
  dt3 = df2 %>% dplyr::select(c(Group, i+1))
  names(dt3)[2] = "gid" 
  cp1 = cutpointr(dt3, gid, Group, method = maximize_metric, metric = sum_sens_spec)
  cp2 = cbind(genes[i], cp1[[2]],cp1[[3]], cp1[[4]],cp1[[5]],
              cp1[[6]], cp1[[7]],cp1[[8]],cp1[[9]], cp1[[10]],
              cp1[[11]],cp1[[12]], cp1[[13]])
  cp = rbind(cp2, cp)
  print(i)
}
dev.off()
colnames(cp) = c("Gene", names(cp1[2]),names(cp1[3]), 
                 names(cp1[4]),names(cp1[5]),names(cp1[6]), names(cp1[7]),
                 names(cp1[8]),names(cp1[9]), names(cp1[10]),
                 names(cp1[11]),names(cp1[12]), names(cp1[13]))
write.csv(cp, "NanoString_Cutpoints.csv")

cp_data = read.csv("NanoString_Cutpoints.csv")
summary(cp_data)


AUC = cp_data %>% rownames_to_column() %>%
  ggplot(aes(x = rowname, y = AUC)) +
  geom_bar(stat = "identity", show.legend = F) + theme_classic() +
  theme(text = element_text(angle = 90)) + xlab("") 
Sensitivity = cp_data %>% rownames_to_column() %>%
  ggplot(aes(x = rowname, y = Sensitivity)) +
  geom_bar(stat = "identity", show.legend = F) + theme_classic() +
  theme(text = element_text(angle = 90)) + xlab("") 
Specificity = cp_data %>% rownames_to_column() %>%
  ggplot(aes(x = rowname, y = Specificity)) +
  geom_bar(stat = "identity", show.legend = F) + theme_classic() +
  theme(text = element_text(angle = 90)) + xlab("") 
Accuracy = cp_data %>% rownames_to_column() %>%
  ggplot(aes(x = rowname, y = Accuracy)) +
  geom_bar(stat = "identity", show.legend = F) + theme_classic() +
  theme(text = element_text(angle = 90)) + xlab("") 

grid.arrange(AUC, Sensitivity, Specificity, Accuracy)


###
cp_data$Gene = as.factor(cp_data$Gene)
cp_data %>% 
  ggplot(aes(x = Gene, group = 1)) +
  geom_path(aes(y = Sensitivity, color = "red")) +
  geom_path(aes(y = Specificity, color = "blue")) +
  geom_path(aes(y = Accuracy, color = "green")) +
  geom_path(aes(y = AUC, color = "pink")) +
  theme_classic() +
  ylab("") +
  theme(axis.text.x = element_text(angle = 90))

cp_data2 = cp_data %>% 
  select(Gene, Sensitivity, Specificity, Accuracy, AUC) %>%
  gather(key = "variable", value = "value", -Gene)

cp_data2 %>% filter(variable == "Sensitivity" | variable == "Specificity") %>% 
  arrange(variable) %>% 
  ggplot(aes(x = Gene, y = value)) +
  geom_path(aes(color = variable, group = variable, linetype = variable)) +
  geom_point(aes(color = variable, group = 1)) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
  

cp_data2 %>% filter(variable == "Sensitivity" | variable == "Specificity") %>%
  arrange(variable) %>% 
  ggplot(aes(x = Gene, y = value)) +
  geom_path(aes(color = variable, group = 1)) +
  geom_point(aes(color = variable, group = 1)) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90))



## for sensitivity and specificity
plot(cp_data$Sensitivity, ylim=c(0.3,1.2),type='o', axes = F, xlab='', ylab = '',,lwd=2)
points(cp_data$Specificity, type='o', lty=3, lwd=2)
axis(1, at=1:30, labels = cp_data$Gene, las=2, font=2, cex=0.7)
yval=seq(0.3,1.2, length.out=4)
axis(2, at=yval, yval, font=2, las=2,cex=0.7)
legend("top", horiz=T, legend=c("Sensitivity","Specificity"), lty=c(1,3), lwd=2, bty="n")

## for AUC and accuracy
plot(cp_data$AUC, ylim=c(0.3,1.2),type='o', axes = F, xlab='', ylab = '',,lwd=2)
points(cp_data$Accuracy, type='o', lty=3, lwd=2)
axis(1, at=1:30, labels = cp_data$Gene, las=2, font=2, cex=0.7)
yval=seq(0.3,1.2, length.out=4)
axis(2, at=yval, yval, font=2, las=2,cex=0.7)
legend("top", horiz=T, legend=c("AUC","Accuracy"), lty=c(1,3), lwd=2, bty="n")



## for AUC and accuracy
par(oma=c(2,1,0.5,0))
temp= cp_data[,c("Gene", "AUC","Accuracy")]
xx<-barplot(t(temp[,2:3]), ylim=c(0,1.2),axes = F, xlab='', ylab = '',beside=T)
xvals=apply(xx, 2, mean)
axis(1, at=xvals, labels = cp_data$Gene, las=2, font=2, cex=0.7)
yval=seq(0,1.2, length.out=6)
axis(2, at=yval, yval, font=2, las=2,cex=0.7)
legend("top", horiz=T, inset=c(0,-0.05),legend=c("AUC","Accuracy"), 
       pch=c(15,22), col=c('black', "black"), bty="n", cex = 1)


par(oma=c(2,1,0.5,0))
temp2= cp_data[,c("Gene", "Sensitivity","Specificity")]
xx<-barplot(t(temp2[,2:3]), ylim=c(0,1.2),axes = F, xlab='', ylab = '',beside=T)
xvals=apply(xx, 2, mean)
axis(1, at=xvals, labels = cp_data$Gene, las=2, font=2, cex=0.7)
yval=seq(0,1.2, length.out=6)
axis(2, at=yval, yval, font=2, las=2,cex=0.7)
legend("top", horiz=T, inset=c(0,-0.05),
       legend=c("Sensitivity","Specificity"), pch=c(15,22), 
       col=c('black', "black"), bty="n", cex = 1)

