#setwd("/Users/mishrap/Documents/GitHub")
library(dplyr)
x=read.csv("TumorCompendium_v11_PolyA_hugo_log2tpm_58581genes_2020-04-09.tsv", header=TRUE, sep="\t", check.names=FALSE, row.names=1)
#Subset columns 
#th_comp <- x %>% dplyr:: select(starts_with("TH"))
#treehouse ids
#th=read.csv("treehouse_neuro_ids", sep="\t", header=FALSE)
#th_t = t(th)



#target = read.table("TARGET_ids")

#Most neuroblastoma target ids start with TARGET-30
#This info comes from "target" file
df_target <- x %>% dplyr:: select(starts_with("TARGET-30"))

high <-read.csv("HighRisk_NBL_changes.csv",header=TRUE, check.names = FALSE)
DF4<- df
#replace the last part of the colnames in df 
df1<-data.frame()
for (i in 1:162){
  col = colnames(df_target)
  
  d =strsplit(col[i], split = "-")
  vecRow=data.frame()
  vecRow = paste(d[[1]][1],(paste(d[[1]][2],d[[1]][3], sep="-")), sep="-")
  df1 <- as.data.frame(rbind(df1,vecRow))
}
df2<-t(df1) 
rownames(df2)<- NULL

colnames(df_target)<-df2
targetl =t(target)
DF3 <- df_target[intersect(names(df_target),high$TARGET.USI)]

#DF3 <- df[intersect(names(df),targetl)]
#DF3 <- df[intersect(names(df),targetl)]
#DF4 missing
#treehouse data
DF5 <- th_comp[intersect(names(th_comp),th_t)]

write.csv(DF3,file="u", quote=FALSE)
library(ggplot2)

ptprc <- DF3["PTPRC",]
ptprc_t <- t(ptprc)
ptprc_t_t <- as.data.frame(t(ptprc_t))
ptprc_th_t<-as.data.frame(t(ptprc_th))


gpc2 <- DF3["GPC2",]
gpc2_t <- as.data.frame(2^t(gpc2))
colnames(gpc2_t) =c("TPM")
gpc2_t <- gpc2_t %>% add_count()

st8 <- DF3["ST8SIA1",]
st8_t <- as.data.frame(2^t(st8))
colnames(st8_t) =c("TPM")
st8_t <- st8_t %>% add_count()
rownames(st8_t) <- NULL
colnames(st8_t ) =c("value", "total_count")
st8_t$variable <- "Neuroblastoma"

b4<- DF3["B4GALNT1",]
b4_t <- as.data.frame(2^t(b4))

colnames(b4_t) =c("TPM")
b4_t <- b4_t %>% add_count()
rownames(b4_t) <- NULL
colnames(b4_t ) =c("value", "total_count")
b4_t$variable <- "Neuroblastoma"

gpc2 <- DF4["GPC2",]
gpc2_th <- DF5["GPC2",]
ptprc <- DF4["PTPRC",]
st8 <- DF4["ST8SIA1",]
st8_th <- DF5["ST8SIA1",]
st8s_th_t = as.data.frame(t(st8_th))
b4<- DF4["B4GALNT1",]
b4_th<- DF5["B4GALNT1",]
b4_th_t<- as.data.frame(t(b4_th))
nkg7_th<-DF5["B4GALNT1",]
nkg7_t<- DF4["B4GALNT1",]
ptprc_th <- DF5["PTPRC",]

#
ptprct = as.data.frame(t(ptprc))
ptprc_th_t = as.data.frame(t(ptprc_th))

gpc2t = as.data.frame(t(gpc2))
gpc2th_t = as.data.frame(t(gpc2_th))
st8s = as.data.frame(t(st8))
nkg7_th_t = as.data.frame(t(nkg7_th))
nkg7_t_t = as.data.frame(t(nkg7_t))

b4gal = as.data.frame(t(b4))
colnames(ptprct) <- "PTPRC"
colnames(gpc2t) <- "GPC2"
colnames(st8s) <- "ST8SIA1"
colnames(b4gal) <- "B4GALNT1"
colnames(Neuro_nkg7)[2] <- "value"
ticksx <- theme(
  axis.text.x = element_blank(),
  axis.ticks.x  = element_blank())

All_gpc2_neuro <- rbind(gpc2th_t, gpc2t)
All_b4_neuro <- rbind(b4_th_t, b4gal)
All_st8_neuro <- rbind(st8s_th_t, st8s)
All_nkg7_neuro <- rbind(nkg7_th_t, nkg7_t_t)
All_ptprc_neuro <- rbind(ptprc_th_t, ptprct)


All_gpc2_neuro$variable ="NeuroBlastoma"
All_b4_neuro$variable = "Neuroblastoma"
All_st8_neuro$variable = "Neuroblastoma"
All_nkg7_neuro$variable = "Neuroblastoma"
All_ptprc_neuro$variable = "Neuroblastoma"
Neuro_gpc<-All_gpc2_neuro[,c(2,1)]
Neuro_b4<-All_b4_neuro[,c(2,1)]
Neuro_st8<-All_st8_neuro[,c(2,1)]
Neuro_nkg7<-All_nkg7_neuro[,c(2,1)]
Neuro_ptprc<-All_ptprc_neuro[,c(2,1)]
colnames(Neuro_ptprc)[2] ="value"


rownames(Neuro_gpc) <-NULL
rownames(Neuro_b4) <-NULL
rownames(Neuro_st8) <-NULL
rownames(Neuro_nkg7) <-NULL
Neuro_gpc <- Neuro_gpc %>% group_by(variable) %>% mutate(total_count = n())
Neuro_b4 <- Neuro_b4 %>% group_by(variable) %>% mutate(total_count = n())
Neuro_st8 <- Neuro_st8 %>% group_by(variable) %>% mutate(total_count = n())
Neuro_nkg7 <- Neuro_nkg7 %>% group_by(variable) %>% mutate(total_count = n())
Neuro_ptprc <- Neuro_ptprc %>% group_by(variable) %>% mutate(total_count = n())

pdf("Treehouse_neuroblastoma.pdf")
gpc2<-ggplot(gpc2t, aes(x=0, y=GPC2)) + 
  geom_boxplot(fill="lightskyblue3") + xlab("GPC2") +ylab("NBL log2 TPM") + theme_classic() + ticksx 

ptp<-ggplot(ptprct, aes( y=PTPRC)) + 
  geom_boxplot(fill="lightskyblue3") + xlab("PTPRC") +ylab("NBL log2 TPM") + theme_classic() + ticksx 
st <- ggplot(st8s, aes(y=ST8SIA1)) + 
  geom_boxplot(fill="lightskyblue3") +ylab("NBL log2 TPM") + xlab("ST8SIA1") + theme_classic() + ticksx 

b4g <- ggplot(b4gal, aes( y=B4GALNT1)) + 
  geom_boxplot(fill="lightskyblue3") + ylab("NBL log2 TPM") + xlab("B4GALNT1") + theme_classic() + ticksx 

plot_grid( gpc2, ptp, st, b4g)

dev.off()

library(data.table)
setnames(df, df$old, dat$new)
names(mt)
colnames(df)<-c(df1)

str_replace_all(colnames(df), "")
x=read.table("tumor_head",sep="\t", header=FALSE)

colnames(target) =c("Target_ids")
tumor_head=read.table("tumor_head", sep="\t", header=FALSE)
colnames(tumor_head) =c("Target_ids")

df1<-data.frame()
#df1 is target ids
for (i in 1:nrow(tumor_head)){
  
  d =strsplit(tumor_head[i,], split = "-")
  vecRow=data.frame()
  vecRow = paste(d[[1]][1],(paste(d[[1]][2],d[[1]][3], sep="-")), sep="-")

  df1 <- as.data.frame(rbind(df1,vecRow))

}
colnames(df1)="Target_ids"

inner_join(target, df1, by=c("Target_ids"))