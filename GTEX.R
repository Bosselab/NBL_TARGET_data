library(dplyr)
library(ggplot2)
library(stringr)
setwd("GTEX")
#Name all the tissues 
tissues=c("adipose","adrenal","artery","bladder","brain","breast","cells","cervix","colon","esophagus","fallopian","heart","kidney","liver","lung","muscle","nerve","ovary","pancreas","pituitary","prostate","skin",
          "small_intestine","spleen","stomach","testis","thyroid","uterus","vagina","whole_blood" )

#Find all the files in the directory
files=list.files(pattern="gene_tpm_2017")


df_file1= read.table(files1, sep="\t", header=TRUE, skip=2, check.names=FALSE)
df_file2= read.table(files2, sep="\t", header=TRUE, skip=2, check.names=FALSE)
#intersect(colnames(df_file1), colnames(df_file2))
GPC2=data.frame()
ST8SIA1=data.frame()
files_gpc2=list.files(pattern="GPC2")
files_b4=list.files(pattern="ST8SIA1")



#Extract the genes from each GTEX tissue and write it to files
for (i in files){
  
  #Just remove the gct and keep the tissue name
  file_name = gsub(".gct", "", (gsub("gene_tpm_2017-06-05_v8_", "", i)))
  print (i)
  df = read.table(i, sep="\t", header=TRUE, skip=2, check.names = FALSE)
  genes = df[df$Description == "GPC2",]
  genes_t = t(genes)
  genes_t = genes_t[-c(1:3),]
  genes_t <- as.numeric(genes_t)
  genes_t <- as.data.frame(genes_t)
  colnames(genes_t) = file_name
  file_name_csv = paste("GPC2",file_name, sep="_")
  write.csv(genes_t, file=file_name_csv,quote=FALSE, row.names=FALSE )
    }
files_st8=list.files(pattern="GPC2_")
gpc2_files=data.frame()
ST8SIA1_files=data.frame()
#Select the files that have this 
for (t in tissues){
  st8= files_st8[grep( t, files, fixed = TRUE)]
  print (st8)
  ST8SIA1_files = rbind(ST8SIA1_files, st8)
  
}

library(reshape2)
#All GPC2 files were pasted
X=read.table("All_GPC2_files", header=TRUE, sep="\t")

X=read.table("all_GPC2", header=TRUE, sep="\t")
X=read.table("all_B4", header=TRUE, sep="\t")
X=read.table("all_ST8", header=TRUE, sep="\t")
X=read.table("all_NKG7", header=TRUE, sep="\t")
X=read.table("all_PTPRC", header=TRUE, sep="\t")

#For each gene
all<-melt(X)
all_b4<-melt(X)
head(X)
all<- all_b4
all$variable <- str_split_fixed(all$variable, pattern="_", n=2)[,1]
all <- na.omit(all)
#all_with_count <- all %>% group_by(variable)  %>% add_count()

all_with_count <- all %>% group_by(variable)  %>% summarise(total_count=n())
all_counts <- left_join(all,all_with_count, by="variable")
all_counts$variable = str_replace(all_counts$variable, "^\\w{1}", toupper)
colnames(Neuro_st8)[2] <- "value"
rownames(gpc2_t)<- NULL
colnames(gpc2_t)[1] = "value"
colnames(gpc2_t)[2] = "total_count"
gpc2_t[3] <- "Neuroblastoma"
colnames(gpc2_t)[3] = "variable"
gpc2_t$variable = paste(gpc2_t$variable, "(", gpc2_t$total_count, ")")
all_counts$variable = paste(all_counts$variable, "(", all_counts$total_count, ")")



Neuro_gpc$value <- 2^Neuro_gpc$value
Neuro_nkg7$value <- 2^Neuro_nkg7$value
Neuro_ptprc$value <- 2^Neuro_ptprc$value

Neuro_b4$value <- 2^Neuro_b4$value
Neuro_b4$variable = paste(Neuro_b4$variable, "(", Neuro_b4$total_count, ")")
Neuro_st8$variable = paste(Neuro_st8$variable, "(", Neuro_st8$total_count, ")")
#filter if needed
Neuro_b4<-Neuro_b4[Neuro_b4$value < 900,]
Neuro_nkg7$variable = paste(Neuro_nkg7$variable, "(", Neuro_nkg7$total_count, ")")

Neuro_ptprc$variable = paste(Neuro_ptprc$variable, "(", Neuro_ptprc$total_count, ")")


b<-rbind(Neuro_gpc, all_counts)
g <-rbind(Neuro_b4, all_counts)
b<-rbind(Neuro_ptprc, all_counts)


b<- rbind(gpc2_t, all_counts)
b1$variable<-gsub("Blastoma", "blastoma", b1$variable)
b$variable <- factor(b$variable, levels = unique(b$variable))


g <-rbind(Neuro_nkg7, all_counts)
b<-g
ggplot(all_counts, aes(x=variable, y=value)) + geom_boxplot(outlier.size=0, fill="pink3", color="gray1") + 
 theme_classic() + theme(axis.text.x =element_text(angle=90,hjust=0.95,vjust=0.2), axis.ticks.x = element_blank()) +
  xlab("") + ylab("TPM")

#remove cells:
b1<- b[-grep('Cells', b$variable), ]
b1$variable<-gsub("Whole", "Whole Blood", b1$variable)
b1$variable<-gsub("Small", "Small Intestine", b1$variable)
b1$variable<-gsub("Minor", "Minor Salivary gland", b1$variable)
b1$variable <- factor(b1$variable, levels = unique(b1$variable))

b1$type="Normal"
b1[grep("Neuro", b1$variable),]$type="High"

pdf("GPC2_TARGET_vs_GTEX_Jun25.pdf")
ggplot(  b1, aes(x=variable, y=value, fill=type, alpha=0.96)) + 
  geom_boxplot(outlier.size = 0, color="gray1") +
  scale_fill_manual(values=c("#12AF10", "gray68")) +

  theme_classic() + theme(axis.text.x =element_text(angle=90, hjust=0.95,vjust=0.2),legend.position = "None", axis.ticks.x = element_blank()) +
  xlab("") + ylab("TPM")
  
dev.off()
 
  

library(tidyverse)
dataset <- do.call("cbind",lapply(".",FUN=function(files){ read.csv(files)}))
list.files(path = ".",
           pattern="GPC2*", 
           full.names = T) %>% 
  map_df(~read_csv(.))

