library(seqinr)
library(ggplot2)
library(data.table)
library(ggrepel)
library(writexl)

# ref<-read.fasta("/media/nacho/Data/MTB_Tests/MTuberculosisH37Rv.fasta")
# files<-list.files("/media/nacho/Data/temp/TB894/TBTypResults/noise/", full.names = TRUE, pattern = ".tsv")
# typing2<-read.csv("/media/nacho/Data/MTB_Tests/Coll_scheme_classification.csv")
# typing<-read.csv("/media/nacho/Data/MTB_Tests/tbdb.barcode.bed",sep = "\t", header = FALSE)

 ref<-read.fasta("/home/docker/CommonFiles/Refs/MTuberculosisH37Rv.fasta")
 files<-list.files("/Data/TBTypResults", full.names = TRUE, pattern = ".tsv")
 typing2<-read.csv("/home/docker/CommonFiles/Schemes/Coll_scheme_classification.csv")
 typing<-read.csv("/home/docker/CommonFiles/Schemes/tbdb.barcode.bed",sep = "\t", header = FALSE)

colnames(typing)<-c("ID","PositionS","Position","lineage","Allele.change","SuperLineage","Gene","SNP")


for (i in 1:length(files)) {
  
df<-fread(files[i])
ref<-as.data.frame(toupper(as.character(ref[[1]])) )
colnames(ref)<-"Reference"
ref$Position<-c(1:nrow(ref))

df<-as.data.frame(df)
df<-df[,c(1:5)]

colnames(df)<-c("Position", "Noise","Reads", "Base1", "Base2")

dummy<-as.data.frame(c(1:nrow(ref)))
colnames(dummy)<-"Position"
dummy$Noise<-NA
dummy$Reads<-0
dummy$Base1<-"N"
dummy$Base2<-"N"
dummy<-dummy[-which(dummy$Position %in% df$Position),]
df<-rbind(df,dummy)

df<-df[order(df$Position), ]

df_typ<-df[which(df$Position %in% typing$Position),]

df_typ<-merge(df_typ, typing, by="Position")
df_typ<-merge(df_typ, ref, by="Position", all.y = FALSE)
df_typ$Index<-c(1:nrow(df_typ))
df_typ$Allele.change<-gsub("./","",df_typ$Allele.change)
df_typ$assignment<-NA
df_typ$assignment[which(df_typ$Base1==df_typ$Allele.change)]<-df_typ$lineage[which(df_typ$Base1==df_typ$Allele.change)]
df_typ$assignment[which(df_typ$Base1==df_typ$Reference)]<-"Reference"
df_typ$assignment[which(is.na(df_typ$assignment))]<-"Unknown"

df_typ$assignment2<-NA
df_typ$assignment2[which(df_typ$Base2==df_typ$Allele.change)]<-df_typ$lineage[which(df_typ$Base2==df_typ$Allele.change)]
df_typ$assignment2[which(df_typ$Base2==df_typ$Reference)]<-"Reference"
df_typ$assignment2[which(is.na(df_typ$assignment2))]<-"Unknown"

df_typ$Position<-as.character(df_typ$Position)
df_typ$Position<-factor(df_typ$Position, levels=df_typ$Position[order(as.numeric(df_typ$Position))])
df_typ$ToAnalyze<-"NO"
df_typ$ToAnalyze[-which(df_typ$assignment %in% c("Reference", "Unknown") )] <-"YES"

df_typ$ToAnalyze2<-"NO"
if(length(which(df_typ$Noise>0.1))>0){
  df_typ$ToAnalyze2[-which(df_typ$assignment2 %in% c("Reference", "Unknown") )] <-"YES"  
  df_typ$ToAnalyze2[which(df_typ$ToAnalyze2 == "YES" & df_typ$Noise <0.10) ] <-"NO" 
  if(length(which(df_typ$ToAnalyze2=="YES" & df_typ$Reads <100 ))>0){
    df_typ$ToAnalyze2[which(df_typ$ToAnalyze2=="YES" & df_typ$Reads <100 ) ] <-"NO" 
  } 
} 




df_typ$ToAnalyzeG<-"NO"
df_typ$ToAnalyzeG[which(df_typ$ToAnalyze =="YES" | df_typ$ToAnalyze2=="YES")]<-"YES"


length(unique(c(df_typ$assignment[which(df_typ$ToAnalyzeG=="YES")], df_typ$assignment2[which(df_typ$ToAnalyzeG=="YES")])))-1 

ggplot(df_typ[which(df_typ$ToAnalyzeG=="YES"),])+
  #geom_line(aes(Position, 1-Noise),alpha=0.4)+
  geom_hline(yintercept = 0.9, col="red", size=0.2,alpha=0.5)+
  geom_point(aes(Position, 1-Noise, col=assignment),size=3,alpha=0.3)+
  geom_text_repel(aes(Position, 1-Noise, label=paste("Depth:",Reads)), max.overlaps = 50)+
  geom_point(aes(Position, Noise, col=assignment2),size=3,alpha=0.3)+
  scale_color_manual(values = c(rainbow(length(unique(c(df_typ$assignment[which(df_typ$ToAnalyzeG=="YES")], df_typ$assignment2[which(df_typ$ToAnalyzeG=="YES")])))-1),"grey"))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(0,1.1)+
  xlab("Informative Positions in Rreference")+
  ylab("Certainty (1-Noise)")+
  ggtitle(paste("Sample:",gsub(".*/","",gsub("_noise.tsv","",files[i])), " | Noise: ",round(mean(df_typ$Noise,na.rm = TRUE)*100,1),"%", " | Coverage:", 
                100-round((length(which(df$Reads==0))*100)/nrow(df),1),
                sep = ""))

ggsave(paste(gsub("_noise.tsv","",files[i]),"_FULL_SNPHits.pdf",sep = ""),width = 10, height = 5 )


ggplot(df_typ[-which(df_typ$assignment %in% c("Reference","Unknown")),])+
  #geom_line(aes(Position, 1-Noise),alpha=0.4)+
  geom_hline(yintercept = 0.9, col="red", size=0.2,alpha=0.5)+
  geom_point(aes(Position, 1-Noise, col=assignment),size=3,alpha=0.3)+
  geom_text_repel(aes(Position, 1-Noise, label=paste("Depth:",Reads)), max.overlaps = 50)+
  scale_color_manual(values = c(rainbow(length(unique(df_typ$assignment))-1),"grey"))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(0,1.1)+

  xlab("Informative Positions in Rreference")+
  ylab("Certainty (1-Noise)")+
  ggtitle(paste("Sample:",gsub(".*/","",gsub("_noise.tsv","",files[i])), " | Noise: ",round(mean(df_typ$Noise,na.rm = TRUE)*100,1),"%", " | Coverage:", 
                100-round((length(which(df$Reads==0))*100)/nrow(df),1),
                sep = ""))
  
ggsave(paste(gsub("_noise.tsv","",files[i]),"_SNPHits.pdf",sep = ""),width = 10, height = 5 )



#df2<-df_typ[-which(df_typ$assignment %in% c("Reference","Unknown")),]
df2<-df_typ[which(df_typ$ToAnalyzeG =="YES"),]
df2$Sample<-paste("Sample:",gsub(".*/","",gsub("_noise.tsv","",files[i])), " | Noise: ",round(mean(df_typ$Noise,na.rm = TRUE)*100,1),"%", " | Coverage:", 
                  100-round((length(which(df$Reads==0))*100)/nrow(df),1),
                  sep = "")
if(!exists("df2out")){
  df2out<-df2
}else{
  df2out<-rbind(df2,df2out)
}

}



samples<-as.data.frame(unique(df2out$Sample))
colnames(samples)<-"Sample"
samples$Coverage<-gsub(".*\\| Coverage:", "", samples$Sample)
samples$Noise<- as.numeric(gsub(" .*","", gsub(".*\\| Noise: ","",gsub("%","", samples$Sample))))/1000
samples$Lineage<-NA
samples$LineageProbability<-NA

for (i in 1:nrow(samples)) {
  lin<-df2out[which(df2out$Sample==samples$Sample[i]),]
  lin.unique<-as.data.frame(unique(lin$lineage))
  colnames(lin.unique)<-"Lineage"
  lin.unique$Support<-0
  lin.unique$Score<-0
  lin.unique$Final<-"NO"
  
  for (k in 1:nrow(lin.unique)) {
    lin.unique$Score[k]<-sum(lin$Reads[which(lin$lineage==lin.unique$Lineage[k])])
  }
  
  for (k in 1:nrow(lin.unique)) {
    Score.to.add<- lin.unique$Score[k]
    index.to.add<- grep(lin.unique$Lineage[k],lin.unique$Lineage)
    lin.to.add<-lin.unique$Lineage[grep(lin.unique$Lineage[k],lin.unique$Lineage)]
    #index.to.add<-index.to.add[-which(lin.to.add==lin.unique$Lineage[k])]
    if(length(index.to.add)>0) lin.unique$Support[index.to.add]<- lin.unique$Support[index.to.add] + Score.to.add
    if(length(grep(lin.unique$Lineage[k],lin.unique$Lineage))==1) lin.unique$Final[k]<-"YES"
  }
  
  lin.unique<-lin.unique[which(lin.unique$Final=="YES"),]
  lin.unique$TotalScore<-lin.unique$Support/sum(lin.unique$Support)  
  
  if(length(which(lin.unique$TotalScore<0.05))>0) lin.unique<-lin.unique[-which(lin.unique$TotalScore<0.05),]
  
  lin.unique$TotalScore<-lin.unique$Support/sum(lin.unique$Support)  
  samples$Lineage[i]<-paste(lin.unique$Lineage[which(lin.unique$Final=="YES")], collapse = "/")
  

  samples$LineageProbability[i] <-   paste(lin.unique$Lineage[which(lin.unique$Final=="YES")], " : ",
                                      round(lin.unique$TotalScore[which(lin.unique$Final=="YES")]*100,2),"%",sep="",collapse = " / ")
  samples$SamplePurity[i]<-mean((1-lin$Noise),na.rm=TRUE)
  
}


samples$Sample<-gsub(" .*","",gsub("Sample:","",samples$Sample))
write_xlsx(samples, "LineageAssignmentScore.xlsx" )

