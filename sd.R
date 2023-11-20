## Standar deviation for methylation betas
norm.beta<-read.table("/user/home/uy23281/scratch/model1/EPIC_mQTL/processed_data/methylation_data/norm.beta.gz",sep="\t",header=T) 
sds<-apply(norm.beta,1,function(x)sd(x,na.rm=T))
summary(sds)
sds<-as.data.frame(sds)
sds$cpg<-rownames(sds)
rownames(sds)<-NULL
colnames(sds)<-c("sd","cpg")
write.table(sds,"/user/home/uy23281/rmd.results/model1/sd.methy.txt",sep="\t",col.names=T,row.names=F,quote=F)

beta<-read.table("/home/uy23281/model4/bcell/EPIC_mQTL/processed_data/methylation_data/norm.beta.mqtl")
sds<-apply(beta,1,function(x)sd(x,na.rm=T))

write.table(sds,"/home/uy23281/model4/bcell/EPIC_mQTL/processed_data/methylation_data/sd.methy_bcell.txt",sep="\t",col.names=T,row.names=F,quote=F)


sd_neu<-read.table("/user/home/uy23281/scratch/rmd.results/model4/neu/sd.methy.txt",sep="\t",header=T)
sd_mono<-read.table("/user/home/uy23281/scratch/rmd.results/model4/mono/sd.methy.txt",sep="\t",header=T)
