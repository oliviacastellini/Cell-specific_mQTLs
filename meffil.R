#Meffil to format and clean IDAt sorted cell type DNAm 
# Open meffil 
library(meffil)
options(mc.cores=10)
# Create sample sheet, set recursive to TRUE so many files can be read across different directories
samplesheet <- meffil.create.samplesheet("/home/uy23281/data/mono_methylation/", recursive=TRUE)

# Read PRECISEADs infofile with sex and diseases to complete sex column in datasheet 
info<-read.table("/home/uy23281/data/SortedCelltype.methylation.258.minInfo.txt",header=T,sep="\t")
info<-subset(info,select=c("OMICID","Gender"))
# Read DNAm samplesheet information PREC 
batch3<-read.table("/home/uy23281/data/mono_methylation/SampleSheet_EPIC_GBarturen_TCelulares_Batch3.csv",header=T,sep=",")
batch3$sentrix<-paste0(batch3$Sentrix_ID,"_",batch3$Sentrix_Position)
batch3<-subset(batch3,select=c("Sample_Name","sentrix"))  

batch4<-read.table("/home/uy23281/data/mono_methylation/SampleSheet_EPIC_GBarturen_TCelulares_Batch4.csv",header=T,sep=",")
batch4$sentrix<-paste0(batch4$Sentrix_ID,"_",batch4$Sentrix_Position)
batch4<-subset(batch4,select=c("Sample_Name","sentrix")) 

batch8<-read.table("/home/uy23281/data/bcell_methylation/Batch8_B/SampleSheet_EPIC_GBarturen_Bcell_Batch8.csv",header=T,sep=";")
batch8$sentrix<-paste0(batch8$Sentrix_ID,"_",batch8$Sentrix_Position)
batch8<-subset(batch8,select=c("Sample_Name","sentrix")) 

batch<-rbind(batch3,batch4)
batch$Sample_Name<-gsub("-M","",batch$Sample_Name)  

# Merge sex column with idat individuals identifiers 
info<-merge(info,batch,by.x="OMICID",by.y="Sample_Name",all.x=F,all.y=T)

#Merge sex column with samplesheet and replace Sample_name column by the identifiers
samplesheet<-merge(samplesheet,info,by.x="Sample_Name",by.y="sentrix")
samplesheet$Sex<-ifelse(samplesheet$Gender =="Female","F","M")
samplesheet$Sample_Name<-samplesheet$OMICID

#Remove extra columns in samplesheet 
samplesheet<-samplesheet[,-8]
samplesheet<-samplesheet[,-7]

#Background correction, dye bias correction, sex prediction and cell count estimates

meffil.list.cell.type.references()
qc.objects <- meffil.qc(samplesheet, cell.type.reference="blood gse35069 complete", verbose=TRUE)
save(qc.objects,file="/home/uy23281/data/bcell_methylation/meffil_output/qc.objects.Robj")
write.table(samplesheet,"/home/uy23281/data/mono_methylation/meffil_output/samplesheet.txt",col.names=T,sep="\t",quote=F)

##Chechk which file is having a reading error if "schedule core not delivered results" error when running meffil.qc is found 
fnames <- samplesheet$Basename[which(sapply(qc.objects, class) == 'try-error')]
x <-sapply(fnames, function(f) { tryCatch(meffil:::read.rg(f),error=function(e) print(f))  })
samplesheet[which(samplesheet$Basename=="/home/uy23281/data/bcell_methylation//Batch8_B/206113360011/Bcell/206113360011_R02C01"),]
samplesheet<-samplesheet[which(samplesheet$Basename!="/home/uy23281/data/bcell_methylation//Batch8_B/206113360011/Bcell/206113360011_R02C01"),]

# Gentypes for ID mismatches
#featureset <- qc.objects[[1]]$featureset
#writeLines(meffil.snp.names(featureset), con="/home/uy23281/data/tcell_methylation/snp-names.txt")

genotypes <- meffil.extract.genotypes("/home/uy23281/data/genotypes.meffil.raw")
genotypes <- genotypes[,match(samplesheet$Sample_Name, colnames(genotypes))]
colnames(genotypes)<-samplesheet$Sample_Name
rownames(genotypes)<-gsub("X","",rownames(genotypes))

#chr:pos need to be change to rsids for meffil to analyze it. Same file for all celltype DNAm
rsids<-read.table("/home/uy23281/data/rsids_chrpos.txt",header=T)
## change snp format to the genotypes snp format which is 1.pos.A.T
rsids$snp<-(gsub(":", ".", gsub("\\.", "", rsids$snp)))
genotypes<-merge(rsids,genotypes,by.x="snp",by.y="row.names",all=F)
rownames(genotypes)<-genotypes$name
genotypes<-genotypes[,-c(1,2)]

#Generates QC report 
      qc.parameters <- meffil.qc.parameters(
	beadnum.samples.threshold             = 0.1,
	detectionp.samples.threshold          = 0.1,
	detectionp.cpgs.threshold             = 0.1, 
	beadnum.cpgs.threshold                = 0.1,
	sex.outlier.sd                        = 5,
	snp.concordance.threshold             = 0.95,
	sample.genotype.concordance.threshold = 0.8
)

qc.summary <- meffil.qc.summary(
	qc.objects,
	parameters = qc.parameters,
	genotypes=genotypes
)

save(qc.summary, file="/home/uy23281/data/bcell_methylation/meffil_output/qcsummary.Robj")

#QC report 
meffil.qc.report(qc.summary, output.file="/home/uy23281/data/bcell_methylation/meffil_output/qc-report.html")

##Remove bad samples 
outlier <- qc.summary$bad.samples
table(outlier$issue)
index <- outlier$issue %in% c("Control probe (dye.bias)", 
                              "Methylated vs Unmethylated",
                              "X-Y ratio outlier",
                              "Low bead numbers",
                              "Detection p-value",
                              "Sex mismatch",
                              "Genotype mismatch",
                              "Control probe (bisulfite1)",
                              "Control probe (bisulfite2)")

outlier <- outlier[index,]
length(qc.objects)
qc.objects <- meffil.remove.samples(qc.objects, outlier$sample.name)
length(qc.objects)
save(qc.objects,file="/home/uy23281/data/bcell_methylation/meffil_output/qc.objects.clean.Robj")

#Rerun summary and report with good samples 
qc.summary <- meffil.qc.summary(qc.objects,parameters=qc.parameters,genotypes=genotypes)
save(qc.summary, file="/home/uy23281/data/bcell_methylation/meffil_output/qcsummary.clean.Robj")
meffil.qc.report(qc.summary, output.file="/home/uy23281/data/bcell_methylation/meffil_output/qc-report.clean.html")


#FUNCTIONAL NORMALIZATION
y <- meffil.plot.pc.fit(qc.objects)
ggsave(y$plot,filename="/home/uy23281/data/bcell_methylation/meffil_output/pc.fit.pdf",height=6,width=6)
pc <- 10 # based on the previous plot 

norm.objects <- meffil.normalize.quantiles(qc.objects, number.pcs=pc)
save(norm.objects,file=paste("/home/uy23281/data/bcell_methylation/meffil_output/norm.obj.pc",pc,".Robj",sep=""))

norm.beta <- meffil.normalize.samples(norm.objects, cpglist.remove=qc.summary$bad.cpgs$name)
save(norm.beta,file=paste("/home/uy23281/data/bcell_methylation/meffil_output/norm.beta.pc",pc,".Robj",sep=""))

#Normalization report 
str(norm.objects[[1]]$samplesheet)

#You change it by running a loop

for (i in 1:length(norm.objects)){
norm.objects[[i]]$samplesheet$Slide<-as.factor(norm.objects[[i]]$samplesheet$Slide)
norm.objects[[i]]$samplesheet$Sex<-as.factor(norm.objects[[i]]$samplesheet$Sex)
norm.objects[[i]]$samplesheet$sentrix_row<-as.factor(norm.objects[[i]]$samplesheet$sentrix_row)
norm.objects[[i]]$samplesheet$sentrix_col<-as.factor(norm.objects[[i]]$samplesheet$sentrix_col)
}

batch_var<-c("Slide", "sentrix_row","sentrix_col","Sex")
norm.parameters <- meffil.normalization.parameters(
	norm.objects,
	variables=batch_var,
	control.pcs=1:10,
	batch.pcs=1:10,
	batch.threshold=0.01
)

pcs <- meffil.methylation.pcs(norm.beta,probe.range=20000)
save(pcs,file="/home/uy23281/data/bcell_methylation/meffil_output/pcs.norm.beta.Robj")
norm.summary <- meffil.normalization.summary(norm.objects, pcs=pcs,parameters=norm.parameters)
meffil.normalization.report(norm.summary, output.file="/home/uy23281/data/bcell_methylation/meffil_output/normalization-report.html")



