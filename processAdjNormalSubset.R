# subset - only select those arrays and type from 'normal' tissue
# should be ~ 66
# file to create a TCGA expression set from excel file and demographic data.
# demographic data we are interested in call from PAM50
# note that this data has already been filtered and imputed
# will only need to filter by IQR during analysis:

#plan: 
# 1) Create expression set with demographic data
# 2) filter data based on median IQR
# 3) find overlap with signature(s)
# 4 - do Creighton method on all including Age signature.
#load libraries
library(Biobase)
library(gdata) #not all data is in excel format (which I hate)
setwd("/Users/mdarcy100/Desktop/MTroester/TCGA/")
#data located in AdjacentNormal

tcga.imputed.raw <- read.delim(file="AdjacentNormal/BRCA.exp.586.txt", header=TRUE, sep="\t",check.names=FALSE, stringsAsFactors=FALSE)
dim(tcga.imputed.raw)
#17815   589
tcga.imputed.exprs <- tcga.imputed.raw[2:nrow(tcga.imputed.raw),4:ncol(tcga.imputed.raw)]

probe.names <- (as.character(tcga.imputed.raw$CLID)[2:nrow(tcga.imputed.raw)])
rownames(tcga.imputed.exprs) <- probe.names

#tcga.imputed.exprs <- data.matrix(tcga.imputed.exprs)
dim(tcga.imputed.exprs) #586 in the exprs set
#17814   586




#read in paired data information
pairs.data <- read.xls("AdjacentNormal/TCGA.breast.upto103_pam50scores.for.Troester.xlsx",sheet="tumor.Normal.Pairs")
dim(pairs.data)
#[1] 124  20
pair.data.ind<-which(pairs.data$Type=='normal')
length(pair.data.ind)
#62
pairs.data<-pairs.data[pair.data.ind,]
pairs.data
colnames(pairs.data)
pairs.data<-pairs.data[,-(2:17)]
colnames(pairs.data)
#[1] "Sample"       "Type"         "Claudin.Call" "TumorCall"   

#subset the data with array name
all.patient.data <- read.xls("AdjacentNormal/TCGA.breast.upto103_pam50scores.for.Troester.xlsx",sheet="TCGA.breast.upto103_pam50scores")
sub.patient.data.ind<-which(all.patient.data[,1]%in% pairs.data$Sample)
length(sub.patient.data.ind)
#[1] 62
sub.patient.data<-all.patient.data[sub.patient.data.ind,]
dim(sub.patient.data)
# [1] 62 20

#get rid of all columns except for array name and sample
sub.patient.data<-sub.patient.data[,-2]
sub.patient.data<-sub.patient.data[,1:2]

#order both datasets by sample
merged.data<-merge(pairs.data,sub.patient.data,by.x="Sample",by.y="Sample")
rownames(merged.data) <-merged.data$ArraySample


####### FIRST EXTRACT OUT THOSE SAMPLES THAT ARE 'normal' to get identifer (Type = normal)
#normals_adj_id<-which(patient.data$Type=='normal')
###### get their name
#normals_sample_id<-patient.data[normals_adj_id,1]

length(intersect(rownames(merged.data),colnames(tcga.imputed.exprs)))
#[1] 60

patient.id.match <- rownames(merged.data)[which(rownames(merged.data) %in% colnames(tcga.imputed.exprs))]
patient.match <- merged.data[patient.id.match,]
array.match<-which(colnames(tcga.imputed.exprs)  %in%  rownames(merged.data))
tcga.imputed.exprs<-tcga.imputed.exprs[,array.match]
dim(tcga.imputed.exprs)
#[1] 17814    60
all(colnames(tcga.imputed.exprs) == rownames(patient.match))
# false - so need to sort
tcga.imputed.exprs<-tcga.imputed.exprs[,order(colnames(tcga.imputed.exprs))]
patient.match<-patient.match[order(rownames(patient.match)),]

all(colnames(tcga.imputed.exprs) == rownames(patient.match))
#returned true


length(unique(rownames(tcga.imputed.exprs)))==length(rownames(tcga.imputed.exprs))
#TRUE

################CREATE expression set
#there are no duplicates - i guess genes have already been collapsed
patient.match <- data.frame(patient.match)
tcga.imputed.exprs <- data.matrix(tcga.imputed.exprs)

tcga.sub.phenoData <- new("AnnotatedDataFrame", data = patient.match)

tcga.sub.eset <- new("ExpressionSet", exprs = data.matrix(tcga.imputed.exprs),
                   phenoData = tcga.sub.phenoData, annotation = "hgug4112a")
dim(tcga.sub.eset)
#Features  Samples 
#   17814      60 

save(tcga.sub.eset,file="AdjacentNormal/tcga_sub_eset.RData")
