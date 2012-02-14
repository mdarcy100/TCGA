library(Biobase)
library(hgug4112a.db)
library(gage)
library(gdata)
library(limma)
source("/Users/mdarcy100/Desktop/MTroester/AgePaper/AgeManuscript_062011/age_paper_scripts2.R")
setwd("/Users/mdarcy100/Desktop/MTroester/TCGA/")

load("AdjacentNormal/tcga_sub_eset.RData")
load("/Users/mdarcy100/Desktop/MTroester/AgePaper/AgeManuscript_062011/sig_fold_rm_eset.RData")
load("/Users/mdarcy100/Desktop/MTroester/AgePaper/AgeManuscript_062011/sig_full_rm_eset.RData") 
dim(tcga.sub.eset)


## Median center data
## first remove rows that have missing data

missing.tcga <- esApply(tcga.sub.eset,1,
                       function(x){sum(is.na(x))})
missing.inds.tcga <- which(missing.tcga > 0)
length(missing.inds.tcga)
#there were 168 rows with missing data (null).  strange because it looked like it was imputed
tcga.sub.eset<-tcga.sub.eset[-missing.inds.tcga,]

medians.tcga <- esApply(tcga.sub.eset,1,median)
tcga.exprs.scale <- t(scale(t(exprs(tcga.sub.eset)),center=medians.tcga,scale=FALSE))

#making sure that no missing data
which(is.na(tcga.exprs.scale[,1])==TRUE)
sum(is.na(exprs(tcga.sub.eset))[1,])
which(is.na(exprs(tcga.sub.eset)[1,])==TRUE)


## Create a new eset
tcga.scale.eset <- tcga.sub.eset
exprs(tcga.scale.eset) <- tcga.exprs.scale
dim(tcga.scale.eset)

## Remove low variability probes with by 50% IQR
iqr.val.tcga.scale <- esApply(tcga.sub.eset,1,function(x){IQR(x,na.rm = TRUE)})
high.iqr.inds.tcga <- which(iqr.val.tcga.scale > median(iqr.val.tcga.scale))
tcga.scale.eset <- tcga.scale.eset[high.iqr.inds.tcga,]
sum(is.na(exprs(tcga.scale.eset)[,1]))
colnames(exprs(tcga.scale.eset))

claudin_low <- read.xls("ErickSignatures/claudin low sig.xls",sheet="Sheet1") #Entrez Group2 Average Direction (807)
dtf <- read.xls("ErickSignatures/DTF sig.xls",sheet="Sheet1") #colnames = Direction, Symbol, Name, Entrez (787 rows)
E2<-read.xls("ErickSignatures/E2 response sig.xls",sheet="Sheet1") #Entrez, Group2, Average, Direction (754 rows)
twist<-read.xls("ErickSignatures/Twist sig.xls",sheet="Sheet1") #Column 1 = symbol, #Column 2 = description, Column 3, Entrez, Direction (838)
tgfb<-read.xls("ErickSignatures/TGFb sig.xls",sheet="Sheet1")#Name2 2 Column3 - description, Column4, Entrez, Direction
active<-read.xls("ErickSignatures/Active.xls",sheet="Sheet1") #NAME, Symbol, Centroid, Direction, Entrez (3518)
data(egSymb)

#convert to entrez id
tcga_entrez<-sym2eg(rownames(exprs(tcga.scale.eset)))
rownames(exprs(tcga.scale.eset))<-tcga_entrez #reset rowname to entrez to have comparable group
###################################################################################################################
#Claudin_Low
#first narrow down TCGA to those are in Claudin low signature 
###################################################################################################################
claudin_low$Entrez<- sapply(claudin_low$Entrez, as.character) # it was read in as numeric not string
claudinlow_OVERLAP.eset <- tcga.scale.eset[rownames(exprs(tcga.scale.eset))%in% claudin_low$Entrez,]
claudinlow_OVERLAP.eset$TumorCall<-trimWhiteSpace(claudinlow_OVERLAP.eset$TumorCall)
dim(claudinlow_OVERLAP.eset)
#Features  Samples 
#     492       60 

length(unique(rownames(exprs(claudinlow_OVERLAP.eset))))
#492

dim(claudin_low)
#[1] 807   4
claudin_low_sub<-claudin_low[claudin_low$Entrez %in% rownames(exprs(claudinlow_OVERLAP.eset)),]
dim(claudin_low_sub)
#[1] 492   4


length(intersect(claudin_low_sub$Entrez,rownames(exprs(claudinlow_OVERLAP.eset))))
length(unique(rownames(exprs(claudinlow_OVERLAP.eset))))
length(unique(claudin_low_sub$Entrez))
#all 492




#order both claudin_low and claudinlow_overlap.eset by entrez genes
#also order claudinlow_overlap.eset by $TumorCall

claudin_low_sub<-claudin_low_sub[order(claudin_low_sub$Entrez),]
claudinlow_OVERLAP.eset.order<-claudinlow_OVERLAP.eset[,order(claudinlow_OVERLAP.eset$TumorCall)]
#order by entrez id
exprs(claudinlow_OVERLAP.eset.order)<-exprs(claudinlow_OVERLAP.eset.order[order(rownames(exprs(claudinlow_OVERLAP.eset.order))),])
rownames(exprs(claudinlow_OVERLAP.eset.order))
length(intersect(claudin_low_sub$Entrez,rownames(exprs(claudinlow_OVERLAP.eset.order))))


#############Lum A
claudin_low_lumA <- claudinlow_OVERLAP.eset.order[,claudinlow_OVERLAP.eset.order$TumorCall %in% c("LumA", " LumA", "LumA ")]
dim(claudin_low_lumA) 
#Features  Samples 
#    492      25
claudinlow_cor_LumA<-cor(exprs(claudin_low_lumA),as.matrix(claudin_low_sub$Direction))
median(claudinlow_cor_LumA)
#-0.20
pdf('Images/Claudin_low/Adj_LumA_Correlations.pdf')
barplot(t(claudinlow_cor_LumA),main="Correlation, type=LumA, N=25, median = -0.20",names.arg=claudin_low_lumA$TumorCall)
dev.off()

 #############Lum B
claudin_low_lumB <- claudinlow_OVERLAP.eset.order[,claudinlow_OVERLAP.eset.order$TumorCall %in% c("LumB", " LumB", "LumB ")] 
dim(claudin_low_lumB) 
#Features  Samples 
#     492      12
claudinlow_cor_LumB<-cor(exprs(claudin_low_lumB),as.matrix(claudin_low_sub$Direction))
median(claudinlow_cor_LumB)
#0.17
pdf('Images/Claudin_low/Adj_LumB_Correlations.pdf')
barplot(t(claudinlow_cor_LumB),main="Correlation, type=LumB, N=12, median = 0.17",names.arg=claudin_low_lumB$TumorCall)
dev.off()

#############Basal
claudin_low_Basal <- claudinlow_OVERLAP.eset.order[,claudinlow_OVERLAP.eset.order$TumorCall %in% c("Basal", "Basal ", " Basal")] 
dim(claudin_low_Basal) 
#Features  Samples 
#  492       9
claudinlow_cor_Basal<-cor(exprs(claudin_low_Basal),as.matrix(claudin_low_sub$Direction))
median(claudinlow_cor_Basal)
# -0.18
pdf('Images/Claudin_low/Adj_Basal_Correlations.pdf')
barplot(t(claudinlow_cor_Basal),main="Correlation, type=Basal, N=9, median = -0.18",names.arg=claudin_low_Basal$TumorCall)
dev.off()

#############Her2
claudin_low_Her2 <- claudinlow_OVERLAP.eset.order[,claudinlow_OVERLAP.eset.order$TumorCall %in% c("Her2", " Her2", "Her2 ")] 
dim(claudin_low_Her2) 
#Features  Samples 
#     492       4
claudinlow_cor_Her2<-cor(exprs(claudin_low_Her2),as.matrix(claudin_low_sub$Direction))
median(claudinlow_cor_Her2)
# -0.15
pdf('Images/Claudin_low/Adj_Her2_Correlations.pdf')
barplot(t(claudinlow_cor_Her2),main="Correlation, type=Her2, N=4, median = -0.15",names.arg=claudin_low_Her2$TumorCall)
dev.off()


############Claudin low
claudin_low_claudin <- claudinlow_OVERLAP.eset.order[,claudinlow_OVERLAP.eset.order$TumorCall %in% c("Claudin")] 
dim(claudin_low_claudin) 
#Features  Samples 
#     492       10
claudinlow_cor_claudin<-cor(exprs(claudin_low_claudin),as.matrix(claudin_low_sub$Direction))
median(claudinlow_cor_claudin)
# 0.56
pdf('Images/Claudin_low/Adj_Claudin_Correlations.pdf')
barplot(t(claudinlow_cor_claudin),main="Correlation, type=Claudin, N=10, median = 0.56",names.arg=claudin_low_claudin$TumorCall)
dev.off()

########### FUll correlations
claudinlow_cor<-cor(exprs(claudinlow_OVERLAP.eset.order),as.matrix(claudin_low_sub$Direction))
pdf('Images/Claudin_low/Full_Correlations.pdf')
barplot(t(claudinlow_cor), main ="Correlation with Claudin low  - sorted by subtype", names.arg=claudinlow_OVERLAP.eset.order$TumorCall)  
dev.off()


length(which(claudinlow_cor_Her2 > 0))
length(which(claudinlow_cor_claudin> 0))
length(which(claudinlow_cor_Basal > 0))
length(which(claudinlow_cor_LumB > 0))
length(which(claudinlow_cor_LumA > 0))


###################################################################################################
#ACTIVE
#first narrow down TCGA to those are in Active signature 
#########################################################################################
# entrez id was read in as numeric type - 
#need to convert to character for meaningful comparison
#there was trailing white space - need to remove

active$Entrez<- sapply(active$Entrez, as.character)
active$Entrez<-trimWhiteSpace(active$Entrez)
intersect(active$Entrez,rownames(exprs(tcga.scale.eset)))
head(active$Entrez)
active<-active[order(active$Entrez),]

active_OVERLAP.eset <- tcga.scale.eset[rownames(exprs(tcga.scale.eset))%in% active$Entrez,]
active_OVERLAP.eset$TumorCall<-trimWhiteSpace(active_OVERLAP.eset$TumorCall)
dim(active_OVERLAP.eset)
#Features  Samples 
#     2619       60 

length(unique(rownames(exprs(active_OVERLAP.eset))))
#2619

dim(active)
#[1] 3518   4
active_sub<-active[active$Entrez %in% rownames(exprs(active_OVERLAP.eset)),]
dim(active_sub)
#[1] 2621   4 #need to collapse a bit


active_col<-tapply(as.vector(active_sub$Direction), factor(active_sub$Entrez), mean)
length(active_col)
#2619
active_col<-as.matrix(active_col)


active_sub<-as.matrix(active_col[rownames(active_col) %in% rownames(exprs(active_OVERLAP.eset)),]) #is not necessary
which(is.na(exprs(active_OVERLAP.eset)[,1])==TRUE)

#order both by entrez genes

active_OVERLAP.eset.order<-active_OVERLAP.eset[,order(active_OVERLAP.eset$TumorCall)]
head(active_OVERLAP.eset.order$TumorCall)
#order by entrez id
exprs(active_OVERLAP.eset.order)<-exprs(active_OVERLAP.eset.order[order(rownames(exprs(active_OVERLAP.eset.order))),])

#make sure sorted
head(rownames(exprs(active_OVERLAP.eset.order)))
length(intersect(rownames(active_sub),rownames(exprs(active_OVERLAP.eset.order))))


#############Lum A
active_lumA <- active_OVERLAP.eset.order[,active_OVERLAP.eset.order$TumorCall %in% c("LumA")] 
dim(active_lumA)
#Features  Samples 
#     2619      25
active_cor_LumA<-cor(exprs(active_lumA),active_sub[,1])
median(active_cor_LumA)
#-0.33
pdf('Images/Active/Adj_LumA_Correlations.pdf')
barplot(t(active_cor_LumA),main="Correlation, type=LumA, N=25, median = -0.33",names.arg=active_lumA$TumorCall)
dev.off()

#############Lum B
active_lumB <- active_OVERLAP.eset.order[,active_OVERLAP.eset.order$TumorCall %in% c("LumB")] 
dim(active_lumB)
#Features  Samples 
#     2619      12
active_cor_LumB<-cor(exprs(active_lumB),active_sub[,1])
median(active_cor_LumB)

#0.08
pdf('Images/Active/Adj_LumB_Correlations.pdf')
barplot(t(active_cor_LumB),main="Correlation, type=LumB, N=12, median = 0.08",names.arg=active_lumB$TumorCall)
dev.off()

#############Basal
active_basal <- active_OVERLAP.eset.order[,active_OVERLAP.eset.order$TumorCall %in% c("Basal")] 
dim(active_basal)
#Features  Samples 
#     2619      9
active_cor_basal<-cor(exprs(active_basal),active_sub[,1])
median(active_cor_basal)

#-0.19
pdf('Images/Active/Adj_Basal_Correlations.pdf')
barplot(t(active_cor_basal),main="Correlation, type=Basal, N=9, median = -0.19",names.arg=active_basal$TumorCall)
dev.off()

#############Her2
active_her2 <- active_OVERLAP.eset.order[,active_OVERLAP.eset.order$TumorCall %in% c("Her2")] 
dim(active_her2)
#Features  Samples 
#     2619      4
active_cor_her2<-cor(exprs(active_her2),active_sub[,1])
median(active_cor_her2)

#-0.23
pdf('Images/Active/Adj_Her2_Correlations.pdf')
barplot(t(active_cor_her2),main="Correlation, type=Her2, N=4, median = -0.23",names.arg=active_her2$TumorCall)
dev.off()


############Claudin low
active_claudin <- active_OVERLAP.eset.order[,active_OVERLAP.eset.order$TumorCall %in% c("Claudin")] 
dim(active_claudin) 
#Features  Samples 
#     2619       10
active_cor_claudin<-cor(exprs(active_claudin),active_sub[,1])
median(active_cor_claudin)
# 0.63
pdf('Images/Active/Adj_Claudin_Correlations.pdf')
barplot(t(active_cor_claudin),main="Correlation, type=Claudin, N=10, median = 0.63",names.arg=active_claudin$TumorCall)
dev.off()

########### FUll correlations
active_cor<-cor(exprs(active_OVERLAP.eset.order),active_sub[,1])
pdf('Images/Active/Full_Correlations.pdf')
barplot(t(active_cor), main ="Correlation with Active  - sorted by subtype", names.arg=active_OVERLAP.eset.order$TumorCall)  
dev.off()


# proportions of +/-
length(which(active_cor_her2 > 0))
length(which(active_cor_her2 < 0))
length(which(active_cor_claudin> 0))
length(which(active_cor_claudin < 0))
length(which(active_cor_basal < 0))
length(which(active_cor_basal > 0))
length(which(active_cor_LumB < 0))
length(which(active_cor_LumB > 0))
length(which(active_cor_LumA < 0))
length(which(active_cor_LumA > 0))

###################################################################################################################
#E2 signature
#first narrow down TCGA to those are in E2 signature 
###################################################################################################################
E2_OVERLAP.eset <- tcga.scale.eset[rownames(exprs(tcga.scale.eset))%in% E2$Entrez,]
E2_OVERLAP.eset$TumorCall<-trimWhiteSpace(E2_OVERLAP.eset$TumorCall)
dim(E2_OVERLAP.eset)
#Features  Samples 
#     366      60
# check to see if there are any repeats
length(unique(rownames(exprs(E2_OVERLAP.eset))))
#366 in this list; don't need to collapse genes

dim(E2)
#[1] 754   4
E2_sub<-E2[E2$Entrez %in% rownames(exprs(E2_OVERLAP.eset)),]
dim(E2_sub)
#[1] 366   4

#make sure the length is the same/complete
length(intersect(E2_sub$Entrez,rownames(exprs(E2_OVERLAP.eset))))
#[1] 366
length(unique(rownames(exprs(E2_OVERLAP.eset))))
length(unique(E2_sub$Entrez))

# entrez id was read in as numeric type (Entrez in first column)
E2_sub$Entrez<- sapply(E2_sub[,1], as.character)
E2_sub<-E2_sub[order(E2_sub$Entrez),]

#order both by entrez genes
#also order overlap.eset by $TumorCall
E2_OVERLAP.eset.order<-E2_OVERLAP.eset[,order(E2_OVERLAP.eset$TumorCall)]

#order by entrez id
exprs(E2_OVERLAP.eset.order)<-exprs(E2_OVERLAP.eset.order[order(rownames(exprs(E2_OVERLAP.eset.order))),])
head(rownames(exprs(E2_OVERLAP.eset.order)))
length(intersect(E2_sub$Entrez,rownames(exprs(E2_OVERLAP.eset.order))))
#366

# do correlation - first extract out the subtypes

#############Lum A
E2_lumA <- E2_OVERLAP.eset.order[,E2_OVERLAP.eset.order$TumorCall %in% c("LumA")] 
dim(E2_lumA)
#Features  Samples 
#     366      32
E2_cor_LumA<-cor(exprs(E2_lumA),as.matrix(E2_sub$Direction))
median(E2_cor_LumA)
#0.03
pdf('Images/E2/LumA_Correlations.pdf')
barplot(t(E2_cor_LumA),main="Correlation, type=LumA, N=25, median = 0.03",names.arg=E2_lumA$TumorCall)
dev.off()

#############Lum B
E2_lumB <- E2_OVERLAP.eset.order[,E2_OVERLAP.eset.order$TumorCall %in% c("LumB")] 
dim(E2_lumB) 
#Features  Samples 
#     482      12
E2_cor_LumB<-cor(exprs(E2_lumB),as.matrix(E2_sub$Direction))
median(E2_cor_LumB)
#0.05
pdf('Images/E2/LumB_Correlations.pdf')
barplot(t(E2_cor_LumB),main="Correlation, type=LumB, N=12, median = 0.05",names.arg=E2_lumB$TumorCall)
dev.off()

#############Basal
E2_Basal <- E2_OVERLAP.eset.order[,E2_OVERLAP.eset.order$TumorCall %in% c("Basal")] 
dim(E2_Basal) 
#Features  Samples 
#  366       9 
E2_cor_Basal<-cor(exprs(E2_Basal),as.matrix(E2_sub$Direction))
median(E2_cor_Basal)
# 0.08
pdf('Images/E2/Basal_Correlations.pdf')
barplot(t(E2_cor_Basal),main="Correlation, type=Basal, N=9, median = 0.08",names.arg=E2_Basal$TumorCall)
dev.off()

#############Her2
E2_Her2 <- E2_OVERLAP.eset.order[,E2_OVERLAP.eset.order$TumorCall %in% c("Her2")] 
dim(E2_Her2) 
#Features  Samples 
#     366       4
E2_cor_Her2<-cor(exprs(E2_Her2),as.matrix(E2_sub$Direction))
median(E2_cor_Her2)
# 0.03
pdf('Images/E2/Her2_Correlations.pdf')
barplot(t(E2_cor_Her2),main="Correlation, type=Her2, N=4, median = 0.03",names.arg=E2_Her2$TumorCall)
dev.off()


############Claudin low
E2_claudin <- E2_OVERLAP.eset.order[,E2_OVERLAP.eset.order$TumorCall %in% c("Claudin")] 
dim(E2_claudin) 
#Features  Samples 
#     2619       10
E2_cor_claudin<-cor(exprs(E2_claudin),as.matrix(E2_sub$Direction))
median(E2_cor_claudin)
# -0.15
pdf('Images/E2/Adj_Claudin_Correlations.pdf')
barplot(t(E2_cor_claudin),main="Correlation, type=Claudin, N=10, median = -0.15",names.arg=E2_claudin$TumorCall)
dev.off()


###############Full correlations
E2_cor<-cor(exprs(E2_OVERLAP.eset.order),as.matrix(E2_sub$Direction))
pdf('Images/E2/Full_Correlations.pdf')
barplot(t(E2_cor), main ="Correlation  - sorted by subtype",names.arg=E2_OVERLAP.eset.order$TumorCall)  
dev.off()

#what proportions are +/- associated

length(which(E2_cor_LumA>0))
length(which(E2_cor_LumB>0))
length(which(E2_cor_claudin>0))
length(which(E2_cor_Basal>0))
length(which(E2_cor_Her2>0))

###################################################################################################################
#  Twist signature
###################################################################################################################
twist<-read.xls("ErickSignatures/Twist sig.xls",sheet="Sheet1") #Column 1 = symbol, #Column 2 = description, Column 3, Entrez, Direction (838)
twist<-twist[,-(1:3)]
twist$Entrez<- sapply(twist$Entrez, as.character)
twist$Entrez<-trimWhiteSpace(twist$Entrez)

twist_OVERLAP.eset <- tcga.scale.eset[rownames(exprs(tcga.scale.eset))%in% twist$Entrez,]
twist_OVERLAP.eset$TumorCall<-trimWhiteSpace(twist_OVERLAP.eset$TumorCall)
dim(twist_OVERLAP.eset)
#Features  Samples 
#     633      586 
# check to see if there are any repeats
length(unique(rownames(exprs(twist_OVERLAP.eset))))
#633 in this list; 

#subset the original dtf
dim(twist)
#[1] 838   2
length(unique(twist$Entrez))
# 834
twist_col<-tapply(as.vector(twist[,2]), factor(twist[,1]), mean)
twist_col<-as.matrix(twist_col)

twist_sub<-as.matrix(twist_col[rownames(twist_col) %in% rownames(exprs(twist_OVERLAP.eset)),])

head(colnames(exprs(twist_OVERLAP.eset)))
head(rownames(exprs(twist_OVERLAP.eset)))

length(twist_sub)
# 633
length(unique(rownames(twist_sub)))
#[1] 633
#there must be a lot of duplicates
#make sure the length is the same/complete
length(intersect(rownames(twist_sub),rownames(exprs(twist_OVERLAP.eset))))
#[1] 633

twist_sub<-twist_sub[order(rownames(twist_sub)),]
twist_sub<-as.matrix(twist_sub)

#order by entrez and TumorCall
twist_OVERLAP.eset.order<-twist_OVERLAP.eset[,order(twist_OVERLAP.eset$TumorCall)]
#order by entrez id
exprs(twist_OVERLAP.eset.order)<-exprs(twist_OVERLAP.eset.order[order(rownames(exprs(twist_OVERLAP.eset.order))),])
head(rownames(exprs(twist_OVERLAP.eset.order)))
length(intersect(rownames(twist_sub),rownames(exprs(twist_OVERLAP.eset.order))))
#633



#############Lum A
twist_lumA <- twist_OVERLAP.eset.order[,twist_OVERLAP.eset.order$TumorCall %in% c("LumA")] 
dim(twist_lumA)
#Features  Samples 
#     633      25
twist_cor_LumA<-cor(exprs(twist_lumA),twist_sub[,1])
median(twist_cor_LumA)
#-0.19
pdf('Images/Twist/LumA_Correlations.pdf')
barplot(t(twist_cor_LumA),main="Correlation, type=LumA, N=25, median = -0.19",names.arg=twist_lumA$TumorCall)
dev.off()

#############Lum B
twist_lumB <- twist_OVERLAP.eset.order[,twist_OVERLAP.eset.order$TumorCall %in% c("LumB")] 
dim(twist_lumB) 
#Features  Samples 
#     633      12
twist_cor_LumB<-cor(exprs(twist_lumB),twist_sub[,1])
median(twist_cor_LumB)
#0.05
pdf('Images/Twist/LumB_Correlations.pdf')
barplot(t(twist_cor_LumB),main="Correlation, type=LumB, N=12, median = 0.05",names.arg=twist_lumB$TumorCall)
dev.off()

#############Basal
twist_Basal <- twist_OVERLAP.eset.order[,twist_OVERLAP.eset.order$TumorCall %in% c("Basal")] 
dim(twist_Basal) 
#Features  Samples 
#  633      9
twist_cor_Basal<-cor(exprs(twist_Basal),twist_sub[,1])
median(twist_cor_Basal)
# -0.14
pdf('Images/Twist/Basal_Correlations.pdf')
barplot(t(twist_cor_Basal),main="Correlation, type=Basal, N=9, median = -0.14",names.arg=twist_Basal$TumorCall)
dev.off()


#############Her2
twist_Her2 <- twist_OVERLAP.eset.order[,twist_OVERLAP.eset.order$TumorCall %in% c("Her2")] 
dim(twist_Her2) 
#Features  Samples 
#     633       4
twist_cor_Her2<-cor(exprs(twist_Her2),twist_sub[,1])
median(twist_cor_Her2)
# -0.08
pdf('Images/Twist/Her2_Correlations.pdf')
barplot(t(twist_cor_Her2),main="Correlation, type=Her2, N=4, median = -0.08",names.arg=twist_Her2$TumorCall)
dev.off()


############Claudin low
twist_claudin <- twist_OVERLAP.eset.order[,twist_OVERLAP.eset.order$TumorCall %in% c("Claudin")] 
dim(twist_claudin) 
#Features  Samples 
#     633       10
twist_cor_claudin<-cor(exprs(twist_claudin),twist_sub[,1])
median(twist_cor_claudin)
# 0.40
pdf('Images/Twist/Adj_Claudin_Correlations.pdf')
barplot(t(active_cor_claudin),main="Correlation, type=Claudin, N=10, median = 0.40",names.arg=twist_claudin$TumorCall)
dev.off()


###############Full correlations
twist_cor<-cor(exprs(twist_OVERLAP.eset.order),twist_sub[,1])
pdf('Images/Twist/Full_Correlations.pdf')
barplot(t(twist_cor), main ="Correlation, Twist - sorted by subtype",names.arg=twist_OVERLAP.eset.order$TumorCall)  
dev.off()


# look at proportions that are concordant
length(which(twist_cor_LumA>0))
length(which(twist_cor_LumB>0))
length(which(twist_cor_Basal>0))
length(which(twist_cor_Her2>0))
length(which(twist_cor_claudin>0))
###################################################################################################################
#  DTF signature
# extract out only the direction and entrez id
# will need to convert entrez id to character
###################################################################################################################
dtf<-dtf[,-(2:3)]
dtf$Entrez<- sapply(dtf$Entrez, as.character)
dtf$Entrez<-trimWhiteSpace(dtf$Entrez)
dtf_OVERLAP.eset <- tcga.scale.eset[rownames(exprs(tcga.scale.eset))%in% dtf$Entrez,]
dtf_OVERLAP.eset$TumorCall<-trimWhiteSpace(dtf_OVERLAP.eset$TumorCall)
dim(dtf_OVERLAP.eset)
#Features  Samples 
#     362      60 
# check to see if there are any repeats
length(unique(rownames(exprs(dtf_OVERLAP.eset))))
#362 

#subset the original dtf
dim(dtf)
#[1] 786   2
length(unique(dtf$Entrez))
# 512

#Direction, Entrez is order of columns, collpsing
dtf_col<-tapply(as.vector(dtf[,1]), factor(dtf[,2]), mean)
dtf_col<-as.matrix(dtf_col)
dtf_sub<-as.matrix(dtf_col[rownames(dtf_col) %in% rownames(exprs(dtf_OVERLAP.eset)),])
length(dtf_sub)
# 362
length(unique(rownames(dtf_sub)))
#[1] 362

#order by entrez id
dtf_sub<-dtf_sub[order(rownames(dtf_sub)),]
dtf_sub<-as.matrix(dtf_sub)

#order both by entrez genes
#also order overlap.eset by $TumorCall
dtf_OVERLAP.eset.order<-dtf_OVERLAP.eset[,order(dtf_OVERLAP.eset$TumorCall)]

#order by entrez id
exprs(dtf_OVERLAP.eset.order)<-exprs(dtf_OVERLAP.eset.order[order(rownames(exprs(dtf_OVERLAP.eset.order))),])
length(intersect(rownames(dtf_sub),rownames(exprs(dtf_OVERLAP.eset.order))))
#362

#############Lum A
dtf_lumA <- dtf_OVERLAP.eset.order[,dtf_OVERLAP.eset.order$TumorCall %in% c("LumA")] 
dim(dtf_lumA)
#Features  Samples 
#     362      25
dtf_cor_LumA<-cor(exprs(dtf_lumA),dtf_sub[,1])
median(dtf_cor_LumA)
#-0.01
pdf('Images/DTF/LumA_Correlations.pdf')
barplot(t(dtf_cor_LumA),main="Correlation, type=LumA, N=25, median = -0.01",names.arg=dtf_lumA$TumorCall)
dev.off()

#############Lum B
dtf_lumB <- dtf_OVERLAP.eset.order[,dtf_OVERLAP.eset.order$TumorCall %in% c("LumB")] 
dim(dtf_lumB)
#Features  Samples 
#     362      12
dtf_cor_LumB<-cor(exprs(dtf_lumB),dtf_sub[,1])
median(dtf_cor_LumB)
#-0.002
pdf('Images/DTF/LumB_Correlations.pdf')
barplot(t(dtf_cor_LumB),main="Correlation, type=LumB, N=12, median = -0.002",names.arg=dtf_lumB$TumorCall)
dev.off()

#############Basal
dtf_Basal <- dtf_OVERLAP.eset.order[,dtf_OVERLAP.eset.order$TumorCall %in% c("Basal")] 
dim(dtf_Basal)
#Features  Samples 
#     362      9
dtf_cor_Basal<-cor(exprs(dtf_Basal),dtf_sub[,1])
median(dtf_cor_Basal)
#0.02
pdf('Images/DTF/Basal_Correlations.pdf')
barplot(t(dtf_cor_Basal),main="Correlation, type=Basal, N=9, median = 0.02",names.arg=dtf_Basal$TumorCall)
dev.off()

#############Her2
dtf_Her2 <- dtf_OVERLAP.eset.order[,dtf_OVERLAP.eset.order$TumorCall %in% c("Her2")] 
dim(dtf_Her2)
#Features  Samples 
#     362      4
dtf_cor_Her2<-cor(exprs(dtf_Her2),dtf_sub[,1])
median(dtf_cor_Her2)
#-0.05
pdf('Images/DTF/Her2_Correlations.pdf')
barplot(t(dtf_cor_Her2),main="Correlation, type=Her2, N=4, median = -0.04")
dev.off()

####################claudin
dtf_claudin <- dtf_OVERLAP.eset.order[,dtf_OVERLAP.eset.order$TumorCall %in% c("Claudin")] 
dim(dtf_claudin) 
#Features  Samples 
#     362       10
dtf_cor_claudin<-cor(exprs(dtf_claudin),dtf_sub[,1])
median(dtf_cor_claudin)
# 0.06
pdf('Images/DTF/Adj_Claudin_Correlations.pdf')
barplot(t(active_cor_claudin),main="Correlation, type=Claudin, N=10, median = 0.06",names.arg=twist_claudin$TumorCall)
dev.off()


###############Full correlations
dtf_cor<-cor(exprs(dtf_OVERLAP.eset.order),dtf_sub[,1])
pdf('Images/DTF/Full_Correlations.pdf')
barplot(t(dtf_cor), main ="Correlation DTF - sorted by subtype",names.arg=dtf_OVERLAP.eset.order$TumorCall)  
dev.off()


length(which(dtf_cor_LumA>0))
length(which(dtf_cor_LumB>0))
length(which(dtf_cor_Basal>0))
length(which(dtf_cor_Her2>0))
length(which(dtf_cor_claudin>0))

###################################################################################################################
#TGFB stuff - will also need to check to see if this has been collapsed and trim whitespace
###################################################################################################################

tgfb<-tgfb[,-(1:3)] #only need entrez and direction

missing.tgfb <- apply(tgfb,1,function(x){sum(is.na(x))})
missing.inds.tgfb <- which(missing.tgfb > 0)
#there was one row with missing data
tgfb<-tgfb[-missing.inds.tgfb,]


tgfb$Entrez<- sapply(tgfb$Entrez, as.character)
tgfb$Entrez<-trimWhiteSpace(tgfb$Entrez)
tgfb_OVERLAP.eset <- tcga.scale.eset[rownames(exprs(tcga.scale.eset))%in% tgfb$Entrez,]
tgfb_OVERLAP.eset$TumorCall<-trimWhiteSpace(tgfb_OVERLAP.eset$TumorCall)

dim(tgfb_OVERLAP.eset)
#Features  Samples 
#     167      60
# check to see if there are any repeats
length(unique(rownames(exprs(tgfb_OVERLAP.eset))))
#167 

dim(tgfb)
#[1] 234   4
length(unique(tgfb$Entrez))
# 217
#Entrez, Direction is order of columns
tgfb_col<-tapply(as.vector(tgfb[,2]), factor(tgfb[,1]), mean)
tgfb_col<-as.matrix(tgfb_col)
tgfb_sub<-as.matrix(tgfb_col[rownames(tgfb_col) %in% rownames(exprs(tgfb_OVERLAP.eset)),])
length(tgfb_sub)
# 167
length(unique(rownames(tgfb_sub)))
#[1] 167

#order by entrez id
tgfb_sub<-tgfb_sub[order(rownames(tgfb_sub)),]
tgfb_sub<-as.matrix(tgfb_sub)


#order both by entrez genes
#also order overlap.eset by $Call
tgfb_OVERLAP.eset.order<-tgfb_OVERLAP.eset[,order(tgfb_OVERLAP.eset$TumorCall)]
#order by entrez id
exprs(tgfb_OVERLAP.eset.order)<-exprs(tgfb_OVERLAP.eset.order[order(rownames(exprs(tgfb_OVERLAP.eset.order))),])
length(intersect(rownames(tgfb_sub),rownames(exprs(tgfb_OVERLAP.eset.order))))



#############Lum A
tgfb_lumA <- tgfb_OVERLAP.eset.order[,tgfb_OVERLAP.eset.order$TumorCall %in% c("LumA")] 
dim(tgfb_lumA)
#Features  Samples 
#     167     25
tgfb_cor_LumA<-cor(exprs(tgfb_lumA),tgfb_sub[,1])
median(tgfb_cor_LumA)
#-0.09
pdf('Images/TGFB/LumA_Correlations.pdf')
barplot(t(tgfb_cor_LumA),main="Correlation, type=LumA, N=25, median = -0.09",names.arg=tgfb_lumA$TumorCall)
dev.off()

#############Lum B
tgfb_lumB <- tgfb_OVERLAP.eset.order[,tgfb_OVERLAP.eset.order$TumorCall %in% c("LumB")] 
dim(tgfb_lumB)
#Features  Samples 
#     167      12
tgfb_cor_LumB<-cor(exprs(tgfb_lumB),tgfb_sub[,1])
median(tgfb_cor_LumB)
#0.05
pdf('Images/TGFB/LumB_Correlations.pdf')
barplot(t(tgfb_cor_LumB),main="Correlation, type=LumB, N=12, median = 0.05",names.arg=tgfb_lumB$TumorCall)
dev.off()

#############Basal
tgfb_Basal <- tgfb_OVERLAP.eset.order[,tgfb_OVERLAP.eset.order$TumorCall %in% c("Basal")] 
dim(tgfb_Basal)
#Features  Samples 
#     167     9
tgfb_cor_Basal<-cor(exprs(tgfb_Basal),tgfb_sub[,1])
median(tgfb_cor_Basal)
#0.004
pdf('Images/TGFB/Basal_Correlations.pdf')
barplot(t(tgfb_cor_Basal),main="Correlation, type=Basal, N=11, median = 0.004",names.arg=tgfb_Basal$TumorCall)
dev.off()

#############Her2
tgfb_Her2 <- tgfb_OVERLAP.eset.order[,tgfb_OVERLAP.eset.order$TumorCall %in% c("Her2")] 
dim(tgfb_Her2)
#Features  Samples 
#    167      4
tgfb_cor_Her2<-cor(exprs(tgfb_Her2),tgfb_sub[,1])
median(tgfb_cor_Her2)
#-0.06
pdf('Images/TGFB/Her2_Correlations.pdf')
barplot(t(tgfb_cor_Her2),main="Correlation, type=Her2, N=4, median = -.06",names.arg=tgfb_Her2$TumorCall)
dev.off()

####################claudin
tgfb_claudin <- tgfb_OVERLAP.eset.order[,tgfb_OVERLAP.eset.order$TumorCall %in% c("Claudin")] 
dim(tgfb_claudin) 
#Features  Samples 
#     167       10
tgfb_cor_claudin<-cor(exprs(tgfb_claudin),tgfb_sub[,1])
median(tgfb_cor_claudin)
# 0.12
pdf('Images/TGFB/Adj_Claudin_Correlations.pdf')
barplot(t(active_cor_claudin),main="Correlation: type=Claudin, N=10, median = 0.12",names.arg=tgfb_claudin$TumorCall)
dev.off()



###############Full correlations
tgfb_cor<-cor(exprs(tgfb_OVERLAP.eset.order),tgfb_sub[,1])
pdf('Images/TGFB/Full_Correlations.pdf')
barplot(t(tgfb_cor), main ="Correlation: TGFB - sorted by subtype", names.arg=tgfb_OVERLAP.eset.order$TumorCall)  
dev.off()


length(which(tgfb_cor_LumA>0))
length(which(tgfb_cor_LumB>0))
length(which(tgfb_cor_Basal>0))
length(which(tgfb_cor_Her2>0))
length(which(tgfb_cor_claudin>0))



###################################################################################################################
############### AGE SIGNATURE
###################################################################################################################

age <- read.delim("/Users/mdarcy100/Desktop/MTroester/AgePaper/AgeManuscript_062011/Signature/AgeFullUpDown_2011-12-13.txt",sep="\t") #Entrez Group2 Average Direction (807)
age<-age[,(2:3)]
age$entrezid<- sapply(age$entrezid, as.character)
age$entrezid<-trimWhiteSpace(age$entrezid)
age_OVERLAP.eset <- tcga.scale.eset[rownames(exprs(tcga.scale.eset))%in% age$entrezid,]
age_OVERLAP.eset$TumorCall<-trimWhiteSpace(age_OVERLAP.eset$TumorCall)
dim(age_OVERLAP.eset)
#445       60 

length(unique(rownames(exprs(age_OVERLAP.eset))))
#445 

dim(age)
#[1] 802   4
length(unique(age$entrezid))
# 719
#Entrez, Direction is order of columns
age_col<-tapply(as.vector(age$direction), factor(age$entrezid), mean)
age_col<-as.matrix(age_col)
age_sub<-as.matrix(age_col[rownames(age_col) %in% rownames(exprs(age_OVERLAP.eset)),])
length(age_sub)
# 445
length(unique(rownames(age_sub)))
#[1] 445

#order by entrez id
age_sub<-age_sub[order(rownames(age_sub)),]
age_sub<-as.matrix(age_sub)


#order both by entrez genes
#also order overlap.eset by $TumorCall
age_OVERLAP.eset.order<-age_OVERLAP.eset[,order(age_OVERLAP.eset$TumorCall)]
#order by entrez id
exprs(age_OVERLAP.eset.order)<-exprs(age_OVERLAP.eset.order[order(rownames(exprs(age_OVERLAP.eset.order))),])
length(intersect(rownames(age_sub),rownames(exprs(age_OVERLAP.eset.order))))
#445

#############Lum A
age_lumA <- age_OVERLAP.eset.order[,age_OVERLAP.eset.order$TumorCall %in% c("LumA")] 
dim(age_lumA)
#Features  Samples 
#     445      25
age_cor_LumA<-cor(exprs(age_lumA),age_sub[,1])
median(age_cor_LumA)
#0.39
pdf('Images/Age/LumA_Correlations.pdf')
barplot(t(age_cor_LumA),main="Correlation, type=LumA, N=32, median = 0.39",names.arg=age_lumA$TumorCall)
dev.off()

#############Lum B
age_lumB <- age_OVERLAP.eset.order[,age_OVERLAP.eset.order$TumorCall %in% c("LumB")] 
dim(age_lumB)
#Features  Samples 
#     445      12
age_cor_LumB<-cor(exprs(age_lumB),age_sub[,1])
median(age_cor_LumB)
#-0.37
pdf('Images/Age/LumB_Correlations.pdf')
barplot(t(age_cor_LumB),main="Correlation, type=LumB, N=12, median = -0.37",names.arg=age_lumB$TumorCall)
dev.off()

#############Basal
age_Basal <- age_OVERLAP.eset.order[,age_OVERLAP.eset.order$TumorCall %in% c("Basal")] 
dim(age_Basal)
#Features  Samples 
#     445      9
age_cor_Basal<-cor(exprs(age_Basal),age_sub[,1])
median(age_cor_Basal)
#0.24
pdf('Images/Age/Basal_Correlations.pdf')
barplot(t(age_cor_Basal),main="Correlation, type=Basal, N=9, median = 0.24",names.arg=age_Basal$TumorCall)
dev.off()

#############Her2
age_Her2 <- age_OVERLAP.eset.order[,age_OVERLAP.eset.order$TumorCall %in% c("Her2")] 
dim(age_Her2)
#Features  Samples 
#    445      4
age_cor_Her2<-cor(exprs(age_Her2),age_sub[,1])
median(age_cor_Her2)
#0.20
pdf('Images/Age/Her2_Correlations.pdf')
barplot(t(age_cor_Her2),main="Correlation, type=Her2, N=4, median = 0.20",names.arg=age_Her2$TumorCall)
dev.off()

####################claudin
age_claudin <- age_OVERLAP.eset.order[,age_OVERLAP.eset.order$TumorCall %in% c("Claudin")] 
dim(age_claudin) 
#Features  Samples 
#     445       10
age_cor_claudin<-cor(exprs(age_claudin),age_sub[,1])
median(age_cor_claudin)
# -0.63
pdf('Images/AGE/Adj_Claudin_Correlations.pdf')
barplot(t(age_cor_claudin),main="Correlation: type=Claudin, N=10, median = -0.63",names.arg=age_claudin$TumorCall)
dev.off()

length(which(age_cor_LumA>0))
length(which(age_cor_LumB>0))
length(which(age_cor_Basal>0))
length(which(age_cor_Her2>0))
length(which(age_cor_claudin>0))

tmp_age<-exprs(age_OVERLAP.eset.order)
colnames(tmp_age)<-age_OVERLAP.eset.order$TumorCall
write.table(tmp_age,file=paste("TCGA_AgeSig_", Sys.Date(),".txt",sep=""),sep="\t",col.names=NA)


###############Full correlations
age_cor<-cor(exprs(age_OVERLAP.eset.order),age_sub[,1])
pdf('Images/Age/Full_Correlations.pdf')
barplot(t(age_cor), main ="Correlation: Age - sorted by subtype", names.arg=age_OVERLAP.eset.order$TumorCall)  
dev.off()

###################################################################################################################
# intersection
###################################################################################################################

#transform to all 1 / -1
age_cor_cast<-as.matrix(apply(age_cor,1,castDirection))
table(age_cor_cast)

tgfb_cor_cast<-as.matrix(apply(tgfb_cor,1,castDirection))
table(tgfb_cor_cast)

dtf_cor_cast<-as.matrix(apply(dtf_cor,1,castDirection))
table(tgfb_cor_cast)

twist_cor_cast<-as.matrix(apply(twist_cor,1,castDirection))
table(tgfb_cor_cast)

E2_cor_cast<-as.matrix(apply(E2_cor,1,castDirection))
table(E2_cor_cast)
dim(E2_sub)

active_cor_cast<-as.matrix(apply(active_cor,1,castDirection))
table(active_cor_cast)

claudinlow_cor_cast<-as.matrix(apply(claudinlow_cor,1,castDirection))
table(claudinlow_cor_cast)

#make sure everything is still sorted properly before doing chi-square test
all(rownames(claudinlow_cor_cast) == rownames(active_cor_cast))
all(rownames(E2_cor_cast) == rownames(active_cor_cast))
all(rownames(twist_cor_cast) == rownames(active_cor_cast))
all(rownames(dtf_cor_cast) == rownames(active_cor_cast))
all(rownames(tgfb_cor_cast) == rownames(active_cor_cast))
all(rownames(age_cor_cast) == rownames(active_cor_cast))


table(claudinlow_cor_cast,active_cor_cast)
chisq.test(claudinlow_cor_cast,active_cor_cast)
#X-squared = 34.3546, df = 1, p-value = 4.593e-09

table(age_cor_cast,active_cor_cast)
chisq.test(age_cor_cast,active_cor_cast)
#X-squared = 36.8428, df = 1, p-value = 1.28e-09

table(tgfb_cor_cast,active_cor_cast)
chisq.test(tgfb_cor_cast,active_cor_cast)
#X-squared = 20.2301, df = 1, p-value = 6.866e-06

table(dtf_cor_cast,active_cor_cast)
chisq.test(dtf_cor_cast,active_cor_cast)
#X-squared = 2.2117, df = 1, p-value = 0.137

table(twist_cor_cast,active_cor_cast)
chisq.test(twist_cor_cast,active_cor_cast)
#X-squared = 38.5791, df = 1, p-value = 5.258e-10

table(E2_cor_cast,active_cor_cast)
chisq.test(E2_cor_cast,active_cor_cast)
#X-squared = 1.8325, df = 1, p-value = 0.1758


correlations<-cbind(active_cor_cast,age_cor_cast,claudinlow_cor_cast,tgfb_cor_cast,dtf_cor_cast,twist_cor_cast,E2_cor_cast)
colnames(correlations)<-c("Active","Age","ClaudinLow","TGFB","DTF","TWIST","E2")
write.table(correlations,file="SignatureCorrelationsNormal.txt",sep="\t",col.names=NA)


entrez.fold<-as.vector(unlist(mget(rownames(exprs(sig.fold.rm.eset)),env=hgug4112aENTREZID)))
entrez.full<-as.vector(unlist(mget(rownames(exprs(sig.full.rm.eset)),env=hgug4112aENTREZID)))
missing.fold<-which(is.na(entrez.fold))
missing.full<-which(is.na(entrez.full))

sig.fold.rm.eset<-sig.fold.rm.eset[-missing.fold,] #now 160 instead of 164
sig.full.rm.eset<-sig.full.rm.eset[-missing.full,]

sum(is.na(entrez.fold))
sum(is.na(entrez.full))

#need to collapse probe id based on entrez id (there are duplicate gene probes)

collapsed.sig.fold.exprs<-collapseIDs(exprs(sig.fold.rm.eset),entrez.fold,"mean")
collapsed.sig.full.exprs<-collapseIDs(exprs(sig.full.rm.eset),entrez.full,"mean")

#find the expression of young in the overlapped genes
exprs.age.fold.young<-RM_PAUL_OVERLAP.eset[,sig.fold.rm.eset$age<30] #young defined as < 30 years of age
dim(exprs.age.fold.young) #20

#find the expression of the old
exprs.age.fold.old<-RM_PAUL_OVERLAP.eset[,sig.fold.rm.eset$age>39] #old defined as > 39 years of age
dim(exprs.age.fold.old) #23
medians.exprs.young <- apply(exprs.age.fold.young,1,median)
medians.exprs.old <- apply(exprs.age.fold.old,1,median)

young_signature<-as.matrix(medians.exprs.young)
old_signature<-as.matrix(medians.exprs.old)

#will need to set to -1 if < 0 or 1 if > 0
age_signature<-as.matrix(medians.exprs.young-medians.exprs.old)
age_signature_cast<-as.matrix(apply(age_signature,1,castDirection)) #cast direction is a function MD wrote (example: -.52 -> -1, .52 ->1)

