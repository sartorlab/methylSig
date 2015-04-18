### R code from vignette source 'methylSig.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: methylSig.Rnw:23-24
###################################################
options(width=90)


###################################################
### code chunk number 2: methylSig.Rnw:69-71
###################################################
library(methylSig)
print(read.table(system.file("extdata", "AML_1.txt", package = "methylSig"),header=T,nrows=6), row.names=F)


###################################################
### code chunk number 3: methylSig.Rnw:80-96
###################################################
fileList = c(system.file("extdata", "AML_1.txt", package = "methylSig"),
             system.file("extdata", "AML_2.txt", package = "methylSig"),
             system.file("extdata", "AML_3.txt", package = "methylSig"),
             system.file("extdata", "AML_4.txt", package = "methylSig"), 
             system.file("extdata", "NBM_1.txt", package = "methylSig"),
             system.file("extdata", "NBM_2.txt", package = "methylSig"),
             system.file("extdata", "NBM_3.txt", package = "methylSig"),
             system.file("extdata", "NBM_4.txt", package = "methylSig"))


sample.id = c("AML1", "AML2", "AML3", "AML4", "NBM1", "NBM2", "NBM3", "NBM4")

treatment = c(1,1,1,1,0,0,0,0)
#### Read Data ####
meth <- methylSigReadData(fileList, sample.ids = sample.id, assembly = "hg18", 
           treatment = treatment, context = "CpG", destranded=TRUE)


###################################################
### code chunk number 4: methylSig.Rnw:101-104
###################################################
meth <- methylSigReadData(fileList, sample.ids = sample.id, assembly = "hg18", 
           treatment = treatment, context = "CpG", destranded=TRUE, 
           minCount=10, maxCount=500, quiet=TRUE)


###################################################
### code chunk number 5: methylSig.Rnw:129-130
###################################################
myDiffSigboth = methylSigCalc(meth, groups=c(1,0), min.per.group=3)


###################################################
### code chunk number 6: methylSig.Rnw:135-137
###################################################
myDiffSigbothDMCs = myDiffSigboth[myDiffSigboth[,"qvalue"] <= 0.05 
                                  & abs(myDiffSigboth[,"meth.diff"])>=25, ] 


###################################################
### code chunk number 7: methylSig.Rnw:144-148
###################################################
### Tiled analysis
methTile = methylSigTile(meth,win.size = 25)

myDiffSigbothTile = methylSigCalc(methTile, groups=c(1,0), min.per.group=3)


###################################################
### code chunk number 8: methylSig.Rnw:155-159
###################################################
###### Variance from sample treatment group "0" only
myDiffSignorm = methylSigCalc(meth, groups=c(1,0), dispersion=0, min.per.group=3)
myDiffSignormTile = methylSigCalc(methTile, groups=c(1,0), 
                                    dispersion=0, min.per.group=3)


###################################################
### code chunk number 9: methylSig.Rnw:166-173
###################################################
###### Variance from both groups and using local information for variance
myDiffSigBothLoc = methylSigCalc(meth, groups=c(1,0), 
          min.per.group=3, local.disp=TRUE, winsize.disp=200)

###### Variance from sample treatment group "0" only and using local information for variance
myDiffSignormLoc = methylSigCalc(meth, groups=c(1,0), dispersion=0, 
          min.per.group=3, local.disp=TRUE, winsize.disp=200)


###################################################
### code chunk number 10: methylSig.Rnw:176-184
###################################################
###### Variance from both groups and using local information for methylation level
myDiffSigBothMLoc = methylSigCalc(meth, groups=c(1,0), 
          min.per.group=3, local.meth=TRUE, winsize.meth=200)

###### Variance from both groups and using local information for methylation level and variance
myDiffSigBothMDLoc = methylSigCalc(meth, groups=c(1,0),  
          min.per.group=3, local.disp=TRUE, winsize.disp=200,
          local.meth=TRUE, winsize.meth=200)


###################################################
### code chunk number 11: methylSig.Rnw:192-199
###################################################
#### Read Data using 5 cores
meth <- methylSigReadData(fileList, sample.ids = sample.id, assembly = "hg18",
           treatment = treatment, context = "CpG", destranded=TRUE,
           num.cores=5, quiet=TRUE)

#### Differential methylation analysis using 5 cores
myDiffSigboth = methylSigCalc(meth, groups=c(1,0), min.per.group=3, num.cores=5)


###################################################
### code chunk number 12: methylSig.Rnw:214-224
###################################################
library("graphics")
cpgInfo = getCpGInfo(system.file("annotation", "cpgi.hg18.bed.txt", 
                                  package = "methylSig"))

myDiffDMCs = myDiffSigboth[myDiffSigboth[,"qvalue"] < 0.05 
                               & abs(myDiffSigboth[,"meth.diff"])>25,]
cpgAnn = cpgAnnotation(cpgInfo,myDiffSigboth)
cpgAnnDmc = cpgAnnotation(cpgInfo, myDiffDMCs)
cpgAnnotationPlot(cpgAnn,main="ALL")
cpgAnnotationPlot(cpgAnnDmc,main="DMCs")


###################################################
### code chunk number 13: CpGALL
###################################################
cpgAnnotationPlot(cpgAnn,main="ALL")


###################################################
### code chunk number 14: CpGDMC
###################################################
cpgAnnotationPlot(cpgAnnDmc,main="DMCs")


###################################################
### code chunk number 15: methylSig.Rnw:258-267
###################################################
refGeneInfo = getRefgeneInfo(system.file("annotation", "refGene.txt",
                                   package = "methylSig"))

refGeneAnn = refGeneAnnotation(refGeneInfo, myDiffSigboth)
refGeneAnnDmc = refGeneAnnotation(refGeneInfo, myDiffDMCs)
refGeneAnnotationPlot(refGeneAnn,main="ALL",
                 priority=c("promoter","cds", "noncoding", "5'utr", "3'utr"))
refGeneAnnotationPlot(refGeneAnnDmc, main="DMC",
                 priority=c("promoter","cds", "noncoding", "5'utr", "3'utr"))


###################################################
### code chunk number 16: RefGeneALL
###################################################
refGeneAnnotationPlot(refGeneAnn,main="ALL",
                 priority=c("promoter","cds", "noncoding", "5'utr", "3'utr"))


###################################################
### code chunk number 17: RefGeneDMC
###################################################
refGeneAnnotationPlot(refGeneAnnDmc, main="DMC",
                 priority=c("promoter","cds", "noncoding", "5'utr", "3'utr"))


###################################################
### code chunk number 18: TFBSEnrichment
###################################################
tfbsInfo = getTFBSInfo(system.file("annotation", "tfbsUniform.txt", 
                                   package = "methylSig"))
DMCIndex = (myDiffSigboth[,"qvalue"] < 0.05 
                 & abs(myDiffSigboth[,"meth.diff"]) > 25)
pvalue = methylSig.tfbsEnrichTest(myDiffSigboth, DMCIndex, tfbsInfo)  


###################################################
### code chunk number 19: methylSig.Rnw:322-324
###################################################
methTileTFs = methylSigTileTFBS(meth, tfbsInfo)
myDiffTFs = methylSigCalc(methTileTFs, groups=c(1,0))


###################################################
### code chunk number 20: DataVisualBroad
###################################################
methylSigPlot(meth, "chr21", c(43000000, 44000000), groups=c(1,0), 
    cpgInfo=cpgInfo,refGeneInfo=refGeneInfo,
    myDiff=myDiffSigboth,tfbsInfo=tfbsInfo,tfbsDense=F,sigQ=0.05)


###################################################
### code chunk number 21: DataVisualNarrow
###################################################
methylSigPlot(meth, "chr21", c(43800000, 43900000), groups=c(1,0), 
    cpgInfo=cpgInfo, refGeneInfo=refGeneInfo,
    myDiff=myDiffSigboth,tfbsInfo=tfbsInfo,tfbsDense=F,sigQ=0.05)


###################################################
### code chunk number 22: methylSig.Rnw:390-391
###################################################
meth1_4 = meth[,1:4]


###################################################
### code chunk number 23: methylSig.Rnw:395-396
###################################################
methSub1_100 = meth[1:100,]


###################################################
### code chunk number 24: methylSig.Rnw:400-402
###################################################
methSubData = meth[1:100,1:2]
methSubData


###################################################
### code chunk number 25: methylSig.Rnw:408-410
###################################################
coverage1 = meth[,"coverage1"]
startTop200 = meth[1:200,"start"]


###################################################
### code chunk number 26: methylSig.Rnw:418-419
###################################################
options(width=80)


###################################################
### code chunk number 27: methylSig.Rnw:421-422
###################################################
myDiffSigboth


###################################################
### code chunk number 28: methylSig.Rnw:429-430
###################################################
myDiff100 = myDiffSigboth[1:100,]


###################################################
### code chunk number 29: methylSig.Rnw:437-438
###################################################
qvalues = myDiffSigboth[,"qvalue"]


###################################################
### code chunk number 30: methylSig.Rnw:445-447
###################################################
myDiffq05D25 = myDiffSigboth[myDiffSigboth[,"qvalue"] < 0.05 
                           & abs(myDiffSigboth[,"meth.diff"]) > 25,]


###################################################
### code chunk number 31: methylSig.Rnw:453-455
###################################################
myDiffp05D25 = myDiffSigboth[myDiffSigboth[,"pvalue"] < 0.05
                           & abs(myDiffSigboth[,"meth.diff"]) > 25,]


###################################################
### code chunk number 32: methylSig.Rnw:463-469
###################################################
methRaw = methylSigReadData(fileList, sample.ids = sample.id,assembly = "hg18", 
             treatment = treatment, context = "CpG", minCount = 0,
             maxCount=Inf, destranded=F, quiet=T)

summary(methRaw[,"numCs1"]/methRaw[,"coverage1"])
summary(methRaw[,"coverage1"])


###################################################
### code chunk number 33: histogram1
###################################################
hist(methRaw[,"numCs1"]/methRaw[,"coverage1"], 
                 main="Histogram of methylation rate for sample 1",
                 xlab="methylation rate")


###################################################
### code chunk number 34: histogram2
###################################################
hist(methRaw[,"coverage1"], main="Histogram of coverage for sample 1",
                 xlab="coverage")


###################################################
### code chunk number 35: methylSig.Rnw:486-496
###################################################
library(gplots)

x = meth[,"numCs"]/meth[, "coverage"]
colnames(x) = meth@sample.ids
rownames(x) =rep(NA, NROW(x))
 
corrALL = cor(x, use="pairwise.complete.obs")
heatmap.2(1-corrALL, na.rm=T, breaks=100, 
          hclustfun = function(x) hclust(x,method="ward"), 
          col="bluered",trace="none", symm=T, keysize=1,density.info="none")


###################################################
### code chunk number 36: methylSig.Rnw:501-509
###################################################
myDiffDMC = myDiffSigboth[myDiffSigboth[,"qvalue"] < 0.05 
                          & abs(myDiffSigboth[,"meth.diff"]) >=25,]
listInMeth = match(myDiffDMC@data.ids, meth@data.ids)
y = x[listInMeth,]
corrDMC = cor(y, use="pairwise.complete.obs")
heatmap.2(1-corrDMC, na.rm=T, breaks=100,
          hclustfun = function(x) hclust(x,method="ward"),
          col="bluered",trace="none", symm=T, keysize=1,density.info="none")


###################################################
### code chunk number 37: heatmap1
###################################################
heatmap.2(1-corrALL, na.rm=T, breaks=100, 
          hclustfun = function(x) hclust(x,method="ward"), 
          col="bluered",trace="none", symm=T, keysize=1,density.info="none")


###################################################
### code chunk number 38: heatmap2
###################################################
heatmap.2(1-corrDMC, na.rm=T, breaks=100,
          hclustfun = function(x) hclust(x,method="ward"),
          col="bluered",trace="none", symm=T, keysize=1,density.info="none")


###################################################
### code chunk number 39: methylSig.Rnw:539-541
###################################################
system("sed -i 's/^+/ /g' methylSig.tex")
system("sed -i 's#ccmb/CoreBA2/sartorlab/yongpark#...#g' methylSig.tex")


