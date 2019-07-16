# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Sat Oct 6 05:36:05 EDT 2018

################################################################
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE6919", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL8300", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- paste0("XXXXXXXXXXXXXXXXXX0000000000000000000000000XXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXX1111111111111111111111111111111111111111111111",
               "11111111111111111XXXX")
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
cat(sml)
gset <- gset[ ,sel]

# log2 transform
ex <- exprs(gset)

# d = t(ex["36553_at",])
# d = (d - min(d))/(max(d)-min(d))
# e = t(ex["38355_at",])
# e = (e - min(e))/(max(e)-min(e))
# f = t(ex["41138_at",])
# f = (f - min(f))/(max(f)-min(f))
# g = t(ex["38182_at",])
# g = (g - min(g))/(max(g)-min(g))
# mat = rbind(d,e,f,g)
# rownames(mat) <- c("ASMTL","DDX3Y","CD99","TXLNGY")
# heatmap(mat)

qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=12625)
tT <- subset(tT, select=c("ID","t","B","logFC","Gene.symbol"))

tT[,"Mag. logFC"] <- abs(tT[,"logFC"])
tT <- tT[order(-tT$`Mag. logFC`),]
top500 <- tT[which(tT$logFC<0),]
top500 <- top500[1:500,]
tT <- tT[order(tT$`Mag. logFC`),]
bot500 <- tT[which(tT$logFC<0),]
bot500 <- bot500[1:500,]
top500[,"Class"] <- 1
bot500[,"Class"] <- 0

table <- rbind(top500,bot500)

# write.table(tT, file=stdout(), row.names=F, sep="\t")


################################################################
#   Boxplot for selected GEO samples
library(Biobase)
library(GEOquery)

# load series and platform data from GEO

gset <- getGEO("GSE6919", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL8300", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# group names for all samples in a series
gsms <- paste0("XXXXXXXXXXXXXXXXXX0000000000000000000000000XXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXX1111111111111111111111111111111111111111111111",
               "11111111111111111XXXX")
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
sml <- paste("G", sml, sep="")  #set group names

# eliminate samples marked as "X"
sel <- which(sml != "GX")
sml <- sml[sel]
gset <- gset[ ,sel]

# order samples by group
ex <- exprs(gset)[ , order(sml)]
sml <- sml[order(sml)]
fl <- as.factor(sml)
labels <- c("Metastatic","Adjacent")

# set parameters and draw the plot
palette(c("#f4dfdf","#dfeaf4", "#AABBCC"))
dev.new(width=4+dim(gset)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
title <- paste ("GSE6919", '/', annotation(gset), " selected samples", sep ='')
boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
legend("topleft", labels, fill=palette(), bty="n")

# extracting top 250 results and saving to a file
ids <- tT[,1,drop = FALSE]
ids <- as.matrix(ids)
nor_met <- ex[ids, ]
nor_met <- as.matrix(nor_met)
WriteMatrixToFile(nor_met, "D:/Academics/Research/JBCB/Tool/Gene-Analysis-Tool/bioinformatics-software-api/Nor_Met.txt", blnRowNames = TRUE, blnColNames = TRUE )

ex_trainer <- ex[table[,1],]

# for (i in 1:nrow(table)) 
# {
#   print(i)
#   if (table[i,"Mag. logFC"] > 2) 
#   {
#     tT[i,"Class"] <- 1
#   }
# }

Class <- table[,"Class"]
ex_trainer <- cbind(ex_trainer, Class)
ex_trainer <- as.data.frame(ex_trainer)