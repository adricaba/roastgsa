
####### load R packages
library(roastgsa)
library(GSEABenchmarkeR)
library(EnrichmentBrowser)
library(preprocessCore)
library(hgu133plus2.db)

## R version 4.1.3 (2022-03-10)
## EnrichmentBrowser_2.24.2
## GSEABenchmarkeR_1.14.0
## preprocessCore_1.56.0

####### create destination folder
maindir <- "/Volumes/biostats/acaballe/GSEA_significance/results4paper/"
resp <- paste0(maindir, "roastgsa_on_benchmarking_data/")
dir.create(resp)


####### load data
## 1,2,6,7,8,9,10,12,14,15,16,17,20 - only GROUP
## 3,4,5,11,13,18,19, GROUP, BLOCK
geo2keggR <- loadEData("geo2kegg", nr.datasets=42)

####### vector to save results
ares <- list()
RESULTS  <- list()

Ns  <- list()

####### Rotations and cores
nrot <- 2000
mccores <- 10

####### MalaCards scores
mala.kegg <- readRDS("/Volumes/biostats/acaballe/GSEA_significance/gseaBenchmarking/KEGG.rds")
d2d.map <- readDataId2diseaseCodeMap("/Volumes/biostats/acaballe/GSEA_significance/gseaBenchmarking/GseId2Disease.txt")

####### Execute roastgsa for every dataset
for(whdata in 1:42){
    print(whdata)

    ## get data
    geo2kegg <- maPreproc(geo2keggR[whdata])
    geo2kegg <- runDE(geo2kegg, de.method="limma", padj.method="flexible")

    cnames <- colnames(colData(geo2kegg[[1]]))
    cnames <- cnames[length(cnames):2]

    kegg.gs <- getGenesets(org="hsa", db="kegg")
    kegg.gs <- kegg.gs[which(sapply(kegg.gs,length)>=10&sapply(kegg.gs,length)<500)]


    ysel <- assays(geo2kegg[[1]])$expr

    ynorm <- normalize.quantiles(ysel)
    rownames(ynorm) <- rownames(ysel)
    colnames(ynorm) <- colnames(ysel)
    ysel  <- ynorm
    N  <- ncol(ysel)
    Ns[whdata]  <- N

    covar <- data.frame(colData(geo2kegg[[1]])[,cnames])
    colnames(covar) <- cnames
    index <- lapply(kegg.gs,function(x) which(rownames(ysel)%in%x) )
    index  <- index[sapply(index,length)>5]


    design <- model.matrix( as.formula(paste0("~ ", paste0(cnames, collapse = "+"))), data = covar)
    terms <- colnames(design)
    contrast <- which(colnames(design) == terms[length(terms)])


    ares[["absmean"]] <- roastgsa(ysel, form = as.formula(paste0("~ ", paste0(cnames, collapse = "+"))), covar = covar,
                              index = index, nrot = nrot, mccores = mccores, set.statistic = "absmean",
                              self.contained = FALSE)


    ares[["maxmean"]] <- roastgsa(ysel, form = as.formula(paste0("~ ", paste0(cnames, collapse = "+"))), covar = covar,
                  index = index, nrot = nrot, mccores = mccores, set.statistic = "maxmean",
                              self.contained = FALSE)


    ares[["mean"]] <- roastgsa(ysel, form = as.formula(paste0("~ ", paste0(cnames, collapse = "+"))), covar = covar,
                  index = index, nrot = nrot, mccores = mccores,  set.statistic = "mean",
                              self.contained = FALSE)


    ares[["med"]] <- roastgsa(ysel, form = as.formula(paste0("~ ", paste0(cnames, collapse = "+"))), covar = covar,
                              index = index, nrot = nrot, mccores = mccores, set.statistic = "median",
                              self.contained = FALSE)

    ares[["ksmean"]] <- roastgsa(ysel, form = as.formula(paste0("~ ", paste0(cnames, collapse = "+"))), covar = covar,
                 index = index, nrot = nrot, mccores = mccores, set.statistic = "ksmean",
                             self.contained = FALSE)

    ares[["ksmax"]] <- roastgsa(ysel, form = as.formula(paste0("~ ", paste0(cnames, collapse = "+"))), covar = covar,
                 index = index, nrot = nrot, mccores = mccores, set.statistic = "ksmax",
                             self.contained = FALSE)

    ares[["meanrank"]] <- roastgsa(ysel, form = as.formula(paste0("~ ", paste0(cnames, collapse = "+"))), covar = covar,
                  index = index, nrot = nrot, mccores = mccores, set.statistic = "mean.rank",
                              self.contained = FALSE )

    ## measures M1 and M2
    RESULTS[[whdata]]  <- lapply(ares, function(x){
    a1  <- x$res
    a1$rank  <- rank(abs(a1[,which(regexpr("pval",names(x$res))>0)[1]]))
    rownames(a1)  <- substring(rownames(a1), 1, 8)
    M1  <- sum((1-a1[rownames(mala.kegg[[d2d.map[names(geo2kegg)]]]),]$rank/nrow(a1)) *
        mala.kegg[[d2d.map[names(geo2kegg)]]]$REL.SCORE, na.rm=TRUE)
    whe  <- which(a1[rownames(mala.kegg[[d2d.map[names(geo2kegg)]]]),]$rank < 50)
    M2  <-  sum((1-a1[rownames(mala.kegg[[d2d.map[names(geo2kegg)]]]),]$rank/nrow(a1))[whe] *
        mala.kegg[[d2d.map[names(geo2kegg)]]]$REL.SCORE[whe], na.rm=TRUE)
    c(M1,M2)
   })
}

saveRDS(RESULTS, file = paste0(resp,"/RESULTS_dec22.RDS"))




#######  barplots
load(paste0(resp,"/RESULTS_dec22.RDS"))

sum1  <- sapply(1:7, function(a) apply(apply(sapply(RESULTS,function(x) -sapply(x,"[",1)),2,rank),1,function(x) sum(x==a)))
colnames(sum1)  <- paste0("rank_",1:7)
M1  <- apply(apply(sapply(RESULTS,function(x) -sapply(x,"[",1)),2,rank),1,mean)

pdf(paste0(resp,"M1_barplot_noweights.pdf"), width =8, height =8)
layout(cbind(1:8,9:16),widths = c(1,0.1), height = c(rep(1,7),0.5))
par(mar=c(0,1.5,1,0))
for(k in 1:7) {
   bp  <-  barplot(t(sum1)[,k],beside = TRUE,xaxt="n", yaxt = "n",ylim= c(-0.3,14), col =k, ylab ="")
    title(ylab=rownames(sum1)[k], line=0, cex.lab=1.7 )
}
par(mar=c(0,1.5,1,0))
bp  <-  barplot(t(sum1)[,k],beside = TRUE,xaxt="n", yaxt = "n",ylim= c(-0.3,14), ylab ="",col="white", border=NA)
text(bp[,1], rep(6,7),paste0("rank_",1:7), cex = 1.5)

par(mar=c(0,0,0,0))
for(k in 1:7){
    plot(0,0, xaxt="n",yaxt="n",xlab="",ylab="", pch="", bty="n")
    text(0,-0.5,round(M1,2)[k], cex=1.7, col =1)
}
plot(0,0, xaxt="n",yaxt="n",xlab="",ylab="", pch="", bty="n")
text(0,-0.3,"avg", cex=1.7)
dev.off()

sum2  <- sapply(1:7, function(a) apply(apply(sapply(RESULTS,function(x) -sapply(x,"[",2)),2,rank),1,function(x) sum(x==a)))
colnames(sum2)  <- paste0("rank_",1:7)
M2 <- apply(apply(sapply(RESULTS,function(x) -sapply(x,"[",2)),2,rank),1,mean)

pdf(paste0(resp,"M2_barplot_noweights.pdf"), width =8, height =8)
layout(cbind(1:8,9:16),widths = c(1,0.1), height = c(rep(1,7),0.5))
par(mar=c(0,1.5,1,0))
for(k in 1:7) {
   bp  <-  barplot(t(sum2)[,k],beside = TRUE,xaxt="n", yaxt = "n",ylim= c(-0.3,13), col =k, ylab ="")
    title(ylab=rownames(sum2)[k], line=0, cex.lab=1.7 )
}
par(mar=c(0,1.5,1,0))
bp  <-  barplot(t(sum2)[,k],beside = TRUE,xaxt="n", yaxt = "n",ylim= c(-0.3,13), ylab ="",col="white", border=NA)
text(bp[,1], rep(6,7),paste0("rank_",1:7), cex = 1.5)

par(mar=c(0,0,0,0))
for(k in 1:7){
    plot(0,0, xaxt="n",yaxt="n",xlab="",ylab="", pch="", bty="n")
    text(0,-0.5,round(M2,2)[k], cex=1.7, col =1)
}
plot(0,0, xaxt="n",yaxt="n",xlab="",ylab="", pch="", bty="n")
text(0,-0.3,"avg", cex=1.7)
dev.off()


###### heatmaps
library(ComplexHeatmap)
library(circlize)

dad <-(apply(sapply(RESULTS,function(x) -sapply(x,"[",2)),2,rank)) -4
colnames(dad) <- rep("",ncol(dad))

pdf(paste0(resp,"heatmapwithalldata_M2_noweights.pdf"), height =5, width=10)
Heatmap((dad),  col = colorRamp2(c(3, 0, -3), c("blue", "white", "red")),
        heatmap_legend_param = list(title = "centred rank") )
dev.off()


dad <-(apply(sapply(RESULTS,function(x) -sapply(x,"[",1)),2,rank)) -4
colnames(dad) <- rep("",ncol(dad))

pdf(paste0(resp,"heatmapwithalldata_M1_noweights.pdf"), height =5, width=10)
Heatmap((dad),  col = colorRamp2(c(3, 0, -3), c("blue", "white", "red")),
                heatmap_legend_param = list(title = "centred rank") )
dev.off()
