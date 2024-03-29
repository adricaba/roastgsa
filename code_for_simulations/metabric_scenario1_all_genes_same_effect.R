#################################################################################################
##############  Scenario 1. Simulate new data from Metabric cov: same effect all genes  #########
#################################################################################################
## R version 4.1.3 (2022-03-10)

####### load R packages
library(roastgsa)
library(parallel)

## function tab4latex
tab4latex <- function(tab, file, dec = 3)
    {
        cat("\n", file = file)
        cat(paste(c(" ",colnames(tab)), collapse = "&", sep = " "),"\\\\", "\n", file = file, append = TRUE)
        for(i in 1:nrow(tab))
            cat(paste(c(rownames(tab)[i],round(tab[i,],dec)), collapse = "&", sep = " "),"\\\\",
                "\n", file = file, append = TRUE)
    }

####### set directory
maindir <- getwd()

resp <- paste0(maindir,"/scenarios_recoveryrates_p/")
dir.create(resp)


####### Rotations and cores
Nite  <- 1000
mc.cores  <- 20
nrot  <- 200


####### load correlation from metabric and hallmarks (v3) index
load(paste0(resp,"cor_hallmarks_metabric.Rdata"))
eS <- eigen(cormat + diag(rep(1,ncol(cormat))), symmetric = TRUE)
ev <- eS$values
p <- dim(cormat)[1]
tom <- eS$vectors %*% diag(sqrt(pmax(ev, 0)), p)
tcov <- cormat + diag(rep(1,ncol(cormat)))

####### Number of samples and testing sets
Nselall <- c(6, 10, 20, 30)
whindall <-  c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_ANGIOGENESIS")


####### Simulations
k <- 1
avgcor2t <- numeric()
power1  <- list()

for(SA in c(1,2,3)){
  whind <- whindall[SA]
  aux = cov2cor(tcov[hallmarks_index[[whind]],hallmarks_index[[whind]]])
  avgcor2t[SA] <- ifelse(SA==3,0, mean(aux[lower.tri(aux)]))
  power1[[whind]]  <- list()
  for(k in 1:length(Nselall)){
    Nsel <- Nselall[k]
    voi <- c(rep(0,Nsel/2),rep(1,Nsel/2))
    print(k)
    index2 <- hallmarks_index

    ale <- mclapply(1:Nite, function(ite){
        set.seed(213 + ite)
        mu <- array(0, dim = c(p,Nsel))
        mu[,(Nsel/2+1):Nsel] <- rnorm(p, 0.2, 0.1)  #
        mu[index2[[whind]],(Nsel/2+1):Nsel] <- rnorm(length(index2[[whind]]), 0.5, 0.1)
        print(ite)
        ysim  <- matrix(rnorm(p * Nsel), Nsel)
        if(SA==3)
            ysim  <- drop(mu) + t(ysim)
        else
            ysim  <- drop(mu) + tom %*%  t(ysim)
        rownames(ysim) <- rownames(cormat)

        des <- model.matrix(as.formula(paste0("~ ivar")), data=data.frame(ivar = as.factor(voi)))
        lmf1  <- lmFit(ysim, des)
        y2 <- ysim  - as.matrix(lmf1$coef ) %*% t(as.matrix((des) ))


        N  <- ncol(y2)

        SC.median  <- roastgsa((ysim), form = as.formula(paste0("~ ivar")), covar = data.frame(ivar = as.factor(voi)),
                               index = index2, nrot = nrot, mccores = 1, self.contained = TRUE,
                               set.statistic = "median")$res[names(index2),]$pval

        CO.median <- roastgsa((ysim), form = as.formula(paste0("~ ivar")), covar = data.frame(ivar = as.factor(voi)),
                              index = index2, nrot = nrot, mccores = 1, self.contained = FALSE,
                              set.statistic = "median")$res[names(index2),]$pval


        SC.absmean  <- roastgsa((ysim), form = as.formula(paste0("~ ivar")), covar = data.frame(ivar = as.factor(voi)),
                                index = index2, nrot = nrot, mccores = 1, self.contained = TRUE,
                                set.statistic = "absmean")$res[names(index2),]$pval

        CO.absmean  <- roastgsa((ysim), form = as.formula(paste0("~ ivar")), covar = data.frame(ivar = as.factor(voi)),
                                index = index2, nrot = nrot, mccores = 1, self.contained = FALSE,
                                set.statistic = "absmean")$res[names(index2),]$pval

        CO.mean.rank  <- roastgsa((ysim), form = as.formula(paste0("~ ivar")), covar = data.frame(ivar = as.factor(voi)),
                                  index = index2, nrot = nrot, mccores = 1, self.contained = FALSE,
                                  set.statistic = "mean.rank")$res[names(index2),]$pval.diff

        (SC.maxmean  <- roastgsa((ysim), form = as.formula(paste0("~ ivar")), covar = data.frame(ivar = as.factor(voi)),
                                index = index2, nrot = nrot, mccores = 1, self.contained = TRUE,
                                set.statistic = "maxmean")$res[names(index2),]$pval)


        CO.maxmean <- roastgsa((ysim), form = as.formula(paste0("~ ivar")), covar = data.frame(ivar = as.factor(voi)),
                               index = index2, nrot = nrot, mccores = 1, self.contained = FALSE,
                               set.statistic = "maxmean")$res[names(index2),]$pval

        CO.ksmean <- roastgsa((ysim), form = as.formula(paste0("~ ivar")), covar = data.frame(ivar = as.factor(voi)),
                            index = index2, nrot = nrot, mccores = 1, self.contained = FALSE,
                            set.statistic = "ksmean")$res[names(index2),]$pval

        CO.ksmax  <- roastgsa((ysim), form = as.formula(paste0("~ ivar")), covar = data.frame(ivar = as.factor(voi)),
                             index = index2, nrot = nrot, mccores = 1, self.contained = FALSE,
                             set.statistic = "ksmax")$res[names(index2),]$pval

        SC.mean <- roastgsa((ysim), form = as.formula(paste0("~ ivar")), covar = data.frame(ivar = as.factor(voi)),
                            index = index2, nrot = nrot, mccores = 1, self.contained = TRUE,
                            set.statistic = "mean")$res[names(index2),]$pval
        CO.mean  <- roastgsa((ysim), form = as.formula(paste0("~ ivar")), covar = data.frame(ivar = as.factor(voi)),
                             index = index2, nrot = nrot, mccores = 1, self.contained = FALSE,
                             set.statistic = "mean")$res[names(index2),]$pval

        list(SC.mean, SC.maxmean, SC.median, SC.absmean,
             CO.mean, CO.maxmean, CO.median, CO.absmean, CO.mean.rank, CO.ksmean, CO.ksmax)

    }, mc.cores = mc.cores)
    power1[[whind]][[k]] <- lapply(1:length(ale[[1]]),function(o) sapply(ale, function(x) x[[o]]))
   save(power1, avgcor2t, file= paste0(resp, "scenario1_same_effect_allgenes.Rdata"))
  }
}


##########
load(paste0(resp, "scenario1_same_effect_allgenes.Rdata"))
nm  <- c("SC.mean", "SC.maxmean",  "SC.median", "SC.absmean",
         "CO.mean", "CO.maxmean", "CO.median", "CO.absmean", "CO.mean.rank", "CO.ksmean", "CO.ksmax")

cat("\n", file = paste0(resp, "scenario1_table.txt"), append = FALSE)

for(K in c(1,2,3)){
    print(K)
    L <- which(names(index2) == whindall[K])
    power2 <- power1[[whindall[K]]]
    avg.tab <-  do.call(cbind,lapply(1:4,function(o) sapply(1:11, function(jk)
        mean(power2[[o]][[jk]][L,]<0.05))))
    rownames(avg.tab) <- nm
    colnames(avg.tab) <- paste0("N.",Nselall)
    tab1  <- avg.tab[nm,] *1000
    rownames(tab1)  <- nm

    for(o in 1:nrow(tab1))
    {
        if(o==5) cat(paste0("A",K," & ",rownames(tab1)[o], "& ", paste0(tab1[o,], collapse =" & "), '\\','\\'), "\n", file = paste0(resp, "scenario1_table.txt"), append = TRUE)
        else cat(paste0("& ",rownames(tab1)[o], "& ", paste0(tab1[o,], collapse =" & "), '\\','\\'), "\n", file = paste0(resp, "scenario1_table.txt"), append = TRUE)
    }
    cat("\\hline ", "\n", file = paste0(resp, "scenario1_table.txt"), append = TRUE)
}

