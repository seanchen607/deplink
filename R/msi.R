##' 'msi' compares the microsatellite instability (MSI) between cancer cell lines with highest and lowest dependencies of a gene set (signature).
##'
##'
##' @title msi
##' @param signature.name Names of a signature (format: character)
##' @param signature Gene names of a signature (format: vector)
##' @param cutoff.freq Cutoff for frequency of cancer cell lines for each cancer type, default 10
##' @param cutoff.percentile Cutoff for percentile of cancer cell lines with highest/lowest dependency, default 0.2
##' @return plot
##' @importFrom stats complete.cases
##' @importFrom wesanderson wes_palette
##' @importFrom cowplot plot_grid
##' @importFrom purrr map
##' @importFrom ggrepel geom_label_repel
##' @import data.table ggpubr ggplot2 
##' @export 
##' @author Xiao Chen
##' @references 1. X Chen, J McGuire, F Zhu, X Xu, Y Li, D Karagiannis, R Dalla-Favera, A Ciccia, J Amengual, C Lu (2020). 
##' Harnessing genetic dependency correlation network to reveal chromatin vulnerability in cancer.  
##' In preparation. 
##' @examples
##' source(system.file("script", "load_libs.R", package = "deplink"))
##' signature.name = "9-1-1"
##' signature = c("RAD9A", "RAD1", "HUS1", "RAD17")
##' msi(signature.name, signature)


## Main
msi <- function(signature.name, 
                signature, 
                cutoff.freq        = 10, 
                cutoff.percentile  = 0.2 
                ) {

    # Primary.Disease.freq = table(dep.t.signature.meta.order$disease)
    Primary.Disease.freq = table(dep.t.meta$tcga_code)
    Primary.Disease.freq.cutoff = Primary.Disease.freq[Primary.Disease.freq >= cutoff.freq]
    TCGA.tumor.target = TCGA.tumor[TCGA.tumor$tcga_code %in% names(Primary.Disease.freq.cutoff),,drop=FALSE]
    head(TCGA.tumor.target)
    dim(TCGA.tumor.target)
    # 1047    1

    dep.t = dep.t[rownames(dep.t) %in% rownames(TCGA.tumor.target),,drop=FALSE]
    head(dep.t)
    dim(dep.t)
    # 458 18333

    dep.t.signature = dep.t[,colnames(dep.t) %in% signature, drop=FALSE]
    dep.t.signature$signature.score = rowMeans(dep.t.signature)*(-1)
    dep.t.signature = dep.t.signature[order(dep.t.signature$signature.score, decreasing=TRUE),]
    head(dep.t.signature)
    dim(dep.t.signature)
    # 558  5

    dep.t.signature.high = dep.t.signature[1:ceiling(nrow(dep.t.signature)*cutoff.percentile),,drop=FALSE]
    head(dep.t.signature.high)
    dim(dep.t.signature.high)
    # 56 12
    dep.t.signature.low = dep.t.signature[(nrow(dep.t.signature)-ceiling(nrow(dep.t.signature)*cutoff.percentile) + 1):nrow(dep.t.signature),,drop=FALSE]
    head(dep.t.signature.low)
    dim(dep.t.signature.low)
    # 56 12
    # write.csv(dep.t.signature, paste0("dep_", signature.name, "_score.csv"))
    # write.csv(dep.t.signature.high, paste0("dep_", signature.name, "_score.high", cutoff.percentile, ".csv"))
    # write.csv(dep.t.signature.low, paste0("dep_", signature.name, "_score.low", cutoff.percentile, ".csv"))
    #############################################

    # Microsatellite instability (MSI)
    MSI.share = MSI[rownames(MSI) %in% rownames(dep.t.signature),,drop=FALSE]
    head(MSI.share)
    dim(MSI.share)
    # 474  11

    MSI.share.signature.high = MSI.share[rownames(MSI.share) %in% rownames(dep.t.signature.high),,drop=FALSE]
    MSI.share.signature.high$signature = paste0(signature.name, ".dep.high")
    head(MSI.share.signature.high)
    dim(MSI.share.signature.high)
    # 56  2
    MSI.share.signature.low = MSI.share[rownames(MSI.share) %in% rownames(dep.t.signature.low),,drop=FALSE]
    MSI.share.signature.low$signature = paste0(signature.name, ".dep.low")
    head(MSI.share.signature.low)
    dim(MSI.share.signature.low)
    # 56  2
    MSI.share.signature.mid = MSI.share[!rownames(MSI.share) %in% c(rownames(dep.t.signature.high), rownames(dep.t.signature.low)),,drop=FALSE]
    MSI.share.signature.mid$signature = paste0(signature.name, ".dep.mid")
    head(MSI.share.signature.mid)
    dim(MSI.share.signature.mid)
    # 442   2

    MSI.share.signature = rbind(MSI.share.signature.high, MSI.share.signature.mid, MSI.share.signature.low)
    head(MSI.share.signature)
    dim(MSI.share.signature)
    # 554   2
    # write.csv(MSI.share.signature, paste0("genome.instability/MSI/dep_", signature.name, "_hml", cutoff.percentile, "_MSI.csv"), quote=TRUE)
    # write.csv(MSI.share.signature, paste0("genome.instability/MSI/dep_", signature.name, "_hml", cutoff.percentile, "_MSI.GDSC.csv"), quote=TRUE)

    group.high = MSI.share.signature[MSI.share.signature$signature %like% "high",,drop=FALSE]
    group.mid = MSI.share.signature[MSI.share.signature$signature %like% "mid",,drop=FALSE]
    group.low = MSI.share.signature[MSI.share.signature$signature %like% "low",,drop=FALSE]
    # t.test(group1[,1], group2[,1])$p.value
    # 0.001361099
    p.hm = paste("p.hm = ", signif(t.test(group.high$MSI, group.mid$MSI)$p.value, 4))
    p.hl = paste("p.hl = ", signif(t.test(group.high$MSI, group.low$MSI)$p.value, 4))
    p.ml = paste("p.ml = ", signif(t.test(group.mid$MSI, group.low$MSI)$p.value, 4))
    MSI.share.signature$signature = gsub(paste0(signature.name, "\\.dep\\.high"), "1.high", MSI.share.signature$signature)
    MSI.share.signature$signature = gsub(paste0(signature.name, "\\.dep\\.mid"), "2.mid", MSI.share.signature$signature)
    MSI.share.signature$signature = gsub(paste0(signature.name, "\\.dep\\.low"), "3.low", MSI.share.signature$signature)

    p <- ggplot(MSI.share.signature, aes(x=signature, y= log10(MSI), fill=signature)) + 
         geom_violin(trim=FALSE, linetype="blank", na.rm=TRUE)+
         geom_boxplot(width=0.05, fill="white", outlier.size=0.1, na.rm=TRUE)+
         annotate("text", x=0.5, y=log10(max(na.omit(MSI.share.signature$MSI))*1.5), parse=FALSE, hjust=0, label = p.hm)+
         annotate("text", x=0.5, y=log10(max(na.omit(MSI.share.signature$MSI))*1.1), parse=FALSE, hjust=0, label = p.hl)+
         annotate("text", x=0.5, y=log10(max(na.omit(MSI.share.signature$MSI))*0.8), parse=FALSE, hjust=0, label = p.ml)+
         labs(title=paste0("MSI of signature [", signature.name, "]"), x=paste0(signature.name, " score"), y = "MSI (log10)") +  theme_classic() + rremove("legend")
    return(p)
}

