##' 'isg' compares the immune signature gene (ISG) between cancer cell lines with highest and lowest dependencies of a gene set (signature).
##'
##'
##' @title isg
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
##' isg(signature.name, signature)


## Main
isg <- function(signature.name, 
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

    # Immune signature gene (ISG)
    ISG.share = ISG[rownames(ISG) %in% rownames(dep.t.signature),,drop=FALSE]
    head(ISG.share)
    dim(ISG.share)
    # 558   9

    ISG.share.signature.high = ISG.share[rownames(ISG.share) %in% rownames(dep.t.signature.high),,drop=FALSE]
    ISG.share.signature.high$signature = paste0(signature.name, ".dep.high")
    head(ISG.share.signature.high)
    dim(ISG.share.signature.high)
    # 56  2
    ISG.share.signature.low = ISG.share[rownames(ISG.share) %in% rownames(dep.t.signature.low),,drop=FALSE]
    ISG.share.signature.low$signature = paste0(signature.name, ".dep.low")
    head(ISG.share.signature.low)
    dim(ISG.share.signature.low)
    # 56  2
    ISG.share.signature.mid = ISG.share[!rownames(ISG.share) %in% c(rownames(dep.t.signature.high), rownames(dep.t.signature.low)),,drop=FALSE]
    ISG.share.signature.mid$signature = paste0(signature.name, ".dep.mid")
    head(ISG.share.signature.mid)
    dim(ISG.share.signature.mid)
    # 442   2

    ISG.share.signature = rbind(ISG.share.signature.high, ISG.share.signature.mid, ISG.share.signature.low)
    head(ISG.share.signature)
    dim(ISG.share.signature)
    # 554   2
    # write.csv(ISG.share.signature, paste0("ISG/dep_", signature.name, "_hml", cutoff.percentile, "_ISG.csv"), quote=TRUE)

    group.high = ISG.share.signature[ISG.share.signature$signature %like% "high",,drop=FALSE]
    group.mid = ISG.share.signature[ISG.share.signature$signature %like% "mid",,drop=FALSE]
    group.low = ISG.share.signature[ISG.share.signature$signature %like% "low",,drop=FALSE]
    # t.test(group1[,1], group2[,1])$p.value
    # 0.001361099
    p.hm = paste("p.hm = ", signif(t.test(group.high$ISG, group.mid$ISG)$p.value, 4))
    p.hl = paste("p.hl = ", signif(t.test(group.high$ISG, group.low$ISG)$p.value, 4))
    p.ml = paste("p.ml = ", signif(t.test(group.mid$ISG, group.low$ISG)$p.value, 4))
    ISG.share.signature$signature = gsub(paste0(signature.name, "\\.dep\\.high"), "1.high", ISG.share.signature$signature)
    ISG.share.signature$signature = gsub(paste0(signature.name, "\\.dep\\.mid"), "2.mid", ISG.share.signature$signature)
    ISG.share.signature$signature = gsub(paste0(signature.name, "\\.dep\\.low"), "3.low", ISG.share.signature$signature)

    p <- ggplot(ISG.share.signature, aes(x=signature, y= ISG, fill=signature)) + 
         geom_violin(trim=FALSE, linetype="blank", na.rm=TRUE)+
         geom_boxplot(width=0.05, fill="white", outlier.size=0.1, na.rm=TRUE)+
         annotate("text", x=0.5, y=max(na.omit(ISG.share.signature$ISG))*1.2, parse=FALSE, hjust=0, label = p.hm)+
         annotate("text", x=0.5, y=max(na.omit(ISG.share.signature$ISG))*1.1, parse=FALSE, hjust=0, label = p.hl)+
         annotate("text", x=0.5, y=max(na.omit(ISG.share.signature$ISG))*1.0, parse=FALSE, hjust=0, label = p.ml)+
         labs(title=paste0("ISG of signature [", signature.name, "]"), x=paste0(signature.name, " score"), y = "ISG") +  theme_classic() + rremove("legend")
    return(p)
}

