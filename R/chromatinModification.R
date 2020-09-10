##' 'chromatinModification' compares the chromatin modification abundance between cancer cell lines with highest and lowest dependencies of a gene set (signature).
##'
##'
##' @title chromatinModification
##' @param signature.name Names of a signature (format: character)
##' @param signature Gene names of a signature (format: vector)
##' @param cutoff.freq Cutoff for frequency of cancer cell lines for each cancer type, default 10
##' @param cutoff.percentile Cutoff for percentile of cancer cell lines with highest/lowest dependency, default 0.2
##' @param cutoff.pvalue Cutoff for p-value of the T-test results, default 0.05
##' @param cutoff.diff Cutoff for difference of the T-test results, default 0.1
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
##' chromatinModification(signature.name, signature)


## Main
chromatinModification <- function(signature.name, 
                                  signature, 
                                  cutoff.freq        = 10, 
                                  cutoff.percentile  = 0.2, 
                                  cutoff.pvalue      = 0.05, 
                                  cutoff.diff        = 0.1
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

    # Chromatin modification
    chromatin.share = chromatin[rownames(chromatin) %in% rownames(dep.t),,drop=FALSE]
    head(chromatin.share)
    dim(chromatin.share)
    # 446  31

    chromatin.share.signature.high = chromatin.share[rownames(chromatin.share) %in% rownames(dep.t.signature.high),,drop=FALSE]
    chromatin.share.signature.high$signature = paste0(signature.name, ".dep.high")
    head(chromatin.share.signature.high)
    dim(chromatin.share.signature.high)
    # 40 44
    chromatin.share.signature.low = chromatin.share[rownames(chromatin.share) %in% rownames(dep.t.signature.low),,drop=FALSE]
    chromatin.share.signature.low$signature = paste0(signature.name, ".dep.low")
    head(chromatin.share.signature.low)
    dim(chromatin.share.signature.low)
    # 46 44
    chromatin.share.signature.mid = chromatin.share[!rownames(chromatin.share) %in% c(rownames(dep.t.signature.high), rownames(dep.t.signature.low)),,drop=FALSE]
    chromatin.share.signature.mid$signature = paste0(signature.name, ".dep.mid")
    head(chromatin.share.signature.mid)
    dim(chromatin.share.signature.mid)
    # 360  44

    chromatin.share.signature = rbind(chromatin.share.signature.high, chromatin.share.signature.mid, chromatin.share.signature.low)
    head(chromatin.share.signature)
    dim(chromatin.share.signature)
    # 446  44
    # write.csv(chromatin.share.signature, paste0("chromatin/dep_", signature.name, "_hml", cutoff.percentile, "_chromatin.csv"), quote=FALSE)

    chromatin.high = chromatin.share.signature[chromatin.share.signature$signature %like% "high",,drop=FALSE]
    head(chromatin.high)
    dim(chromatin.high)
    # 40 44
    chromatin.mid = chromatin.share.signature[chromatin.share.signature$signature %like% "mid",,drop=FALSE]
    chromatin.low = chromatin.share.signature[chromatin.share.signature$signature %like% "low",,drop=FALSE]

    # signature.high
    chromatin.high.p = apply(chromatin.share.signature[2:(ncol(chromatin.share.signature)-1)], 2, function(x) signif(t.test(x[1:nrow(chromatin.high)], x[(nrow(chromatin.high)+1):nrow(chromatin.share.signature)])$p.value,5))
    head(chromatin.high.p)
    length(chromatin.high.p)
    # 42
    chromatin.high.q = signif(p.adjust(chromatin.high.p, "BH"),5)
    head(chromatin.high.q)
    length(chromatin.high.q)
    # 42
    chromatin.high.mean1 = apply(chromatin.share.signature[2:(ncol(chromatin.share.signature)-1)],  2, function(x) signif(mean(na.omit(x[1:nrow(chromatin.high)])), 5))
    chromatin.high.mean2 = apply(chromatin.share.signature[2:(ncol(chromatin.share.signature)-1)],  2, function(x) signif(mean(na.omit(x[(nrow(chromatin.high)+1):nrow(chromatin.share.signature)])), 5))
    chromatin.high.diff = apply(chromatin.share.signature[2:(ncol(chromatin.share.signature)-1)], 2, function(x) signif(mean(na.omit(x[1:nrow(chromatin.high)])) - mean(na.omit(x[(nrow(chromatin.high)+1):nrow(chromatin.share.signature)])), 5))
    head(chromatin.high.diff)
    length(chromatin.high.diff)
    # 42

    chromatin.high.deg = data.frame(chromatin.high.mean1, chromatin.high.mean2, chromatin.high.diff, chromatin.high.p, chromatin.high.q, row.names=names(chromatin.high.diff))
    colnames(chromatin.high.deg) = c("chromatin.high", "chromatin.other", "diff.high2other", "pvalue", "qvalue")
    head(chromatin.high.deg)
    dim(chromatin.high.deg)
    # 42  5
    # write.csv(chromatin.high.deg, paste0("chromatin/dep_", signature.name, "_hml", cutoff.percentile, "_chromatin.high.deg.csv"), quote=FALSE)

    # signature.low
    chromatin.low.p = apply(chromatin.share.signature[2:(ncol(chromatin.share.signature)-1)], 2, function(x) signif(t.test(x[1:(nrow(chromatin.share.signature)-nrow(chromatin.low))], x[(nrow(chromatin.share.signature)-nrow(chromatin.low)+1):nrow(chromatin.share.signature)])$p.value,5))
    head(chromatin.low.p)
    length(chromatin.low.p)
    # 42
    chromatin.low.q = signif(p.adjust(chromatin.low.p, "BH"),5)
    head(chromatin.low.q)
    length(chromatin.low.q)
    # 42
    chromatin.low.mean1 = apply(chromatin.share.signature[2:(ncol(chromatin.share.signature)-1)],  2, function(x) signif(mean(na.omit(x[(nrow(chromatin.share.signature)-nrow(chromatin.low)+1):nrow(chromatin.share.signature)])), 5))
    chromatin.low.mean2 = apply(chromatin.share.signature[2:(ncol(chromatin.share.signature)-1)],  2, function(x) signif(mean(na.omit(x[1:(nrow(chromatin.share.signature)-nrow(chromatin.low))])), 5))
    chromatin.low.diff = apply(chromatin.share.signature[2:(ncol(chromatin.share.signature)-1)], 2, function(x) signif(mean(na.omit(x[(nrow(chromatin.share.signature)-nrow(chromatin.low)+1):nrow(chromatin.share.signature)])) - mean(na.omit(x[1:(nrow(chromatin.share.signature)-nrow(chromatin.low))])), 5))
    head(chromatin.low.diff)
    length(chromatin.low.diff)
    # 42

    chromatin.low.deg = data.frame(chromatin.low.mean1, chromatin.low.mean2, chromatin.low.diff, chromatin.low.p, chromatin.low.q, row.names=names(chromatin.low.diff))
    colnames(chromatin.low.deg) = c("chromatin.low", "chromatin.other", "diff.low2other", "pvalue", "qvalue")
    head(chromatin.low.deg)
    dim(chromatin.low.deg)
    # 42  5
    # write.csv(chromatin.low.deg, paste0("chromatin/dep_", signature.name, "_hml", cutoff.percentile, "_chromatin.low.deg.csv"), quote=FALSE)

    # cutoff.pvalue = 1
    # cutoff.diff = 0.1
    set.seed(42)
    p1 <- ggplot(data = chromatin.high.deg, mapping = aes(x = diff.high2other, y = (-1)*log(pvalue, 10))) + 
    geom_point(size=1, color= ifelse(chromatin.high.deg$pvalue < cutoff.pvalue & chromatin.high.deg$diff.high2other > cutoff.diff, "red", ifelse(chromatin.high.deg$pvalue < cutoff.pvalue & chromatin.high.deg$diff.high2other < cutoff.diff*(-1), "blue", "grey60")))+ 
    # xlim(-1,1) + 
    # ylim(0,3) + 
    geom_hline(yintercept = (-1)*log(cutoff.pvalue, 10), linetype="dashed", colour="grey30", size=0.2) +
    geom_vline(xintercept = cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = (-1)*cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    annotate("text", x=min(na.omit(chromatin.high.deg$diff.high2other)), y=max((-1)*log(chromatin.high.deg$pvalue, 10), (-1)*log(chromatin.low.deg$pvalue, 10))*1.1, parse=FALSE, label = paste0("Signature.high cell lines: ", nrow(chromatin.share.signature.high)), color = "red", hjust = 0) + 
    geom_label_repel(aes(label=ifelse(chromatin.high.deg$pvalue < cutoff.pvalue & abs(chromatin.high.deg$diff.high2other) > cutoff.diff, as.character(rownames(chromatin.high.deg)), "")), size = 2, color = ifelse(chromatin.high.deg$diff.high2other > 0, "red", "blue"), segment.size=0.2) +
    labs(x="Modification difference", y="-log10(p value)", title=paste0("Signature [", signature.name, "] High vs. others"))
    p1 <- p1 + theme_classic() + rremove("legend")
    # p1

    p2 <- ggplot(data = chromatin.low.deg, mapping = aes(x = diff.low2other, y = (-1)*log(pvalue, 10))) + 
    geom_point(size=1, color= ifelse(chromatin.low.deg$pvalue < cutoff.pvalue & chromatin.low.deg$diff.low2other > cutoff.diff, "red", ifelse(chromatin.low.deg$pvalue < cutoff.pvalue & chromatin.low.deg$diff.low2other < cutoff.diff*(-1), "blue", "grey60")))+ 
    # xlim(-1,1) + 
    # ylim(0,3) + 
    geom_hline(yintercept = (-1)*log(cutoff.pvalue, 10), linetype="dashed", colour="grey30", size=0.2) +
    geom_vline(xintercept = cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = (-1)*cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    annotate("text", x=min(na.omit(chromatin.low.deg$diff.low2other)), y=max((-1)*log(chromatin.high.deg$pvalue, 10), (-1)*log(chromatin.low.deg$pvalue, 10))*1.1, parse=FALSE, label = paste0("Signature.low cell lines: ", nrow(chromatin.share.signature.low)), color = "red", hjust = 0) + 
    geom_label_repel(aes(label=ifelse(chromatin.low.deg$pvalue < cutoff.pvalue & abs(chromatin.low.deg$diff.low2other) > cutoff.diff, as.character(rownames(chromatin.low.deg)), "")), size = 2, color = ifelse(chromatin.low.deg$diff.low2other > 0, "red", "blue"), segment.size=0.2) +
    labs(x="Modification difference", y="-log10(p value)", title=paste0("Signature [", signature.name, "] Low vs. others"))
    p2 <- p2 + theme_classic() + rremove("legend")
    # p2

    # Arranging the plot using cowplot
    p = suppressWarnings(plot_grid(p1, p2, ncol = 2, align = "hv", rel_widths = c(1,1), rel_heights = c(1,1)))
    # p <- ggplotGrob(p)
    return(p)
}

