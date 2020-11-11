##' 'hallmark' compares the Hallmark signatures between cancer cell lines with highest and lowest dependencies of a gene set (signature).
##'
##'
##' @title hallmark
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
##' hallmark(signature.name, signature)


## Main
hallmark <- function(signature.name, 
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

    # Hallmark
    hallmark.share.signature.high = hallmark.share[rownames(hallmark.share) %in% rownames(dep.t.signature.high),,drop=FALSE]
    hallmark.share.signature.high$signature = paste0(signature.name, ".dep.high")
    head(hallmark.share.signature.high)
    dim(hallmark.share.signature.high)
    # 40 44
    hallmark.share.signature.low = hallmark.share[rownames(hallmark.share) %in% rownames(dep.t.signature.low),,drop=FALSE]
    hallmark.share.signature.low$signature = paste0(signature.name, ".dep.low")
    head(hallmark.share.signature.low)
    dim(hallmark.share.signature.low)
    # 46 44
    hallmark.share.signature.mid = hallmark.share[!rownames(hallmark.share) %in% c(rownames(dep.t.signature.high), rownames(dep.t.signature.low)),,drop=FALSE]
    hallmark.share.signature.mid$signature = paste0(signature.name, ".dep.mid")
    head(hallmark.share.signature.mid)
    dim(hallmark.share.signature.mid)
    # 360  44

    hallmark.share.signature = rbind(hallmark.share.signature.high, hallmark.share.signature.mid, hallmark.share.signature.low)
    head(hallmark.share.signature)
    dim(hallmark.share.signature)
    # 446  44
    # write.csv(hallmark.share.signature, paste0("hallmark/dep_", signature.name, "_hml", cutoff.percentile, "_hallmark.csv"), quote=FALSE)

    hallmark.high = hallmark.share.signature[hallmark.share.signature$signature %like% "high",,drop=FALSE]
    head(hallmark.high)
    dim(hallmark.high)
    # 40 44
    hallmark.mid = hallmark.share.signature[hallmark.share.signature$signature %like% "mid",,drop=FALSE]
    hallmark.low = hallmark.share.signature[hallmark.share.signature$signature %like% "low",,drop=FALSE]

    # signature.high
    hallmark.high.p = apply(hallmark.share.signature[1:(ncol(hallmark.share.signature)-1)], 2, function(x) signif(t.test(x[1:nrow(hallmark.high)], x[(nrow(hallmark.high)+1):nrow(hallmark.share.signature)])$p.value,5))
    head(hallmark.high.p)
    length(hallmark.high.p)
    # 42
    hallmark.high.q = signif(p.adjust(hallmark.high.p, "BH"),5)
    head(hallmark.high.q)
    length(hallmark.high.q)
    # 42
    hallmark.high.mean1 = apply(hallmark.share.signature[1:(ncol(hallmark.share.signature)-1)],  2, function(x) signif(mean(na.omit(x[1:nrow(hallmark.high)])), 5))
    hallmark.high.mean2 = apply(hallmark.share.signature[1:(ncol(hallmark.share.signature)-1)],  2, function(x) signif(mean(na.omit(x[(nrow(hallmark.high)+1):nrow(hallmark.share.signature)])), 5))
    hallmark.high.diff = apply(hallmark.share.signature[1:(ncol(hallmark.share.signature)-1)], 2, function(x) signif(mean(na.omit(x[1:nrow(hallmark.high)])) - mean(na.omit(x[(nrow(hallmark.high)+1):nrow(hallmark.share.signature)])), 5))
    head(hallmark.high.diff)
    length(hallmark.high.diff)
    # 42

    hallmark.high.deg = data.frame(hallmark.high.mean1, hallmark.high.mean2, hallmark.high.diff, hallmark.high.p, hallmark.high.q, row.names=names(hallmark.high.diff))
    colnames(hallmark.high.deg) = c("hallmark.high", "hallmark.other", "diff.high2other", "pvalue", "qvalue")
    head(hallmark.high.deg)
    dim(hallmark.high.deg)
    # 42  5
    # write.csv(hallmark.high.deg, paste0("hallmark/dep_", signature.name, "_hml", cutoff.percentile, "_hallmark.high.deg.csv"), quote=FALSE)

    # signature.low
    hallmark.low.p = apply(hallmark.share.signature[1:(ncol(hallmark.share.signature)-1)], 2, function(x) signif(t.test(x[1:(nrow(hallmark.share.signature)-nrow(hallmark.low))], x[(nrow(hallmark.share.signature)-nrow(hallmark.low)+1):nrow(hallmark.share.signature)])$p.value,5))
    head(hallmark.low.p)
    length(hallmark.low.p)
    # 42
    hallmark.low.q = signif(p.adjust(hallmark.low.p, "BH"),5)
    head(hallmark.low.q)
    length(hallmark.low.q)
    # 42
    hallmark.low.mean1 = apply(hallmark.share.signature[1:(ncol(hallmark.share.signature)-1)],  2, function(x) signif(mean(na.omit(x[(nrow(hallmark.share.signature)-nrow(hallmark.low)+1):nrow(hallmark.share.signature)])), 5))
    hallmark.low.mean2 = apply(hallmark.share.signature[1:(ncol(hallmark.share.signature)-1)],  2, function(x) signif(mean(na.omit(x[1:(nrow(hallmark.share.signature)-nrow(hallmark.low))])), 5))
    hallmark.low.diff = apply(hallmark.share.signature[1:(ncol(hallmark.share.signature)-1)], 2, function(x) signif(mean(na.omit(x[(nrow(hallmark.share.signature)-nrow(hallmark.low)+1):nrow(hallmark.share.signature)])) - mean(na.omit(x[1:(nrow(hallmark.share.signature)-nrow(hallmark.low))])), 5))
    head(hallmark.low.diff)
    length(hallmark.low.diff)
    # 42

    hallmark.low.deg = data.frame(hallmark.low.mean1, hallmark.low.mean2, hallmark.low.diff, hallmark.low.p, hallmark.low.q, row.names=names(hallmark.low.diff))
    colnames(hallmark.low.deg) = c("hallmark.low", "hallmark.other", "diff.low2other", "pvalue", "qvalue")
    head(hallmark.low.deg)
    dim(hallmark.low.deg)
    # 42  5
    # write.csv(hallmark.low.deg, paste0("hallmark/dep_", signature.name, "_hml", cutoff.percentile, "_hallmark.low.deg.csv"), quote=FALSE)

    # cutoff.pvalue = 0.05
    # cutoff.diff = 0
    set.seed(42)
    p1 <- ggplot(data = hallmark.high.deg, mapping = aes(x = diff.high2other, y = (-1)*log(pvalue, 10))) + 
    geom_point(size=1, color= ifelse(hallmark.high.deg$pvalue < cutoff.pvalue & hallmark.high.deg$diff.high2other > cutoff.diff, "red", ifelse(hallmark.high.deg$pvalue < cutoff.pvalue & hallmark.high.deg$diff.high2other < cutoff.diff*(-1), "blue", "grey60")))+ 
    # xlim(-0.2,0.2) + 
    # ylim(0,4) + 
    geom_hline(yintercept = (-1)*log(cutoff.pvalue, 10), linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = (-1)*cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    annotate("text", x=min(hallmark.high.deg$diff.high2other), y=max(na.omit((-1)*log(hallmark.high.deg$pvalue, 10)), na.omit((-1)*log(hallmark.low.deg$pvalue, 10)))*1.1, parse=FALSE, label = paste0("Signature.high cell lines: ", nrow(hallmark.share.signature.high)), color = "red", hjust = 0) + 
    geom_label_repel(aes(label=ifelse(hallmark.high.deg$pvalue < cutoff.pvalue & abs(hallmark.high.deg$diff.high2other) > cutoff.diff, as.character(rownames(hallmark.high.deg)), "")), size = 2, color = ifelse(hallmark.high.deg$diff.high2other > 0, "red", "blue"), segment.size=0.2) +
    labs(x="Hallmark score difference", y="-log10(p value)", title=paste0("Signature [", signature.name, "] High vs. others"))
    p1 <- p1 + theme_classic() + rremove("legend")
    # p1

    p2 <- ggplot(data = hallmark.low.deg, mapping = aes(x = diff.low2other, y = (-1)*log(pvalue, 10))) + 
    geom_point(size=1, color= ifelse(hallmark.low.deg$pvalue < cutoff.pvalue & hallmark.low.deg$diff.low2other > cutoff.diff, "red", ifelse(hallmark.low.deg$pvalue < cutoff.pvalue & hallmark.low.deg$diff.low2other < cutoff.diff*(-1), "blue", "grey60")))+ 
    # xlim(-0.2,0.2) + 
    # ylim(0,4) + 
    geom_hline(yintercept = (-1)*log(cutoff.pvalue, 10), linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = (-1)*cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    annotate("text", x=min(hallmark.low.deg$diff.low2other), y=max(na.omit( (-1)*log(hallmark.high.deg$pvalue, 10)), na.omit((-1)*log(hallmark.low.deg$pvalue, 10)))*1.1, parse=FALSE, label = paste0("Signature.low cell lines: ", nrow(hallmark.share.signature.low)), color = "red", hjust = 0) + 
    geom_label_repel(aes(label=ifelse(hallmark.low.deg$pvalue < cutoff.pvalue & abs(hallmark.low.deg$diff.low2other) > cutoff.diff, as.character(rownames(hallmark.low.deg)), "")), size = 2, color = ifelse(hallmark.low.deg$diff.low2other > 0, "red", "blue"), segment.size=0.2) +
    labs(x="Hallmark score difference", y="-log10(p value)", title=paste0("Signature [", signature.name, "] Low vs. others"))
    p2 <- p2 + theme_classic() + rremove("legend")
    # p2

    # Arranging the plot using cowplot
    p = suppressWarnings(plot_grid(p1, p2, ncol = 2, align = "hv", rel_widths = c(1,1), rel_heights = c(1,1)))
    # p <- ggplotGrob(p)
    return(p)
}

