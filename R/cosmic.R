##' 'cosmic' compares the COSMIC signatures between cancer cell lines with highest and lowest dependencies of a gene set (signature).
##'
##'
##' @title cosmic
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
##' cosmic(signature.name, signature)


## Main
cosmic <- function(signature.name, 
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

    # Mutation - COSMIC
    cosmic.share.signature.high = cosmic.share[rownames(cosmic.share) %in% rownames(dep.t.signature.high),,drop=FALSE]
    cosmic.share.signature.high$signature = paste0(signature.name, ".dep.high")
    head(cosmic.share.signature.high)
    dim(cosmic.share.signature.high)
    # 40 44
    cosmic.share.signature.low = cosmic.share[rownames(cosmic.share) %in% rownames(dep.t.signature.low),,drop=FALSE]
    cosmic.share.signature.low$signature = paste0(signature.name, ".dep.low")
    head(cosmic.share.signature.low)
    dim(cosmic.share.signature.low)
    # 46 44
    cosmic.share.signature.mid = cosmic.share[!rownames(cosmic.share) %in% c(rownames(dep.t.signature.high), rownames(dep.t.signature.low)),,drop=FALSE]
    cosmic.share.signature.mid$signature = paste0(signature.name, ".dep.mid")
    head(cosmic.share.signature.mid)
    dim(cosmic.share.signature.mid)
    # 360  44

    cosmic.share.signature = rbind(cosmic.share.signature.high, cosmic.share.signature.mid, cosmic.share.signature.low)
    head(cosmic.share.signature)
    dim(cosmic.share.signature)
    # 446  44
    # write.csv(cosmic.share.signature, paste0("cosmic/dep_", signature.name, "_hml", cutoff.percentile, "_cosmic.csv"), quote=FALSE)

    cosmic.high = cosmic.share.signature[cosmic.share.signature$signature %like% "high",,drop=FALSE]
    head(cosmic.high)
    dim(cosmic.high)
    # 40 44
    cosmic.mid = cosmic.share.signature[cosmic.share.signature$signature %like% "mid",,drop=FALSE]
    cosmic.low = cosmic.share.signature[cosmic.share.signature$signature %like% "low",,drop=FALSE]

    # signature.high
    cosmic.high.p = apply(cosmic.share.signature[1:(ncol(cosmic.share.signature)-1)], 2, function(x) signif(t.test(x[1:nrow(cosmic.high)], x[(nrow(cosmic.high)+1):nrow(cosmic.share.signature)])$p.value,5))
    head(cosmic.high.p)
    length(cosmic.high.p)
    # 42
    cosmic.high.q = signif(p.adjust(cosmic.high.p, "BH"),5)
    head(cosmic.high.q)
    length(cosmic.high.q)
    # 42
    cosmic.high.mean1 = apply(cosmic.share.signature[1:(ncol(cosmic.share.signature)-1)],  2, function(x) signif(mean(na.omit(x[1:nrow(cosmic.high)])), 5))
    cosmic.high.mean2 = apply(cosmic.share.signature[1:(ncol(cosmic.share.signature)-1)],  2, function(x) signif(mean(na.omit(x[(nrow(cosmic.high)+1):nrow(cosmic.share.signature)])), 5))
    cosmic.high.diff = apply(cosmic.share.signature[1:(ncol(cosmic.share.signature)-1)], 2, function(x) signif(mean(na.omit(x[1:nrow(cosmic.high)])) - mean(na.omit(x[(nrow(cosmic.high)+1):nrow(cosmic.share.signature)])), 5))
    head(cosmic.high.diff)
    length(cosmic.high.diff)
    # 42

    cosmic.high.deg = data.frame(cosmic.high.mean1, cosmic.high.mean2, cosmic.high.diff, cosmic.high.p, cosmic.high.q, row.names=names(cosmic.high.diff))
    colnames(cosmic.high.deg) = c("cosmic.high", "cosmic.other", "diff.high2other", "pvalue", "qvalue")
    head(cosmic.high.deg)
    dim(cosmic.high.deg)
    # 42  5
    # write.csv(cosmic.high.deg, paste0("cosmic/dep_", signature.name, "_hml", cutoff.percentile, "_cosmic.high.deg.csv"), quote=FALSE)

    # signature.low
    cosmic.low.p = apply(cosmic.share.signature[1:(ncol(cosmic.share.signature)-1)], 2, function(x) signif(t.test(x[1:(nrow(cosmic.share.signature)-nrow(cosmic.low))], x[(nrow(cosmic.share.signature)-nrow(cosmic.low)+1):nrow(cosmic.share.signature)])$p.value,5))
    head(cosmic.low.p)
    length(cosmic.low.p)
    # 42
    cosmic.low.q = signif(p.adjust(cosmic.low.p, "BH"),5)
    head(cosmic.low.q)
    length(cosmic.low.q)
    # 42
    cosmic.low.mean1 = apply(cosmic.share.signature[1:(ncol(cosmic.share.signature)-1)],  2, function(x) signif(mean(na.omit(x[(nrow(cosmic.share.signature)-nrow(cosmic.low)+1):nrow(cosmic.share.signature)])), 5))
    cosmic.low.mean2 = apply(cosmic.share.signature[1:(ncol(cosmic.share.signature)-1)],  2, function(x) signif(mean(na.omit(x[1:(nrow(cosmic.share.signature)-nrow(cosmic.low))])), 5))
    cosmic.low.diff = apply(cosmic.share.signature[1:(ncol(cosmic.share.signature)-1)], 2, function(x) signif(mean(na.omit(x[(nrow(cosmic.share.signature)-nrow(cosmic.low)+1):nrow(cosmic.share.signature)])) - mean(na.omit(x[1:(nrow(cosmic.share.signature)-nrow(cosmic.low))])), 5))
    head(cosmic.low.diff)
    length(cosmic.low.diff)
    # 42

    cosmic.low.deg = data.frame(cosmic.low.mean1, cosmic.low.mean2, cosmic.low.diff, cosmic.low.p, cosmic.low.q, row.names=names(cosmic.low.diff))
    colnames(cosmic.low.deg) = c("cosmic.low", "cosmic.other", "diff.low2other", "pvalue", "qvalue")
    head(cosmic.low.deg)
    dim(cosmic.low.deg)
    # 42  5
    # write.csv(cosmic.low.deg, paste0("cosmic/dep_", signature.name, "_hml", cutoff.percentile, "_cosmic.low.deg.csv"), quote=FALSE)

    # cutoff.pvalue = 0.05
    # cutoff.diff = 0
    set.seed(42)
    p1 <- ggplot(data = cosmic.high.deg, mapping = aes(x = diff.high2other, y = (-1)*log(pvalue, 10))) + 
    geom_point(size=1, color= ifelse(cosmic.high.deg$pvalue < cutoff.pvalue & cosmic.high.deg$diff.high2other > 100*cutoff.diff, "red", ifelse(cosmic.high.deg$pvalue < cutoff.pvalue & cosmic.high.deg$diff.high2other < 100*cutoff.diff*(-1), "blue", "grey60")))+ 
    # xlim(-0.2,0.2) + 
    # ylim(0,4) + 
    geom_hline(yintercept = (-1)*log(cutoff.pvalue, 10), linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = 100*cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = (-1)*100*cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    annotate("text", x=min(cosmic.high.deg$diff.high2other), y=max(na.omit((-1)*log(cosmic.high.deg$pvalue, 10)), na.omit((-1)*log(cosmic.low.deg$pvalue, 10)))*1.1, parse=FALSE, label = paste0("Signature.high cell lines: ", nrow(cosmic.share.signature.high)), color = "red", hjust = 0) + 
    geom_label_repel(aes(label=ifelse(cosmic.high.deg$pvalue < cutoff.pvalue & abs(cosmic.high.deg$diff.high2other) > 100*cutoff.diff, as.character(rownames(cosmic.high.deg)), "")), size = 2, color = ifelse(cosmic.high.deg$diff.high2other > 0, "red", "blue"), segment.size=0.2) +
    labs(x="COSMIC difference", y="-log10(p value)", title=paste0("Signature [", signature.name, "] High vs. others"))
    p1 <- p1 + theme_classic() + rremove("legend")
    # p1

    p2 <- ggplot(data = cosmic.low.deg, mapping = aes(x = diff.low2other, y = (-1)*log(pvalue, 10))) + 
    geom_point(size=1, color= ifelse(cosmic.low.deg$pvalue < cutoff.pvalue & cosmic.low.deg$diff.low2other > 100*cutoff.diff, "red", ifelse(cosmic.low.deg$pvalue < cutoff.pvalue & cosmic.low.deg$diff.low2other < 100*cutoff.diff*(-1), "blue", "grey60")))+ 
    # xlim(-0.2,0.2) + 
    # ylim(0,4) + 
    geom_hline(yintercept = (-1)*log(cutoff.pvalue, 10), linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = 100*cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = (-1)*100*cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    annotate("text", x=min(cosmic.low.deg$diff.low2other), y=max(na.omit( (-1)*log(cosmic.high.deg$pvalue, 10)), na.omit((-1)*log(cosmic.low.deg$pvalue, 10)))*1.1, parse=FALSE, label = paste0("Signature.low cell lines: ", nrow(cosmic.share.signature.low)), color = "red", hjust = 0) + 
    geom_label_repel(aes(label=ifelse(cosmic.low.deg$pvalue < cutoff.pvalue & abs(cosmic.low.deg$diff.low2other) > 100*cutoff.diff, as.character(rownames(cosmic.low.deg)), "")), size = 2, color = ifelse(cosmic.low.deg$diff.low2other > 0, "red", "blue"), segment.size=0.2) +
    labs(x="COSMIC difference", y="-log10(p value)", title=paste0("Signature [", signature.name, "] Low vs. others"))
    p2 <- p2 + theme_classic() + rremove("legend")
    # p2

    # Arranging the plot using cowplot
    p = suppressWarnings(plot_grid(p1, p2, ncol = 2, align = "hv", rel_widths = c(1,1), rel_heights = c(1,1)))
    # p <- ggplotGrob(p)
    return(p)
}

