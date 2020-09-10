##' 'dependency' compares the genetic dependency between cancer cell lines with highest and lowest dependencies of a gene set (signature).
##'
##'
##' @title dependency
##' @param signature.name Names of a signature (format: character)
##' @param signature Gene names of a signature (format: vector)
##' @param cutoff.freq Cutoff for frequency of cancer cell lines for each cancer type, default 10
##' @param cutoff.percentile Cutoff for percentile of cancer cell lines with highest/lowest dependency, default 0.2
##' @param cutoff.qvalue Cutoff for q-value of the T-test results, default 0.1
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
##' dependency(signature.name, signature)


## Main
dependency <- function(signature.name, 
                       signature, 
                       cutoff.freq        = 10, 
                       cutoff.percentile  = 0.2, 
                       cutoff.qvalue      = 0.1, 
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

    # Dependency
    # dependency.share = dep.t[,!colnames(dep.t) %in% signature,drop=FALSE]
    dependency.share = dep.t
    head(dependency.share)
    dim(dependency.share)
    # 558 17536

    dependency.share.signature.high = dependency.share[rownames(dependency.share) %in% rownames(dep.t.signature.high),,drop=FALSE]
    dependency.share.signature.high$signature = paste0(signature.name, ".dep.high")
    head(dependency.share.signature.high)
    dim(dependency.share.signature.high)
    # 27 267
    dependency.share.signature.low = dependency.share[rownames(dependency.share) %in% rownames(dep.t.signature.low),,drop=FALSE]
    dependency.share.signature.low$signature = paste0(signature.name, ".dep.low")
    head(dependency.share.signature.low)
    dim(dependency.share.signature.low)
    # 33 267
    dependency.share.signature.mid = dependency.share[!rownames(dependency.share) %in% c(rownames(dep.t.signature.high), rownames(dep.t.signature.low)),,drop=FALSE]
    dependency.share.signature.mid$signature = paste0(signature.name, ".dep.mid")
    head(dependency.share.signature.mid)
    dim(dependency.share.signature.mid)
    # 267 267

    dependency.share.signature = rbind(dependency.share.signature.high, dependency.share.signature.mid, dependency.share.signature.low)
    head(dependency.share.signature)
    dim(dependency.share.signature)
    # 327 267
    # write.csv(dependency.share.signature, paste0("dependency/dep_", signature.name, "_hml", cutoff.percentile, "_dependency.csv"), quote=FALSE)

    dependency.high = dependency.share.signature[dependency.share.signature$signature %like% "high",,drop=FALSE]
    head(dependency.high)
    dim(dependency.high)
    # 27 267
    dependency.mid = dependency.share.signature[dependency.share.signature$signature %like% "mid",,drop=FALSE]
    dependency.low = dependency.share.signature[dependency.share.signature$signature %like% "low",,drop=FALSE]

    # signature.high
    dependency.high.p = apply(dependency.share.signature[1:(ncol(dependency.share.signature)-1)], 2, function(x) signif(t.test(x[1:nrow(dependency.high)], x[(nrow(dependency.high)+1):nrow(dependency.share.signature)])$p.value,5))
    head(dependency.high.p)
    length(dependency.high.p)
    # 266
    dependency.high.q = signif(p.adjust(dependency.high.p, "BH"),5)
    head(dependency.high.q)
    length(dependency.high.q)
    # 266
    dependency.high.mean1 = apply(dependency.share.signature[1:(ncol(dependency.share.signature)-1)],  2, function(x) signif(mean(na.omit(x[1:nrow(dependency.high)])), 5))
    dependency.high.mean2 = apply(dependency.share.signature[1:(ncol(dependency.share.signature)-1)],  2, function(x) signif(mean(na.omit(x[(nrow(dependency.high)+1):nrow(dependency.share.signature)])), 5))
    dependency.high.diff = apply(dependency.share.signature[1:(ncol(dependency.share.signature)-1)], 2, function(x) signif(mean(na.omit(x[1:nrow(dependency.high)])) - mean(na.omit(x[(nrow(dependency.high)+1):nrow(dependency.share.signature)])), 5))
    head(dependency.high.diff)
    length(dependency.high.diff)
    # 266

    dependency.high.deg = data.frame(dependency.high.mean1, dependency.high.mean2, dependency.high.diff, dependency.high.p, dependency.high.q, row.names=names(dependency.high.diff))
    colnames(dependency.high.deg) = c("dependency.high", "dependency.other", "diff.high2other", "pvalue", "qvalue")
    head(dependency.high.deg)
    dim(dependency.high.deg)
    # 4338    5
    # write.csv(dependency.high.deg, paste0("dependency/dep_", signature.name, "_hml", cutoff.percentile, "_dependency.high.deg.csv"), quote=TRUE)

    # signature.low
    dependency.low.p = apply(dependency.share.signature[1:(ncol(dependency.share.signature)-1)], 2, function(x) signif(t.test(x[1:(nrow(dependency.share.signature)-nrow(dependency.low))], x[(nrow(dependency.share.signature)-nrow(dependency.low)+1):nrow(dependency.share.signature)])$p.value,5))
    head(dependency.low.p)
    length(dependency.low.p)
    # 266
    dependency.low.q = signif(p.adjust(dependency.low.p, "BH"),5)
    head(dependency.low.q)
    length(dependency.low.q)
    # 266
    dependency.low.mean1 = apply(dependency.share.signature[1:(ncol(dependency.share.signature)-1)],  2, function(x) signif(mean(na.omit(x[(nrow(dependency.share.signature)-nrow(dependency.low)+1):nrow(dependency.share.signature)])), 5))
    dependency.low.mean2 = apply(dependency.share.signature[1:(ncol(dependency.share.signature)-1)],  2, function(x) signif(mean(na.omit(x[1:(nrow(dependency.share.signature)-nrow(dependency.low))])), 5))
    dependency.low.diff = apply(dependency.share.signature[1:(ncol(dependency.share.signature)-1)], 2, function(x) signif(mean(na.omit(x[(nrow(dependency.share.signature)-nrow(dependency.low)+1):nrow(dependency.share.signature)])) - mean(na.omit(x[1:(nrow(dependency.share.signature)-nrow(dependency.low))])), 5))
    head(dependency.low.diff)
    length(dependency.low.diff)
    # 266

    dependency.low.deg = data.frame(dependency.low.mean1, dependency.low.mean2, dependency.low.diff, dependency.low.p, dependency.low.q, row.names=names(dependency.low.diff))
    colnames(dependency.low.deg) = c("dependency.low", "dependency.other", "diff.low2other", "pvalue", "qvalue")
    head(dependency.low.deg)
    dim(dependency.low.deg)
    # 4338    5
    # write.csv(dependency.low.deg, paste0("dependency/dep_", signature.name, "_hml", cutoff.percentile, "_dependency.low.deg.csv"), quote=TRUE)

    # cutoff.qvalue = 0.1
    # cutoff.diff = 0.1
    set.seed(42)
    p1 <- ggplot(data = dependency.high.deg, mapping = aes(x = diff.high2other, y = (-1)*log(qvalue, 10))) + 
    geom_point(size=1, color= ifelse(dependency.high.deg$qvalue < cutoff.qvalue & dependency.high.deg$diff.high2other > cutoff.diff, "blue", ifelse(dependency.high.deg$qvalue < cutoff.qvalue & dependency.high.deg$diff.high2other < cutoff.diff*(-1), "red", "grey60")))+ 
    # xlim(-0.5,0.5) + 
    # ylim(0,4) + 
    geom_hline(yintercept = (-1)*log(cutoff.qvalue, 10), linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = (-1)*cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    annotate("text", x=min(na.omit(dependency.high.deg$diff.high2other), na.omit(dependency.high.deg$diff.low2other)), y=max(na.omit((-1)*log(dependency.low.deg$qvalue, 10)), na.omit((-1)*log(dependency.high.deg$qvalue, 10)))*1.1, parse=FALSE, label = paste0("Signature.high cell lines: ", nrow(dependency.share.signature.high)), color = "red", hjust = 0) + 
    geom_label_repel(aes(label=ifelse(dependency.high.deg$qvalue < cutoff.qvalue & abs(dependency.high.deg$diff.high2other) > cutoff.diff, as.character(rownames(dependency.high.deg)), "")), size = 2, color = ifelse(dependency.high.deg$diff.high2other > 0, "blue", "red"), segment.size=0.2) +
    labs(x="Dependency difference", y="-log10(q value)", title=paste0("Signature [", signature.name, "] High vs. others"))
    p1 <- p1 + rremove("legend") + theme_classic()
    # p1

    p2 <- ggplot(data = dependency.low.deg, mapping = aes(x = diff.low2other, y = (-1)*log(qvalue, 10))) + 
    geom_point(size=1, color= ifelse(dependency.low.deg$qvalue < cutoff.qvalue & dependency.low.deg$diff.low2other > cutoff.diff, "blue", ifelse(dependency.low.deg$qvalue < cutoff.qvalue & dependency.low.deg$diff.low2other < cutoff.diff*(-1), "red", "grey60")))+ 
    # xlim(-0.5,0.5) + 
    # ylim(0,4) + 
    geom_hline(yintercept = (-1)*log(cutoff.qvalue, 10), linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = (-1)*cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    annotate("text", x=min(na.omit(dependency.high.deg$diff.high2other), na.omit(dependency.high.deg$diff.low2other)), y=max(na.omit((-1)*log(dependency.low.deg$qvalue, 10)), na.omit((-1)*log(dependency.high.deg$qvalue, 10)))*1.1, parse=FALSE, label = paste0("Signature.low cell lines: ", nrow(dependency.share.signature.low)), color = "red", hjust = 0) + 
    geom_label_repel(aes(label=ifelse(dependency.low.deg$qvalue < cutoff.qvalue & abs(dependency.low.deg$diff.low2other) > cutoff.diff, as.character(rownames(dependency.low.deg)), "")), size = 2, color = ifelse(dependency.low.deg$diff.low2other > 0, "blue", "red"), segment.size=0.2) +
    labs(x="Dependency difference", y="-log10(q value)", title=paste0("Signature [", signature.name, "] Low vs. others"))
    p2 <- p2 + rremove("legend") + theme_classic()
    # p2

    # Arranging the plot using cowplot
    p = suppressWarnings(plot_grid(p1, p2, ncol = 2, align = "hv", rel_widths = c(1,1), rel_heights = c(1,1)))
    # p <- ggplotGrob(p)
    return(p)
}

