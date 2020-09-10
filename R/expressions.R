##' 'expressions' compares the gene expression between cancer cell lines with highest and lowest dependencies of a gene set (signature).
##'
##'
##' @title expressions
##' @param signature.name Names of a signature (format: character)
##' @param signature Gene names of a signature (format: vector)
##' @param cutoff.freq Cutoff for frequency of cancer cell lines for each cancer type, default 10
##' @param cutoff.percentile Cutoff for percentile of cancer cell lines with highest/lowest dependency, default 0.2
##' @param cutoff.qvalue Cutoff for q-value of the T-test results, default 0.1
##' @param cutoff.fc Cutoff for fold changes of the T-test results, default 2
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
##' expressions(signature.name, signature)


## Main
expressions <- function(signature.name, 
                        signature, 
                        cutoff.freq        = 10, 
                        cutoff.percentile  = 0.2, 
                        cutoff.qvalue      = 0.1, 
                        cutoff.fc          = 2
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

    # Expression
    expression.share = exp.TPM.t
    head(expression.share)
    dim(expression.share)
    # 554 16950

    expression.share.signature.high = expression.share[rownames(expression.share) %in% rownames(dep.t.signature.high),,drop=FALSE]
    expression.share.signature.high$signature = paste0(signature.name, ".dep.high")
    head(expression.share.signature.high)
    dim(expression.share.signature.high)
    # 27 267
    expression.share.signature.low = expression.share[rownames(expression.share) %in% rownames(dep.t.signature.low),,drop=FALSE]
    expression.share.signature.low$signature = paste0(signature.name, ".dep.low")
    head(expression.share.signature.low)
    dim(expression.share.signature.low)
    # 33 267
    expression.share.signature.mid = expression.share[!rownames(expression.share) %in% c(rownames(dep.t.signature.high), rownames(dep.t.signature.low)),,drop=FALSE]
    expression.share.signature.mid$signature = paste0(signature.name, ".dep.mid")
    head(expression.share.signature.mid)
    dim(expression.share.signature.mid)
    # 267 267

    expression.share.signature = rbind(expression.share.signature.high, expression.share.signature.mid, expression.share.signature.low)
    head(expression.share.signature)
    dim(expression.share.signature)
    # 327 267
    # write.csv(expression.share.signature, paste0("expression/dep_", signature.name, "_hml", cutoff.percentile, "_expression.csv"), quote=FALSE)

    expression.high = expression.share.signature[expression.share.signature$signature %like% "high",,drop=FALSE]
    head(expression.high)
    dim(expression.high)
    # 27 267
    expression.mid = expression.share.signature[expression.share.signature$signature %like% "mid",,drop=FALSE]
    expression.low = expression.share.signature[expression.share.signature$signature %like% "low",,drop=FALSE]

    # signature.high
    expression.high.p = apply(expression.share.signature[1:(ncol(expression.share.signature)-1)], 2, function(x) signif(t.test(x[1:nrow(expression.high)], x[(nrow(expression.high)+1):nrow(expression.share.signature)])$p.value,5))
    head(expression.high.p)
    length(expression.high.p)
    # 266
    expression.high.q = signif(p.adjust(expression.high.p, "BH"),5)
    head(expression.high.q)
    length(expression.high.q)
    # 266
    expression.high.mean1 = apply(expression.share.signature[1:(ncol(expression.share.signature)-1)],  2, function(x) signif(mean(na.omit(x[1:nrow(expression.high)])), 5))
    expression.high.mean2 = apply(expression.share.signature[1:(ncol(expression.share.signature)-1)],  2, function(x) signif(mean(na.omit(x[(nrow(expression.high)+1):nrow(expression.share.signature)])), 5))
    expression.high.fc = apply(expression.share.signature[1:(ncol(expression.share.signature)-1)], 2, function(x) signif(mean(na.omit(x[1:nrow(expression.high)])) / mean(na.omit(x[(nrow(expression.high)+1):nrow(expression.share.signature)])), 5))
    head(expression.high.fc)
    length(expression.high.fc)
    # 266

    expression.high.deg = data.frame(expression.high.mean1, expression.high.mean2, expression.high.fc, expression.high.p, expression.high.q, row.names=names(expression.high.fc))
    colnames(expression.high.deg) = c("expression.high", "expression.other", "fc.high2other", "pvalue", "qvalue")
    head(expression.high.deg)
    dim(expression.high.deg)
    # 4338    5
    # write.csv(expression.high.deg, paste0("expression/dep_", signature.name, "_hml", cutoff.percentile, "_expression.high.deg.csv"), quote=TRUE)

    # signature.low
    expression.low.p = apply(expression.share.signature[1:(ncol(expression.share.signature)-1)], 2, function(x) signif(t.test(x[1:(nrow(expression.share.signature)-nrow(expression.low))], x[(nrow(expression.share.signature)-nrow(expression.low)+1):nrow(expression.share.signature)])$p.value,5))
    head(expression.low.p)
    length(expression.low.p)
    # 266
    expression.low.q = signif(p.adjust(expression.low.p, "BH"),5)
    head(expression.low.q)
    length(expression.low.q)
    # 266
    expression.low.mean1 = apply(expression.share.signature[1:(ncol(expression.share.signature)-1)],  2, function(x) signif(mean(na.omit(x[(nrow(expression.share.signature)-nrow(expression.low)+1):nrow(expression.share.signature)])), 5))
    expression.low.mean2 = apply(expression.share.signature[1:(ncol(expression.share.signature)-1)],  2, function(x) signif(mean(na.omit(x[1:(nrow(expression.share.signature)-nrow(expression.low))])), 5))
    expression.low.fc = apply(expression.share.signature[1:(ncol(expression.share.signature)-1)], 2, function(x) signif(mean(na.omit(x[(nrow(expression.share.signature)-nrow(expression.low)+1):nrow(expression.share.signature)])) / mean(na.omit(x[1:(nrow(expression.share.signature)-nrow(expression.low))])), 5))
    head(expression.low.fc)
    length(expression.low.fc)
    # 266

    expression.low.deg = data.frame(expression.low.mean1, expression.low.mean2, expression.low.fc, expression.low.p, expression.low.q, row.names=names(expression.low.fc))
    colnames(expression.low.deg) = c("expression.low", "expression.other", "fc.low2other", "pvalue", "qvalue")
    head(expression.low.deg)
    dim(expression.low.deg)
    # 4338    5
    # write.csv(expression.low.deg, paste0("expression/dep_", signature.name, "_hml", cutoff.percentile, "_expression.low.deg.csv"), quote=TRUE)

    # expression.high.deg = read.csv(paste0("expression/dep_", signature.name, "_hml", cutoff.percentile, "_expression.high.deg.csv"), header=TRUE, row.names=1)
    # expression.low.deg = read.csv(paste0("expression/dep_", signature.name, "_hml", cutoff.percentile, "_expression.low.deg.csv"), header=TRUE, row.names=1)
    # cutoff.qvalue = 0.1
    # cutoff.fc = 2
    set.seed(42)
    p1 <- ggplot(data = expression.high.deg, mapping = aes(x = log(fc.high2other, 2), y = (-1)*log(qvalue, 10))) + 
    geom_point(size=0.5, color= ifelse(expression.high.deg$qvalue < cutoff.qvalue & expression.high.deg$fc.high2other > cutoff.fc, "red", ifelse(expression.high.deg$qvalue < cutoff.qvalue & expression.high.deg$fc.high2other < 1/cutoff.fc, "blue", "grey60")))+ 
    # xlim(-1,1) + 
    # ylim(0,4) + 
    geom_hline(yintercept = (-1)*log(cutoff.qvalue, 10), linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = log(cutoff.fc,2), linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = (-1)*log(cutoff.fc,2), linetype="dashed", colour="grey30", size=0.2) + 
    annotate("text", x=min(log(expression.high.deg$fc.high2other[is.finite(expression.high.deg$fc.high2other) & expression.high.deg$fc.high2other > 0], 2))*0.9, y=max(na.omit((-1)*log(expression.low.deg$qvalue, 10)), na.omit((-1)*log(expression.high.deg$qvalue, 10)))*1.1, parse=FALSE, label = paste0("Signature.high cell lines: ", nrow(expression.share.signature.high)), color = "red", hjust = 0) + 
    # geom_label_repel(aes(label=ifelse(expression.high.deg$qvalue < (cutoff.qvalue) & abs(log(expression.high.deg$fc.high2other,2)) > log(cutoff.fc,2), as.character(rownames(expression.high.deg)), "")), size = 2, color = ifelse(expression.high.deg$fc.high2other < 0, "blue", "red"), segment.size=0.2) +
    labs(x="Expression fold change (log2)", y="-log10(q value)", title=paste0("Signature [", signature.name, "] High vs. others"))
    p1 <- p1 + rremove("legend") + theme_classic()
    # p1

    p2 <- ggplot(data = expression.low.deg, mapping = aes(x = log(fc.low2other, 2), y = (-1)*log(qvalue, 10))) + 
    geom_point(size=0.5, color= ifelse(expression.low.deg$qvalue < cutoff.qvalue & expression.low.deg$fc.low2other > cutoff.fc, "red", ifelse(expression.low.deg$qvalue < cutoff.qvalue & expression.low.deg$fc.low2other < 1/cutoff.fc, "blue", "grey60")))+ 
    # xlim(-1,1) + 
    # ylim(0,4) + 
    geom_hline(yintercept = (-1)*log(cutoff.qvalue, 10), linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = log(cutoff.fc,2), linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = (-1)*log(cutoff.fc,2), linetype="dashed", colour="grey30", size=0.2) + 
    annotate("text", x=min(log(expression.low.deg$fc.low2other[is.finite(expression.low.deg$fc.low2other) & expression.low.deg$fc.low2other > 0], 2))*0.9, y=max(na.omit((-1)*log(expression.low.deg$qvalue, 10)), na.omit((-1)*log(expression.high.deg$qvalue, 10)))*1.1, parse=FALSE, label = paste0("Signature.low cell lines: ", nrow(expression.share.signature.low)), color = "red", hjust = 0) + 
    # geom_label_repel(aes(label=ifelse(expression.low.deg$qvalue < (cutoff.qvalue) & abs(log(expression.low.deg$fc.low2other,2)) > log(cutoff.fc,2), as.character(rownames(expression.low.deg)), "")), size = 2, color = ifelse(expression.low.deg$fc.low2other < 0, "blue", "red"), segment.size=0.2) +
    labs(x="Expression fold change (log2)", y="-log10(q value)", title=paste0("Signature [", signature.name, "] Low vs. others"))
    p2 <- p2 + rremove("legend") + theme_classic()
    # p2

    # Arranging the plot using cowplot
    p = suppressWarnings(plot_grid(p1, p2, ncol = 2, align = "hv", rel_widths = c(1,1), rel_heights = c(1,1)))
    # p <- ggplotGrob(p)
    return(p)
}

