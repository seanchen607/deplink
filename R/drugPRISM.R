##' 'drugPRISM' compares the drug sensitivity (PRISM dataset) between cancer cell lines with highest and lowest dependencies of a gene set (signature).
##'
##'
##' @title drugPRISM
##' @param signature.name Names of a signature (format: character)
##' @param signature Gene names of a signature (format: vector)
##' @param cutoff.freq Cutoff for frequency of cancer cell lines for each cancer type, default 10
##' @param cutoff.percentile Cutoff for percentile of cancer cell lines with highest/lowest dependency, default 0.2
##' @param cutoff.pvalue Cutoff for p-value of the T-test results, default 0.01
##' @param cutoff.diff Cutoff for difference of the T-test results, default 0.7
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
##' drugPRISM(signature.name, signature)


## Main
drugPRISM <- function(signature.name, 
                      signature, 
                      cutoff.freq        = 10, 
                      cutoff.percentile  = 0.2, 
                      cutoff.pvalue      = 0.01, 
                      cutoff.diff        = 0.7
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

    # Drug sensitivity - PRISM
    drug2.share = drug2[rownames(drug2) %in% rownames(dep.t),,drop=FALSE]
    head(drug2.share)
    dim(drug2.share)
    # 359 4686

    drug2.share.signature.high = drug2.share[rownames(drug2.share) %in% rownames(dep.t.signature.high),,drop=FALSE]
    drug2.share.signature.high$signature = paste0(signature.name, ".dep.high")
    head(drug2.share.signature.high)
    dim(drug2.share.signature.high)
    # 23 4687
    drug2.share.signature.low = drug2.share[rownames(drug2.share) %in% rownames(dep.t.signature.low),,drop=FALSE]
    drug2.share.signature.low$signature = paste0(signature.name, ".dep.low")
    head(drug2.share.signature.low)
    dim(drug2.share.signature.low)
    # 39 4687
    drug2.share.signature.mid = drug2.share[!rownames(drug2.share) %in% c(rownames(dep.t.signature.high), rownames(dep.t.signature.low)),,drop=FALSE]
    drug2.share.signature.mid$signature = paste0(signature.name, ".dep.mid")
    head(drug2.share.signature.mid)
    dim(drug2.share.signature.mid)
    # 297 4687

    drug2.share.signature = rbind(drug2.share.signature.high, drug2.share.signature.mid, drug2.share.signature.low)
    head(drug2.share.signature)
    dim(drug2.share.signature)
    # 359 4687
    # write.csv(drug2.share.signature, paste0("drug.PRISM/dep_", signature.name, "_hml", cutoff.percentile, "_drug.csv"), quote=FALSE)

    drug2.high = drug2.share.signature[drug2.share.signature$signature %like% "high",,drop=FALSE]
    head(drug2.high)
    dim(drug2.high)
    # 27 267
    drug2.mid = drug2.share.signature[drug2.share.signature$signature %like% "mid",,drop=FALSE]
    drug2.low = drug2.share.signature[drug2.share.signature$signature %like% "low",,drop=FALSE]

    # signature.high
    drug2.high.p = apply(drug2.share.signature[1:(ncol(drug2.share.signature)-1)], 2, function(x) ifelse(length(x[1:nrow(drug2.high)][!is.na(x[1:nrow(drug2.high)])])>1, signif(t.test(x[1:nrow(drug2.high)], x[(nrow(drug2.high)+1):nrow(drug2.share.signature)])$p.value,5), NA))
    head(drug2.high.p)
    length(drug2.high.p)
    # 266
    drug2.high.q = signif(p.adjust(drug2.high.p, "BH"),5)
    head(drug2.high.q)
    length(drug2.high.q)
    # 266
    drug2.high.mean1 = apply(drug2.share.signature[1:(ncol(drug2.share.signature)-1)],  2, function(x) signif(mean(na.omit(x[1:nrow(drug2.high)])), 5))
    drug2.high.mean2 = apply(drug2.share.signature[1:(ncol(drug2.share.signature)-1)],  2, function(x) signif(mean(na.omit(x[(nrow(drug2.high)+1):nrow(drug2.share.signature)])), 5))
    drug2.high.diff = apply(drug2.share.signature[1:(ncol(drug2.share.signature)-1)], 2, function(x) signif(mean(na.omit(x[1:nrow(drug2.high)])) - mean(na.omit(x[(nrow(drug2.high)+1):nrow(drug2.share.signature)])), 5))
    head(drug2.high.diff)
    length(drug2.high.diff)
    # 266

    drug2.high.deg = data.frame(drug2.high.mean1, drug2.high.mean2, drug2.high.diff, drug2.high.p, drug2.high.q, row.names=names(drug2.high.diff))
    colnames(drug2.high.deg) = c("drug2.high", "drug2.other", "diff.high2other", "pvalue", "qvalue")
    head(drug2.high.deg)
    dim(drug2.high.deg)
    # 266   5

    drug2.high.deg = merge(drug2.high.deg, drug2.meta, by="row.names", all.x=TRUE)
    rownames(drug2.high.deg) = drug2.high.deg[,1]
    drug2.high.deg = drug2.high.deg[,-1]
    head(drug2.high.deg)
    dim(drug2.high.deg)
    # 266   9
    # write.csv(drug2.high.deg, paste0("drug.PRISM/dep_", signature.name, "_hml", cutoff.percentile, "_drug.high.deg.csv"), quote=TRUE)

    # signature.low
    drug2.low.p = apply(drug2.share.signature[1:(ncol(drug2.share.signature)-1)], 2, function(x) ifelse(length(x[(nrow(drug2.share.signature)-nrow(drug2.low)+1):nrow(drug2.share.signature)][!is.na(x[(nrow(drug2.share.signature)-nrow(drug2.low)+1):nrow(drug2.share.signature)])])>1, signif(t.test(x[1:(nrow(drug2.share.signature)-nrow(drug2.low))], x[(nrow(drug2.share.signature)-nrow(drug2.low)+1):nrow(drug2.share.signature)])$p.value,5), NA))
    head(drug2.low.p)
    length(drug2.low.p)
    # 266
    drug2.low.q = signif(p.adjust(drug2.low.p, "BH"),5)
    head(drug2.low.q)
    length(drug2.low.q)
    # 266
    drug2.low.mean1 = apply(drug2.share.signature[1:(ncol(drug2.share.signature)-1)],  2, function(x) signif(mean(na.omit(x[(nrow(drug2.share.signature)-nrow(drug2.low)):(nrow(drug2.share.signature)-1)])), 5))
    drug2.low.mean2 = apply(drug2.share.signature[1:(ncol(drug2.share.signature)-1)],  2, function(x) signif(mean(na.omit(x[1:(nrow(drug2.share.signature)-nrow(drug2.low)-1)])), 5))
    drug2.low.diff = apply(drug2.share.signature[1:(ncol(drug2.share.signature)-1)], 2, function(x) signif(mean(na.omit(x[(nrow(drug2.share.signature)-nrow(drug2.low)+1):nrow(drug2.share.signature)])) - mean(na.omit(x[1:(nrow(drug2.share.signature)-nrow(drug2.low))])), 5))
    head(drug2.low.diff)
    length(drug2.low.diff)
    # 266

    drug2.low.deg = data.frame(drug2.low.mean1, drug2.low.mean2, drug2.low.diff, drug2.low.p, drug2.low.q, row.names=names(drug2.low.diff))
    colnames(drug2.low.deg) = c("drug2.low", "drug2.other", "diff.low2other", "pvalue", "qvalue")
    head(drug2.low.deg)
    dim(drug2.low.deg)
    # 266   5

    drug2.low.deg = merge(drug2.low.deg, drug2.meta, by="row.names", all.x=TRUE)
    rownames(drug2.low.deg) = drug2.low.deg[,1]
    drug2.low.deg = drug2.low.deg[,-1]
    head(drug2.low.deg)
    dim(drug2.low.deg)
    # 266   9
    # write.csv(drug2.low.deg, paste0("drug.PRISM/dep_", signature.name, "_hml", cutoff.percentile, "_drug.low.deg.csv"), quote=TRUE)

    drug2.high.deg$moa = gsub(",\\(.+\\)", "", drug2.high.deg$moa)
    drug2.low.deg$moa = gsub(",\\(.+\\)", "", drug2.low.deg$moa)

    # cutoff.pvalue = 0.01
    # cutoff.diff = 0.2
    set.seed(42)
    p1 <- ggplot(data = drug2.high.deg, mapping = aes(x = diff.high2other, y = (-1)*log(pvalue, 10))) + 
    geom_point(size=1, color= ifelse(drug2.high.deg$pvalue < cutoff.pvalue & drug2.high.deg$diff.high2other > cutoff.diff, "blue", ifelse(drug2.high.deg$pvalue < cutoff.pvalue & drug2.high.deg$diff.high2other < cutoff.diff*(-1), "red", "grey60")))+ 
    # xlim(-1.5,1.5) + 
    # ylim(0,4) + 
    geom_hline(yintercept = (-1)*log(cutoff.pvalue, 10), linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = (-1)*cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    annotate("text", x=min(na.omit(drug2.high.deg$diff.high2other), na.omit(drug2.high.deg$diff.low2other)), y=max(na.omit((-1)*log(drug2.high.deg$pvalue, 10)), na.omit((-1)*log(drug2.low.deg$pvalue, 10)))*1.1, parse=FALSE, label = paste0("Signature.high cell lines: ", nrow(drug2.share.signature.high)), color = "red", hjust = 0) + 
    geom_label_repel(aes(label=ifelse(drug2.high.deg$pvalue < cutoff.pvalue & abs(drug2.high.deg$diff.high2other) > cutoff.diff, as.character(drug2.high.deg$name), "")), size = 2, color = ifelse(drug2.high.deg$diff.high2other > 0, "blue", "red"), segment.size=0.2) +
    labs(x="Sensitivity difference", y="-log10(p value)", title=paste0("Signature [", signature.name, "] High vs. others"))
    p1 <- p1 + theme_classic() + rremove("legend")
    # p1

    p2 <- ggplot(data = drug2.low.deg, mapping = aes(x = diff.low2other, y = (-1)*log(pvalue, 10))) + 
    geom_point(size=1, color= ifelse(drug2.low.deg$pvalue < cutoff.pvalue & drug2.low.deg$diff.low2other > cutoff.diff, "blue", ifelse(drug2.low.deg$pvalue < cutoff.pvalue & drug2.low.deg$diff.low2other < cutoff.diff*(-1), "red", "grey60")))+ 
    # xlim(-1.5,1.5) + 
    # ylim(0,4) + 
    geom_hline(yintercept = (-1)*log(cutoff.pvalue, 10), linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = (-1)*cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    annotate("text", x=min(na.omit(drug2.high.deg$diff.high2other), na.omit(drug2.high.deg$diff.low2other)), y=max(na.omit((-1)*log(drug2.high.deg$pvalue, 10)), na.omit((-1)*log(drug2.low.deg$pvalue, 10)))*1.1, parse=FALSE, label = paste0("Signature.low cell lines: ", nrow(drug2.share.signature.low)), color = "red", hjust = 0) + 
    geom_label_repel(aes(label=ifelse(drug2.low.deg$pvalue < cutoff.pvalue & abs(drug2.low.deg$diff.low2other) > cutoff.diff, as.character(drug2.low.deg$name), "")), size = 2, color = ifelse(drug2.low.deg$diff.low2other > 0, "blue", "red"), segment.size=0.2) +
    labs(x="Sensitivity difference", y="-log10(p value)", title=paste0("Signature [", signature.name, "] Low vs. others"))
    p2 <- p2 + theme_classic() + rremove("legend")
    # p2

    p3 <- ggplot(data = drug2.high.deg, mapping = aes(x = diff.high2other, y = (-1)*log(pvalue, 10))) + 
    geom_point(size=1, color= ifelse(drug2.high.deg$pvalue < cutoff.pvalue & drug2.high.deg$diff.high2other > cutoff.diff, "blue", ifelse(drug2.high.deg$pvalue < cutoff.pvalue & drug2.high.deg$diff.high2other < cutoff.diff*(-1), "red", "grey60")))+ 
    # xlim(-1.5,1.5) + 
    # ylim(0,4) + 
    geom_hline(yintercept = (-1)*log(cutoff.pvalue, 10), linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = (-1)*cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    annotate("text", x=min(na.omit(drug2.high.deg$diff.high2other), na.omit(drug2.high.deg$diff.low2other)), y=max(na.omit((-1)*log(drug2.high.deg$pvalue, 10)), na.omit((-1)*log(drug2.low.deg$pvalue, 10)))*1.1, parse=FALSE, label = paste0("Signature.high cell lines: ", nrow(drug2.share.signature.high)), color = "red", hjust = 0) + 
    geom_label_repel(aes(label=ifelse(drug2.high.deg$pvalue < cutoff.pvalue & abs(drug2.high.deg$diff.high2other) > cutoff.diff, as.character(drug2.high.deg$moa), "")), size = 2, color = ifelse(drug2.high.deg$diff.high2other > 0, "blue", "red"), segment.size=0.2) +
    labs(x="Sensitivity difference", y="-log10(p value)", title=paste0("Signature [", signature.name, "] High vs. others"))
    p3 <- p3 + theme_classic() + rremove("legend")
    # p3

    p4 <- ggplot(data = drug2.low.deg, mapping = aes(x = diff.low2other, y = (-1)*log(pvalue, 10))) + 
    geom_point(size=1, color= ifelse(drug2.low.deg$pvalue < cutoff.pvalue & drug2.low.deg$diff.low2other > cutoff.diff, "blue", ifelse(drug2.low.deg$pvalue < cutoff.pvalue & drug2.low.deg$diff.low2other < cutoff.diff*(-1), "red", "grey60")))+ 
    # xlim(-1.5,1.5) + 
    # ylim(0,4) + 
    geom_hline(yintercept = (-1)*log(cutoff.pvalue, 10), linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = (-1)*cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    annotate("text", x=min(na.omit(drug2.high.deg$diff.high2other), na.omit(drug2.high.deg$diff.low2other)), y=max(na.omit((-1)*log(drug2.high.deg$pvalue, 10)), na.omit((-1)*log(drug2.low.deg$pvalue, 10)))*1.1, parse=FALSE, label = paste0("Signature.low cell lines: ", nrow(drug2.share.signature.low)), color = "red", hjust = 0) + 
    geom_label_repel(aes(label=ifelse(drug2.low.deg$pvalue < cutoff.pvalue & abs(drug2.low.deg$diff.low2other) > cutoff.diff, as.character(drug2.low.deg$moa), "")), size = 2, color = ifelse(drug2.low.deg$diff.low2other > 0, "blue", "red"), segment.size=0.2) +
    labs(x="Sensitivity difference", y="-log10(p value)", title=paste0("Signature [", signature.name, "] Low vs. others"))
    p4 <- p4 + theme_classic() + rremove("legend")
    # p4

    # Arranging the plot using cowplot
    p = suppressWarnings(plot_grid(p1, p2, p3, p4, ncol = 2, align = "hv", rel_widths = c(1,1), rel_heights = c(1,1)))
    # p <- ggplotGrob(p)
    return(p)
}

