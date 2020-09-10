##' 'drugGDSC' compares the drug sensitivity (GDSC dataset) between cancer cell lines with highest and lowest dependencies of a gene set (signature).
##'
##'
##' @title drugGDSC
##' @param signature.name Names of a signature (format: character)
##' @param signature Gene names of a signature (format: vector)
##' @param cutoff.freq Cutoff for frequency of cancer cell lines for each cancer type, default 10
##' @param cutoff.percentile Cutoff for percentile of cancer cell lines with highest/lowest dependency, default 0.2
##' @param cutoff.pvalue Cutoff for p-value of the T-test results, default 0.05
##' @param cutoff.diff Cutoff for difference of the T-test results, default 0.5
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
##' drugGDSC(signature.name, signature)


## Main
drugGDSC <- function(signature.name, 
                     signature, 
                     cutoff.freq        = 10, 
                     cutoff.percentile  = 0.2, 
                     cutoff.pvalue      = 0.05, 
                     cutoff.diff        = 0.5
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

    # Drug sensitivity - GDSC
    drug.share = drug[rownames(drug) %in% rownames(dep.t),,drop=FALSE]
    head(drug.share)
    dim(drug.share)
    # 327 266

    drug.share.signature.high = drug.share[rownames(drug.share) %in% rownames(dep.t.signature.high),,drop=FALSE]
    drug.share.signature.high$signature = paste0(signature.name, ".dep.high")
    head(drug.share.signature.high)
    dim(drug.share.signature.high)
    # 27 267
    drug.share.signature.low = drug.share[rownames(drug.share) %in% rownames(dep.t.signature.low),,drop=FALSE]
    drug.share.signature.low$signature = paste0(signature.name, ".dep.low")
    head(drug.share.signature.low)
    dim(drug.share.signature.low)
    # 33 267
    drug.share.signature.mid = drug.share[!rownames(drug.share) %in% c(rownames(dep.t.signature.high), rownames(dep.t.signature.low)),,drop=FALSE]
    drug.share.signature.mid$signature = paste0(signature.name, ".dep.mid")
    head(drug.share.signature.mid)
    dim(drug.share.signature.mid)
    # 267 267

    drug.share.signature = rbind(drug.share.signature.high, drug.share.signature.mid, drug.share.signature.low)
    head(drug.share.signature)
    dim(drug.share.signature)
    # 327 267
    # write.csv(drug.share.signature, paste0("drug.GDSC/dep_", signature.name, "_hml", cutoff.percentile, "_drug.csv"), quote=FALSE)

    drug.high = drug.share.signature[drug.share.signature$signature %like% "high",,drop=FALSE]
    head(drug.high)
    dim(drug.high)
    # 27 267
    drug.mid = drug.share.signature[drug.share.signature$signature %like% "mid",,drop=FALSE]
    drug.low = drug.share.signature[drug.share.signature$signature %like% "low",,drop=FALSE]

    # signature.high
    drug.high.p = apply(drug.share.signature[1:(ncol(drug.share.signature)-1)], 2, function(x) ifelse(length(x[1:nrow(drug.high)][!is.na(x[1:nrow(drug.high)])])>1, signif(t.test(x[1:nrow(drug.high)], x[(nrow(drug.high)+1):nrow(drug.share.signature)])$p.value,5), NA))
    head(drug.high.p)
    length(drug.high.p)
    # 266
    drug.high.q = signif(p.adjust(drug.high.p, "BH"),5)
    head(drug.high.q)
    length(drug.high.q)
    # 266
    drug.high.mean1 = apply(drug.share.signature[1:(ncol(drug.share.signature)-1)],  2, function(x) signif(mean(na.omit(x[1:nrow(drug.high)])), 5))
    drug.high.mean2 = apply(drug.share.signature[1:(ncol(drug.share.signature)-1)],  2, function(x) signif(mean(na.omit(x[(nrow(drug.high)+1):nrow(drug.share.signature)])), 5))
    drug.high.diff = apply(drug.share.signature[1:(ncol(drug.share.signature)-1)], 2, function(x) signif(mean(na.omit(x[1:nrow(drug.high)])) - mean(na.omit(x[(nrow(drug.high)+1):nrow(drug.share.signature)])), 5))
    head(drug.high.diff)
    length(drug.high.diff)
    # 266

    drug.high.deg = data.frame(drug.high.mean1, drug.high.mean2, drug.high.diff, drug.high.p, drug.high.q, row.names=names(drug.high.diff))
    colnames(drug.high.deg) = c("drug.high", "drug.other", "diff.high2other", "pvalue", "qvalue")
    head(drug.high.deg)
    dim(drug.high.deg)
    # 266   5

    drug.high.deg = merge(drug.high.deg, drug.meta, by="row.names", all.x=TRUE)
    rownames(drug.high.deg) = drug.high.deg[,1]
    drug.high.deg = drug.high.deg[,-1]
    head(drug.high.deg)
    dim(drug.high.deg)
    # 266   9
    # write.csv(drug.high.deg, paste0("drug.GDSC/dep_", signature.name, "_hml", cutoff.percentile, "_drug.high.deg.csv"), quote=TRUE)

    # signature.low
    drug.low.p = apply(drug.share.signature[1:(ncol(drug.share.signature)-1)], 2, function(x) ifelse(length(x[1:(nrow(drug.share.signature)-nrow(drug.low))][!is.na(x[1:(nrow(drug.share.signature)-nrow(drug.low))])])>1, signif(t.test(x[1:(nrow(drug.share.signature)-nrow(drug.low))], x[(nrow(drug.share.signature)-nrow(drug.low)+1):nrow(drug.share.signature)])$p.value,5), NA))
    head(drug.low.p)
    length(drug.low.p)
    # 266
    drug.low.q = signif(p.adjust(drug.low.p, "BH"),5)
    head(drug.low.q)
    length(drug.low.q)
    # 266
    drug.low.mean1 = apply(drug.share.signature[1:(ncol(drug.share.signature)-1)],  2, function(x) signif(mean(na.omit(x[(nrow(drug.share.signature)-nrow(drug.low)):(nrow(drug.share.signature)-1)])), 5))
    drug.low.mean2 = apply(drug.share.signature[1:(ncol(drug.share.signature)-1)],  2, function(x) signif(mean(na.omit(x[1:(nrow(drug.share.signature)-nrow(drug.low)-1)])), 5))
    drug.low.diff = apply(drug.share.signature[1:(ncol(drug.share.signature)-1)], 2, function(x) signif(mean(na.omit(x[(nrow(drug.share.signature)-nrow(drug.low)+1):nrow(drug.share.signature)])) - mean(na.omit(x[1:(nrow(drug.share.signature)-nrow(drug.low))])), 5))
    head(drug.low.diff)
    length(drug.low.diff)
    # 266

    drug.low.deg = data.frame(drug.low.mean1, drug.low.mean2, drug.low.diff, drug.low.p, drug.low.q, row.names=names(drug.low.diff))
    colnames(drug.low.deg) = c("drug.low", "drug.other", "diff.low2other", "pvalue", "qvalue")
    head(drug.low.deg)
    dim(drug.low.deg)
    # 266   5

    drug.low.deg = merge(drug.low.deg, drug.meta, by="row.names", all.x=TRUE)
    rownames(drug.low.deg) = drug.low.deg[,1]
    drug.low.deg = drug.low.deg[,-1]
    head(drug.low.deg)
    dim(drug.low.deg)
    # 266   9
    # write.csv(drug.low.deg, paste0("drug.GDSC/dep_", signature.name, "_hml", cutoff.percentile, "_drug.low.deg.csv"), quote=TRUE)

    drug.high.deg$DRUG_NAME = gsub(" \\(.+\\)", "", drug.high.deg$DRUG_NAME)
    drug.low.deg$DRUG_NAME = gsub(" \\(.+\\)", "", drug.low.deg$DRUG_NAME)

    # cutoff.pvalue = 0.05
    # cutoff.diff = 0.6
    set.seed(42)
    p1 <- ggplot(data = drug.high.deg, mapping = aes(x = diff.high2other, y = (-1)*log(pvalue, 10))) + 
    geom_point(size=1, color= ifelse(drug.high.deg$pvalue < cutoff.pvalue & drug.high.deg$diff.high2other > cutoff.diff, "blue", ifelse(drug.high.deg$pvalue < cutoff.pvalue & drug.high.deg$diff.high2other < cutoff.diff*(-1), "red", "grey60")))+ 
    # xlim(-1.5,1.5) + 
    # ylim(0,4) + 
    geom_hline(yintercept = (-1)*log(cutoff.pvalue, 10), linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = (-1)*cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    annotate("text", x=min(na.omit(drug.high.deg$diff.high2other), na.omit(drug.high.deg$diff.low2other)), y=max(na.omit((-1)*log(drug.high.deg$pvalue, 10)), na.omit((-1)*log(drug.low.deg$pvalue, 10)))*1.1, parse=FALSE, label = paste0("Signature.high cell lines: ", nrow(drug.share.signature.high)), color = "red", hjust = 0) + 
    geom_label_repel(aes(label=ifelse(drug.high.deg$pvalue < cutoff.pvalue & abs(drug.high.deg$diff.high2other) > cutoff.diff, as.character(drug.high.deg$DRUG_NAME), "")), size = 2, color = ifelse(drug.high.deg$diff.high2other > 0, "blue", "red"), segment.size=0.2) +
    labs(x="Sensitivity difference", y="-log10(p value)", title=paste0("Signature [", signature.name, "] High vs. others"))
    p1 <- p1 + theme_classic() + rremove("legend")
    # p1

    p2 <- ggplot(data = drug.low.deg, mapping = aes(x = diff.low2other, y = (-1)*log(pvalue, 10))) + 
    geom_point(size=1, color= ifelse(drug.low.deg$pvalue < cutoff.pvalue & drug.low.deg$diff.low2other > cutoff.diff, "blue", ifelse(drug.low.deg$pvalue < cutoff.pvalue & drug.low.deg$diff.low2other < cutoff.diff*(-1), "red", "grey60")))+ 
    # xlim(-1.5,1.5) + 
    # ylim(0,4) + 
    geom_hline(yintercept = (-1)*log(cutoff.pvalue, 10), linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = (-1)*cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    annotate("text", x=min(na.omit(drug.high.deg$diff.high2other), na.omit(drug.high.deg$diff.low2other)), y=max(na.omit((-1)*log(drug.high.deg$pvalue, 10)), na.omit((-1)*log(drug.low.deg$pvalue, 10)))*1.1, parse=FALSE, label = paste0("Signature.low cell lines: ", nrow(drug.share.signature.low)), color = "red", hjust = 0) + 
    geom_label_repel(aes(label=ifelse(drug.low.deg$pvalue < cutoff.pvalue & abs(drug.low.deg$diff.low2other) > cutoff.diff, as.character(drug.low.deg$DRUG_NAME), "")), size = 2, color = ifelse(drug.low.deg$diff.low2other > 0, "blue", "red"), segment.size=0.2) +
    labs(x="Sensitivity difference", y="-log10(p value)", title=paste0("Signature [", signature.name, "] Low vs. others"))
    p2 <- p2 + theme_classic() + rremove("legend")
    # p2

    p3 <- ggplot(data = drug.high.deg, mapping = aes(x = diff.high2other, y = (-1)*log(pvalue, 10))) + 
    geom_point(size=1, color= ifelse(drug.high.deg$pvalue < cutoff.pvalue & drug.high.deg$diff.high2other > cutoff.diff, "blue", ifelse(drug.high.deg$pvalue < cutoff.pvalue & drug.high.deg$diff.high2other < cutoff.diff*(-1), "red", "grey60")))+ 
    # xlim(-1.5,1.5) + 
    # ylim(0,4) + 
    geom_hline(yintercept = (-1)*log(cutoff.pvalue, 10), linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = (-1)*cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    annotate("text", x=min(na.omit(drug.high.deg$diff.high2other), na.omit(drug.high.deg$diff.low2other)), y=max(na.omit((-1)*log(drug.high.deg$pvalue, 10)), na.omit((-1)*log(drug.low.deg$pvalue, 10)))*1.1, parse=FALSE, label = paste0("Signature.high cell lines: ", nrow(drug.share.signature.high)), color = "red", hjust = 0) + 
    geom_label_repel(aes(label=ifelse(drug.high.deg$pvalue < cutoff.pvalue & abs(drug.high.deg$diff.high2other) > cutoff.diff, as.character(drug.high.deg$TARGET_PATHWAY), "")), size = 2, color = ifelse(drug.high.deg$diff.high2other > 0, "blue", "red"), segment.size=0.2) +
    labs(x="Sensitivity difference", y="-log10(p value)", title=paste0("Signature [", signature.name, "] High vs. others"))
    p3 <- p3 + theme_classic() + rremove("legend")
    # p3

    p4 <- ggplot(data = drug.low.deg, mapping = aes(x = diff.low2other, y = (-1)*log(pvalue, 10))) + 
    geom_point(size=1, color= ifelse(drug.low.deg$pvalue < cutoff.pvalue & drug.low.deg$diff.low2other > cutoff.diff, "blue", ifelse(drug.low.deg$pvalue < cutoff.pvalue & drug.low.deg$diff.low2other < cutoff.diff*(-1), "red", "grey60")))+ 
    # xlim(-1.5,1.5) + 
    # ylim(0,4) + 
    geom_hline(yintercept = (-1)*log(cutoff.pvalue, 10), linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = (-1)*cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    annotate("text", x=min(na.omit(drug.high.deg$diff.high2other), na.omit(drug.high.deg$diff.low2other)), y=max(na.omit((-1)*log(drug.high.deg$pvalue, 10)), na.omit((-1)*log(drug.low.deg$pvalue, 10)))*1.1, parse=FALSE, label = paste0("Signature.low cell lines: ", nrow(drug.share.signature.low)), color = "red", hjust = 0) + 
    geom_label_repel(aes(label=ifelse(drug.low.deg$pvalue < cutoff.pvalue & abs(drug.low.deg$diff.low2other) > cutoff.diff, as.character(drug.low.deg$TARGET_PATHWAY), "")), size = 2, color = ifelse(drug.low.deg$diff.low2other > 0, "blue", "red"), segment.size=0.2) +
    labs(x="Sensitivity difference", y="-log10(p value)", title=paste0("Signature [", signature.name, "] Low vs. others"))
    p4 <- p4 + theme_classic() + rremove("legend")
    # p4

    # Arranging the plot using cowplot
    p = suppressWarnings(plot_grid(p1, p2, p3, p4, ncol = 2, align = "hv", rel_widths = c(1,1), rel_heights = c(1,1)))
    # p <- ggplotGrob(p)
    return(p)
}

