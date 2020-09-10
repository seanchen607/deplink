##' 'mutations' compares the genetic mutations between cancer cell lines with highest and lowest dependencies of a gene set (signature).
##'
##'
##' @title mutations
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
##' mutations(signature.name, signature)


## Main
mutations <- function(signature.name, 
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

    # Mutation
    mutation.signature.high = mutation[rownames(mutation) %in% rownames(dep.t.signature.high),,drop=FALSE]
    mutation.signature.high.gene = mutation.signature.high[,colSums(mutation.signature.high)>1]
    mutation.signature.high$signature = paste0(signature.name, ".dep.high")
    head(mutation.signature.high)
    dim(mutation.signature.high)
    # 40 44
    mutation.signature.low = mutation[rownames(mutation) %in% rownames(dep.t.signature.low),,drop=FALSE]
    mutation.signature.low.gene = mutation.signature.low[,colSums(mutation.signature.low)>1]
    mutation.signature.low$signature = paste0(signature.name, ".dep.low")
    head(mutation.signature.low)
    dim(mutation.signature.low)
    # 46 44
    mutation.signature.mid = mutation[!rownames(mutation) %in% c(rownames(dep.t.signature.high), rownames(dep.t.signature.low)),,drop=FALSE]
    mutation.signature.mid.gene = mutation.signature.mid[,colSums(mutation.signature.mid)>1]
    mutation.signature.mid$signature = paste0(signature.name, ".dep.mid")
    head(mutation.signature.mid)
    dim(mutation.signature.mid)
    # 360  44

    gene.union = Reduce(union, list(colnames(mutation.signature.high.gene), colnames(mutation.signature.low.gene), colnames(mutation.signature.mid.gene)))
    length(gene.union)
    # 9422

    mutation.signature = rbind(mutation.signature.high, mutation.signature.mid, mutation.signature.low)
    mutation.signature = mutation.signature[,colnames(mutation.signature) %in% c(gene.union, "signature"),drop=FALSE]
    head(mutation.signature)
    dim(mutation.signature)
    # 554 9423
    # write.csv(mutation.signature, paste0("mutation/dep_", signature.name, "_hml", cutoff.percentile, "_mutation.csv"), quote=FALSE)

    mutation.high = mutation.signature[mutation.signature$signature %like% "high",,drop=FALSE]
    head(mutation.high)
    dim(mutation.high)
    # 40 44
    mutation.mid = mutation.signature[mutation.signature$signature %like% "mid",,drop=FALSE]
    mutation.low = mutation.signature[mutation.signature$signature %like% "low",,drop=FALSE]

    # signature.high
    mutation.high.p = apply(mutation.signature[1:(ncol(mutation.signature)-1)], 2, function(x) signif(t.test(na.omit(x[1:nrow(mutation.high)]), na.omit(x[(nrow(mutation.high)+1):nrow(mutation.signature)]))$p.value,5))
    head(mutation.high.p)
    length(mutation.high.p)
    # 42
    mutation.high.q = signif(p.adjust(mutation.high.p, "BH"),5)
    head(mutation.high.q)
    length(mutation.high.q)
    # 42
    mutation.high.mean1 = apply(mutation.signature[1:(ncol(mutation.signature)-1)],  2, function(x) signif(mean(na.omit(x[1:nrow(mutation.high)])), 5))
    mutation.high.mean2 = apply(mutation.signature[1:(ncol(mutation.signature)-1)],  2, function(x) signif(mean(na.omit(x[(nrow(mutation.high)+1):nrow(mutation.signature)])), 5))
    mutation.high.diff = apply(mutation.signature[1:(ncol(mutation.signature)-1)], 2, function(x) signif(mean(na.omit(x[1:nrow(mutation.high)])) - mean(na.omit(x[(nrow(mutation.high)+1):nrow(mutation.signature)])), 5))
    head(mutation.high.diff)
    length(mutation.high.diff)
    # 42

    mutation.high.deg = data.frame(mutation.high.mean1, mutation.high.mean2, mutation.high.diff, mutation.high.p, mutation.high.q, row.names=names(mutation.high.diff))
    colnames(mutation.high.deg) = c("mutation.high", "mutation.other", "diff.high2other", "pvalue", "qvalue")
    head(mutation.high.deg)
    dim(mutation.high.deg)
    # 42  5
    # write.csv(mutation.high.deg, paste0("mutation/dep_", signature.name, "_hml", cutoff.percentile, "_mutation.high.deg.csv"), quote=FALSE)

    # signature.low
    mutation.low.p = apply(mutation.signature[1:(ncol(mutation.signature)-1)], 2, function(x) signif(t.test(x[1:(nrow(mutation.signature)-nrow(mutation.low))], x[(nrow(mutation.signature)-nrow(mutation.low)+1):nrow(mutation.signature)])$p.value,5))
    head(mutation.low.p)
    length(mutation.low.p)
    # 42
    mutation.low.q = signif(p.adjust(mutation.low.p, "BH"),5)
    head(mutation.low.q)
    length(mutation.low.q)
    # 42
    mutation.low.mean1 = apply(mutation.signature[1:(ncol(mutation.signature)-1)],  2, function(x) signif(mean(na.omit(x[(nrow(mutation.signature)-nrow(mutation.low)+1):nrow(mutation.signature)])), 5))
    mutation.low.mean2 = apply(mutation.signature[1:(ncol(mutation.signature)-1)],  2, function(x) signif(mean(na.omit(x[1:(nrow(mutation.signature)-nrow(mutation.low))])), 5))
    mutation.low.diff = apply(mutation.signature[1:(ncol(mutation.signature)-1)], 2, function(x) signif(mean(na.omit(x[(nrow(mutation.signature)-nrow(mutation.low)+1):nrow(mutation.signature)])) - mean(na.omit(x[1:(nrow(mutation.signature)-nrow(mutation.low))])), 5))
    head(mutation.low.diff)
    length(mutation.low.diff)
    # 42

    mutation.low.deg = data.frame(mutation.low.mean1, mutation.low.mean2, mutation.low.diff, mutation.low.p, mutation.low.q, row.names=names(mutation.low.diff))
    colnames(mutation.low.deg) = c("mutation.low", "mutation.other", "diff.low2other", "pvalue", "qvalue")
    head(mutation.low.deg)
    dim(mutation.low.deg)
    # 42  5
    # write.csv(mutation.low.deg, paste0("mutation/dep_", signature.name, "_hml", cutoff.percentile, "_mutation.low.deg.csv"), quote=FALSE)

    # signature.high2low
    mutation.high2low.p = apply(mutation.signature[1:(ncol(mutation.signature)-1)], 2, function(x) signif(t.test(na.omit(x[1:nrow(mutation.high)]), na.omit(x[(nrow(mutation.signature)-nrow(mutation.low)+1):nrow(mutation.signature)]))$p.value,5))
    head(mutation.high2low.p)
    length(mutation.high2low.p)
    # 42
    mutation.high2low.q = signif(p.adjust(mutation.high2low.p, "BH"),5)
    head(mutation.high2low.q)
    length(mutation.high2low.q)
    # 42
    mutation.high2low.mean1 = apply(mutation.signature[1:(ncol(mutation.signature)-1)],  2, function(x) signif(mean(na.omit(x[1:nrow(mutation.high)])), 5))
    mutation.high2low.mean2 = apply(mutation.signature[1:(ncol(mutation.signature)-1)],  2, function(x) signif(mean(na.omit(x[(nrow(mutation.signature)-nrow(mutation.low)+1):nrow(mutation.signature)])), 5))
    mutation.high2low.diff = mutation.high2low.mean1 - mutation.high2low.mean2
    head(mutation.high2low.diff)
    length(mutation.high2low.diff)
    # 42

    mutation.high2low.deg = data.frame(mutation.high2low.mean1, mutation.high2low.mean2, mutation.high2low.diff, mutation.high2low.p, mutation.high2low.q, row.names=names(mutation.high2low.diff))
    colnames(mutation.high2low.deg) = c("mutation.high", "mutation.low", "diff.high2low", "pvalue", "qvalue")
    head(mutation.high2low.deg)
    dim(mutation.high2low.deg)
    # 42  5
    # write.csv(mutation.high2low.deg, paste0("mutation/dep_", signature.name, "_hml", cutoff.percentile, "_mutation.high2low.deg.csv"), quote=FALSE)

    # mutation.high.deg = mutation.high.deg[mutation.high.deg$mutation.high > 0 & mutation.high.deg$mutation.other > 0,,drop=FALSE]
    # mutation.low.deg = mutation.low.deg[mutation.low.deg$mutation.low > 0 & mutation.low.deg$mutation.other > 0,,drop=FALSE]

    # cutoff.pvalue = 0.05
    # cutoff.diff = 0.1
    set.seed(42)
    p1 <- ggplot(data = mutation.high.deg, mapping = aes(x = diff.high2other, y = (-1)*log(pvalue, 10))) + 
    geom_point(size=1, color= ifelse(mutation.high.deg$pvalue < cutoff.pvalue & mutation.high.deg$diff.high2other > cutoff.diff, "red", ifelse(mutation.high.deg$pvalue < cutoff.pvalue & mutation.high.deg$diff.high2other < cutoff.diff*(-1), "blue", "grey60")))+ 
    xlim(-0.2,0.2) + 
    # ylim(0,4) + 
    geom_hline(yintercept = (-1)*log(cutoff.pvalue, 10), linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = (-1)*cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    annotate("text", x=-0.2, y=max((-1)*log(mutation.high.deg$pvalue, 10))*1.1, parse=FALSE, label = paste0("Signature.high cell lines: ", nrow(mutation.signature.high)), color = "red", hjust = 0) + 
    geom_label_repel(aes(label=ifelse(mutation.high.deg$pvalue < cutoff.pvalue & abs(mutation.high.deg$diff.high2other) > cutoff.diff, as.character(rownames(mutation.high.deg)), "")), size = 2, color = ifelse(mutation.high.deg$diff.high2other > 0, "red", "blue"), segment.size=0.2) +
    labs(x="Mutation difference", y="-log10(p value)", title=paste0("Signature [", signature.name, "] High vs. others"))
    p1 <- p1 + theme_classic() + rremove("legend")
    # p1

    p2 <- ggplot(data = mutation.low.deg, mapping = aes(x = diff.low2other, y = (-1)*log(pvalue, 10))) + 
    geom_point(size=1, color= ifelse(mutation.low.deg$pvalue < cutoff.pvalue & mutation.low.deg$diff.low2other > cutoff.diff, "red", ifelse(mutation.low.deg$pvalue < cutoff.pvalue & mutation.low.deg$diff.low2other < cutoff.diff*(-1), "blue", "grey60")))+ 
    xlim(-0.2,0.2) + 
    # ylim(0,4) + 
    geom_hline(yintercept = (-1)*log(cutoff.pvalue, 10), linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = (-1)*cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    annotate("text", x=-0.2, y=max((-1)*log(mutation.low.deg$pvalue, 10))*1.1, parse=FALSE, label = paste0("Signature.low cell lines: ", nrow(mutation.signature.low)), color = "red", hjust = 0) + 
    geom_label_repel(aes(label=ifelse(mutation.low.deg$pvalue < cutoff.pvalue & abs(mutation.low.deg$diff.low2other) > cutoff.diff, as.character(rownames(mutation.low.deg)), "")), size = 2, color = ifelse(mutation.low.deg$diff.low2other > 0, "red", "blue"), segment.size=0.2) +
    labs(x="Mutation difference", y="-log10(p value)", title=paste0("Signature [", signature.name, "] Low vs. others"))
    p2 <- p2 + theme_classic() + rremove("legend")
    # p2

    p3 <- ggplot(data = mutation.high2low.deg, mapping = aes(x = diff.high2low, y = (-1)*log(pvalue, 10))) + 
    geom_point(size=1, color= ifelse(mutation.high2low.deg$pvalue < cutoff.pvalue & mutation.high2low.deg$diff.high2low > cutoff.diff, "red", ifelse(mutation.high2low.deg$pvalue < cutoff.pvalue & mutation.high2low.deg$diff.high2low < cutoff.diff*(-1), "blue", "grey60")))+ 
    xlim(-0.2,0.2) + 
    # ylim(0,4) + 
    geom_hline(yintercept = (-1)*log(cutoff.pvalue, 10), linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = (-1)*cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    annotate("text", x=-0.2, y=max((-1)*log(mutation.high2low.deg$pvalue, 10))*1.1, parse=FALSE, label = paste0("Signature.high/low cell lines: ", nrow(mutation.signature.high), "/", nrow(mutation.signature.low)), color = "red", hjust = 0) + 
    geom_label_repel(aes(label=ifelse(mutation.high2low.deg$pvalue < cutoff.pvalue & abs(mutation.high2low.deg$diff.high2low) > cutoff.diff, as.character(rownames(mutation.high2low.deg)), "")), size = 2, color = ifelse(mutation.high2low.deg$diff.high2low > 0, "red", "blue"), segment.size=0.2) +
    labs(x="Mutation difference", y="-log10(p value)", title=paste0("Signature [", signature.name, "] High vs. Low"))
    p3 <- p3 + theme_classic() + rremove("legend")
    # p3

    # Arranging the plot using cowplot
    p = suppressWarnings(plot_grid(p1, p2, p3, ncol = 3, align = "hv", rel_widths = c(1,1), rel_heights = c(1,1)))
    # p <- ggplotGrob(p)
    return(p)
}

