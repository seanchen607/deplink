##' 'cancertypeLow' displays the cancer type component of cell lines with low dependencies of a gene set (signature).
##'
##'
##' @title cancertypeLow
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
##' cancertypeLow(signature.name, signature)


## Main
cancertypeLow <- function(signature.name, 
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

    # Cancer type information
    meta.signature.high = meta[rownames(meta) %in% rownames(dep.t.signature.high),,drop=FALSE]
    head(meta.signature.high)
    dim(meta.signature.high)
    # 56  8
    meta.signature.low = meta[rownames(meta) %in% rownames(dep.t.signature.low),,drop=FALSE]
    head(meta.signature.low)
    dim(meta.signature.low)
    # 56  8
    # write.csv(meta.signature.high, paste0("meta/dep_", signature.name, "_score.high_meta.csv"))
    # write.csv(meta.signature.low, paste0("meta/dep_", signature.name, "_score.low_meta.csv"))

    # Pie Chart from data frame with Appended Sample Sizes
    # low
    mytable <- table(meta.signature.low$disease)
    lbls <- paste(names(mytable), " : ", mytable, sep="")
    lbls = gsub(".+ : 0", NA, as.list(lbls))
    color = rainbow(length(mytable))
    names(color) = names(mytable)
    
    p = ggplot(data = as.data.frame(mytable[mytable>0]), aes(x = "", y = Freq, fill = Var1))+
    geom_bar(stat = "identity", col = "grey30")+
    coord_polar("y", start = 0) +
    geom_text(aes(label = Freq), position = position_stack(vjust = 0.5), col = "black") +
    theme_void() + 
    labs(fill = "Cancer Type") +
    scale_fill_manual(values = color)

    return(p)
}

