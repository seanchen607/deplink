##' 'cancertypeLandscape' displays the landscape of cancer type component of cell lines with different dependencies of a gene set (signature).
##'
##'
##' @title cancertypeLandscape
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
##' cancertypeLandscape(signature.name, signature)


## Main
cancertypeLandscape <- function(signature.name, 
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

    # dot plot
    dep.t.signature.meta = merge(dep.t.signature, meta, by="row.names", all=FALSE)
    rownames(dep.t.signature.meta) = dep.t.signature.meta[,1]
    dep.t.signature.meta = dep.t.signature.meta[,-1]
    head(dep.t.signature.meta)
    dim(dep.t.signature.meta)
    # 558  20
    # write.csv(dep.t.signature.meta, paste0("meta/dep_", signature.name, "_score.meta.csv"))

    dep.t.signature.meta = merge(dep.t.signature.meta, TCGA.tumor, by="row.names", all=FALSE)
    rownames(dep.t.signature.meta) = dep.t.signature.meta[,1]
    dep.t.signature.meta = dep.t.signature.meta[,-1]
    head(dep.t.signature.meta)
    dim(dep.t.signature.meta)
    # 517  21

    dep.t.signature.meta.order = dep.t.signature.meta[order(dep.t.signature.meta$signature.score),]
    head(dep.t.signature.meta.order)
    dim(dep.t.signature.meta.order)
    # 517  21
    # write.csv(dep.t.signature.meta.order, paste0("meta/dep_", signature.name, "_score_CancerType.TCGA.csv"))

    # Primary.Disease.freq = table(dep.t.signature.meta.order$disease)
    Primary.Disease.freq = table(dep.t.signature.meta.order$tcga_code)
    Primary.Disease.freq.cutoff = Primary.Disease.freq[Primary.Disease.freq > cutoff.freq]
    colors = as.list(wes_palette(length(Primary.Disease.freq.cutoff), name = "Darjeeling1", type = "continuous"))
    plotlist = list()
    for (i in 1:length(Primary.Disease.freq.cutoff)) {
        cancer.type = names(Primary.Disease.freq.cutoff)[i]
        # dep.t.signature.meta.order.subset = dep.t.signature.meta.order[dep.t.signature.meta.order$disease == cancer.type,,drop=FALSE]
        dep.t.signature.meta.order.subset = dep.t.signature.meta.order[dep.t.signature.meta.order$tcga_code == cancer.type,,drop=FALSE]
        dep.t.signature.meta.order.subset$order = seq(1:nrow(dep.t.signature.meta.order.subset))
        head(dep.t.signature.meta.order.subset)
        dim(dep.t.signature.meta.order.subset)
        # 558  20
        cancer.type.name = gsub(" .+", "", cancer.type)
        cancer.type.name = gsub("\\/.+", "", cancer.type.name)
        # p = paste0("p", i)
        if (i == 1) {
            p = ggplot(data = dep.t.signature.meta.order.subset, mapping = aes(x = order, y = signature.score)) + 
            geom_point(size=0.5, color= unlist(colors)[i])+ 
            # xlim(-1,1) + 
            ylim(min(dep.t.signature.meta$signature.score),max(dep.t.signature.meta$signature.score)) + 
            # geom_hline(yintercept = min(dep.t.signature.high$signature.score), linetype="dashed", colour="grey30", size=0.2) +
            # geom_hline(yintercept = max(dep.t.signature.low$signature.score), linetype="dashed", colour="grey30", size=0.2) +
            labs(x= cancer.type.name, y="Signature score")+
            theme_classic() + rremove("legend") +
            theme(axis.title.x=element_text(size=9), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_text(size=9, color=NA), axis.line.x=element_blank(), axis.text.y=element_text(size=9, color=NA), axis.line.y=element_line(color=NA))
        } else {
            p = ggplot(data = dep.t.signature.meta.order.subset, mapping = aes(x = order, y = signature.score)) + 
            geom_point(size=0.5, color= unlist(colors)[i])+ 
            # xlim(-1,1) + 
            ylim(min(dep.t.signature.meta$signature.score),max(dep.t.signature.meta$signature.score)) + 
            # geom_hline(yintercept = min(dep.t.signature.high$signature.score), linetype="dashed", colour="grey30", size=0.2) +
            # geom_hline(yintercept = max(dep.t.signature.low$signature.score), linetype="dashed", colour="grey30", size=0.2) +
            labs(x= cancer.type.name, y="Signature score")+
            theme_classic() + rremove("legend") +
            theme(axis.title.x=element_text(size=9), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.line.y=element_blank())
        }
        plotlist[[i]] = p
    }
    # Arranging the plot using cowplot
    # paste(as.list(paste0("p", seq(1:length(Primary.Disease.freq.cutoff)))), collapse = ",")
    # plotlist = map(paste0("p", seq(1:length(Primary.Disease.freq.cutoff))), get)
    p = ggarrange(plotlist = plotlist, ncol = length(Primary.Disease.freq.cutoff), nrow = 1, widths = c(2, rep(1,length(Primary.Disease.freq.cutoff)-1)), heights = c(1,1,1,1,1))

    p0 = ggplot(data = dep.t.signature.meta.order.subset, mapping = aes(x = order, y = signature.score)) + 
    geom_point(size=NA, color= unlist(colors)[i])+ 
    # xlim(-1,1) + 
    ylim(min(dep.t.signature.meta$signature.score),max(dep.t.signature.meta$signature.score)) + 
    geom_hline(yintercept = min(dep.t.signature.high$signature.score), linetype="dashed", colour="grey30", size=0.2) +
    geom_hline(yintercept = max(dep.t.signature.low$signature.score), linetype="dashed", colour="grey30", size=0.2) +
    labs(x= cancer.type.name, y="Signature score")+
    theme_classic() + rremove("legend") +
    theme(axis.title.x=element_text(size=9, color=NA), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_text(size=9), axis.line.x=element_line(), axis.text.y=element_text(size=9), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA))

    p = suppressWarnings(p + annotation_custom(grob = ggplotGrob(p0)))
    return(p)
}

