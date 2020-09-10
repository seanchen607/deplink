##' 'deplink' compares the genetic/epigenetic features between cancer cell lines with different dependencies of a gene set (signature).
##'
##'
##' @title deplink
##' @param signature.name Names of a signature (format: character)
##' @param signature Gene names of a signature (format: vector)
##' @param outputDir Output directory, default root directory (format: character)
##' @param cutoff.freq Cutoff for frequency of cancer cell lines for each cancer type, default 10
##' @param cutoff.percentile Cutoff for percentile of cancer cell lines with highest/lowest dependency, default 0.2
##' @param cutoff.pvalue Cutoff for p-value of the T-test results, default 0.05
##' @param cutoff.qvalue Cutoff for q-value of the T-test results, default 0.1
##' @param cutoff.diff Cutoff for difference of the T-test results, default 0.1
##' @param cutoff.fc Cutoff for fold changes of the T-test results, default 2
##' @return message
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
##' deplink(signature.name, signature)


## Main
deplink <- function(signature.name, 
                    signature, 
                    outputDir          = "~", 
                    cutoff.freq        = 10, 
                    cutoff.percentile  = 0.2, 
                    cutoff.pvalue      = 0.05, 
                    cutoff.qvalue      = 0.1, 
                    cutoff.diff        = 0.1,
                    cutoff.fc          = 2
                    ) {
    message("*** deplink analyzing ...")
    setwd(file.path(outputDir))

    #############################################

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

    dir.create(file.path(outputDir, signature.name), showWarnings = FALSE)
    setwd(file.path(outputDir, signature.name))

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
    write.csv(dep.t.signature, paste0("dep_", signature.name, "_score.csv"))
    write.csv(dep.t.signature.high, paste0("dep_", signature.name, "_score.high", cutoff.percentile, ".csv"))
    write.csv(dep.t.signature.low, paste0("dep_", signature.name, "_score.low", cutoff.percentile, ".csv"))

    # Cancer type information
    dir.create(file.path(outputDir, signature.name, "meta"), showWarnings = FALSE)

    meta.signature.high = meta[rownames(meta) %in% rownames(dep.t.signature.high),,drop=FALSE]
    head(meta.signature.high)
    dim(meta.signature.high)
    # 56  8
    meta.signature.low = meta[rownames(meta) %in% rownames(dep.t.signature.low),,drop=FALSE]
    head(meta.signature.low)
    dim(meta.signature.low)
    # 56  8
    write.csv(meta.signature.high, paste0("meta/dep_", signature.name, "_score.high_meta.csv"))
    write.csv(meta.signature.low, paste0("meta/dep_", signature.name, "_score.low_meta.csv"))

    # Pie Chart from data frame with Appended Sample Sizes
    # high
    mytable <- table(meta.signature.high$disease)
    lbls <- paste(names(mytable), " : ", mytable, sep="")
    lbls = gsub(".+ : 0", NA, as.list(lbls))
    color = rainbow(length(mytable))
    names(color) = names(mytable)
    
    p1 = ggplot(data = as.data.frame(mytable[mytable>0]), aes(x = "", y = Freq, fill = Var1))+
    geom_bar(stat = "identity", col = "grey30")+
    coord_polar("y", start = 0) +
    geom_text(aes(label = Freq), position = position_stack(vjust = 0.5), col = "black") +
    theme_void() + 
    labs(fill = "Cancer Type") +
    scale_fill_manual(values = color)
    ggsave(paste0("meta/dep_", signature.name, "_score.high_pie.Primary.pdf"), p1, width=8, height=6)

    # low
    mytable <- table(meta.signature.low$disease)
    lbls <- paste(names(mytable), " : ", mytable, sep="")
    lbls = gsub(".+ : 0", NA, as.list(lbls))
    color = rainbow(length(mytable))
    names(color) = names(mytable)
    
    p2 = ggplot(data = as.data.frame(mytable[mytable>0]), aes(x = "", y = Freq, fill = Var1))+
    geom_bar(stat = "identity", col = "grey30")+
    coord_polar("y", start = 0) +
    geom_text(aes(label = Freq), position = position_stack(vjust = 0.5), col = "black") +
    theme_void() + 
    labs(fill = "Cancer Type") +
    scale_fill_manual(values = color)
    ggsave(paste0("meta/dep_", signature.name, "_score.low_pie.Primary.pdf"), p2, width=8, height=6)

    primary.pair = intersect(meta.signature.high$disease, meta.signature.low$disease)
    head(primary.pair)
    length(primary.pair)
    # 14
    meta.signature.high.pair = meta.signature.high[meta.signature.high$disease %in% primary.pair & !meta.signature.high$disease %in% c("Leukemia", "Lymphoma", "Myeloma"),,drop=FALSE]
    meta.signature.high.pair = meta.signature.high.pair[order(meta.signature.high.pair$disease),,drop=FALSE]
    head(meta.signature.high.pair)
    dim(meta.signature.high.pair)
    # 32  8
    meta.signature.low.pair = meta.signature.low[meta.signature.low$disease %in% primary.pair & !meta.signature.low$disease %in% c("Leukemia", "Lymphoma", "Myeloma"),,drop=FALSE]
    meta.signature.low.pair = meta.signature.low.pair[order(meta.signature.low.pair$disease),,drop=FALSE]
    head(meta.signature.low.pair)
    dim(meta.signature.low.pair)
    # 33  8
    write.csv(meta.signature.high.pair, paste0("meta/dep_", signature.name, "_score.high_meta.pair.csv"))
    write.csv(meta.signature.low.pair, paste0("meta/dep_", signature.name, "_score.low_meta.pair.csv"))

    # dot plot
    dep.t.signature.meta = merge(dep.t.signature, meta, by="row.names", all=FALSE)
    rownames(dep.t.signature.meta) = dep.t.signature.meta[,1]
    dep.t.signature.meta = dep.t.signature.meta[,-1]
    head(dep.t.signature.meta)
    dim(dep.t.signature.meta)
    # 558  20
    write.csv(dep.t.signature.meta, paste0("meta/dep_", signature.name, "_score.meta.csv"))

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
    write.csv(dep.t.signature.meta.order, paste0("meta/dep_", signature.name, "_score_CancerType.TCGA.csv"))

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
    ggsave(paste0("meta/dep_", signature.name, "_score_CancerType.TCGA.dotplot.pdf"), p, width=9, height=2)

    message("[01/14] Analysis-cancer.type: done!")
    #############################################

    # Chromatin modification
    dir.create(file.path(outputDir, signature.name, "chromatin"), showWarnings = FALSE)

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
    write.csv(chromatin.share.signature, paste0("chromatin/dep_", signature.name, "_hml", cutoff.percentile, "_chromatin.csv"), quote=FALSE)

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
    write.csv(chromatin.high.deg, paste0("chromatin/dep_", signature.name, "_hml", cutoff.percentile, "_chromatin.high.deg.csv"), quote=FALSE)

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
    write.csv(chromatin.low.deg, paste0("chromatin/dep_", signature.name, "_hml", cutoff.percentile, "_chromatin.low.deg.csv"), quote=FALSE)

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
    p <- ggplotGrob(p)
    ggsave(paste0("chromatin/dep_", signature.name, "_hml", cutoff.percentile, "_chromatin.high.low.deg.p", cutoff.pvalue, ".pdf"), p, width=8.5, height=4.5)

    message("[02/14] Analysis-chromatin.modification: done!")
    #############################################

    # Mutation - COSMIC
    dir.create(file.path(outputDir, signature.name, "COSMIC"), showWarnings = FALSE)

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
    write.csv(cosmic.share.signature, paste0("cosmic/dep_", signature.name, "_hml", cutoff.percentile, "_cosmic.csv"), quote=FALSE)

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
    write.csv(cosmic.high.deg, paste0("cosmic/dep_", signature.name, "_hml", cutoff.percentile, "_cosmic.high.deg.csv"), quote=FALSE)

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
    write.csv(cosmic.low.deg, paste0("cosmic/dep_", signature.name, "_hml", cutoff.percentile, "_cosmic.low.deg.csv"), quote=FALSE)

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
    p <- ggplotGrob(p)
    ggsave(paste0("cosmic/dep_", signature.name, "_hml", cutoff.percentile, "_cosmic.high.low.deg.pdf"), p, width=8.5, height=4.5)

    message("[03/14] Analysis-cosmic: done!")
    #############################################
 
    # Mutation
    dir.create(file.path(outputDir, signature.name, "mutation"), showWarnings = FALSE)

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
    write.csv(mutation.high.deg, paste0("mutation/dep_", signature.name, "_hml", cutoff.percentile, "_mutation.high.deg.csv"), quote=FALSE)

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
    write.csv(mutation.low.deg, paste0("mutation/dep_", signature.name, "_hml", cutoff.percentile, "_mutation.low.deg.csv"), quote=FALSE)

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
    write.csv(mutation.high2low.deg, paste0("mutation/dep_", signature.name, "_hml", cutoff.percentile, "_mutation.high2low.deg.csv"), quote=FALSE)

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
    p <- ggplotGrob(p)
    ggsave(paste0("mutation/dep_", signature.name, "_hml", cutoff.percentile, "_mutation.high.low.deg.p", cutoff.pvalue,".pdf"), p, width=12.5, height=4.5)

    message("[04/14] Analysis-mutation: done!")
    #############################################
 
    # Drug sensitivity - GDSC
    dir.create(file.path(outputDir, signature.name, "drug.GDSC"), showWarnings = FALSE)

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
    write.csv(drug.high.deg, paste0("drug.GDSC/dep_", signature.name, "_hml", cutoff.percentile, "_drug.high.deg.csv"), quote=TRUE)

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
    write.csv(drug.low.deg, paste0("drug.GDSC/dep_", signature.name, "_hml", cutoff.percentile, "_drug.low.deg.csv"), quote=TRUE)

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
    p <- ggplotGrob(p)
    ggsave(paste0("drug.GDSC/dep_", signature.name, "_hml", cutoff.percentile, "_drug.high.low.deg.pdf"), p, width=8.5, height=8.5, limitsize=FALSE, device="pdf")

    message("[05/14] Analysis-drug.GDSC: done!")
    #############################################

    # Drug sensitivity - PRISM
    dir.create(file.path(outputDir, signature.name, "drug.PRISM"), showWarnings = FALSE)
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
    write.csv(drug2.high.deg, paste0("drug.PRISM/dep_", signature.name, "_hml", cutoff.percentile, "_drug.high.deg.csv"), quote=TRUE)

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
    write.csv(drug2.low.deg, paste0("drug.PRISM/dep_", signature.name, "_hml", cutoff.percentile, "_drug.low.deg.csv"), quote=TRUE)

    drug2.high.deg$moa = gsub(",\\(.+\\)", "", drug2.high.deg$moa)
    drug2.low.deg$moa = gsub(",\\(.+\\)", "", drug2.low.deg$moa)

    # cutoff.pvalue = 0.01
    # cutoff.diff = 0.2
    set.seed(42)
    p1 <- ggplot(data = drug2.high.deg, mapping = aes(x = diff.high2other, y = (-1)*log(pvalue, 10))) + 
    geom_point(size=1, color= ifelse(drug2.high.deg$pvalue < cutoff.pvalue & drug2.high.deg$diff.high2other > 3*cutoff.diff, "blue", ifelse(drug2.high.deg$pvalue < cutoff.pvalue & drug2.high.deg$diff.high2other < 3*cutoff.diff*(-1), "red", "grey60")))+ 
    # xlim(-1.5,1.5) + 
    # ylim(0,4) + 
    geom_hline(yintercept = (-1)*log(cutoff.pvalue, 10), linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = 3*cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = (-1)*3*cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    annotate("text", x=min(na.omit(drug2.high.deg$diff.high2other), na.omit(drug2.high.deg$diff.low2other)), y=max(na.omit((-1)*log(drug2.high.deg$pvalue, 10)), na.omit((-1)*log(drug2.low.deg$pvalue, 10)))*1.1, parse=FALSE, label = paste0("Signature.high cell lines: ", nrow(drug2.share.signature.high)), color = "red", hjust = 0) + 
    geom_label_repel(aes(label=ifelse(drug2.high.deg$pvalue < cutoff.pvalue & abs(drug2.high.deg$diff.high2other) > 3*cutoff.diff, as.character(drug2.high.deg$name), "")), size = 2, color = ifelse(drug2.high.deg$diff.high2other > 0, "blue", "red"), segment.size=0.2) +
    labs(x="Sensitivity difference", y="-log10(p value)", title=paste0("Signature [", signature.name, "] High vs. others"))
    p1 <- p1 + theme_classic() + rremove("legend")
    # p1

    p2 <- ggplot(data = drug2.low.deg, mapping = aes(x = diff.low2other, y = (-1)*log(pvalue, 10))) + 
    geom_point(size=1, color= ifelse(drug2.low.deg$pvalue < cutoff.pvalue & drug2.low.deg$diff.low2other > 3*cutoff.diff, "blue", ifelse(drug2.low.deg$pvalue < cutoff.pvalue & drug2.low.deg$diff.low2other < 3*cutoff.diff*(-1), "red", "grey60")))+ 
    # xlim(-1.5,1.5) + 
    # ylim(0,4) + 
    geom_hline(yintercept = (-1)*log(cutoff.pvalue, 10), linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = 3*cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = (-1)*3*cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    annotate("text", x=min(na.omit(drug2.high.deg$diff.high2other), na.omit(drug2.high.deg$diff.low2other)), y=max(na.omit((-1)*log(drug2.high.deg$pvalue, 10)), na.omit((-1)*log(drug2.low.deg$pvalue, 10)))*1.1, parse=FALSE, label = paste0("Signature.low cell lines: ", nrow(drug2.share.signature.low)), color = "red", hjust = 0) + 
    geom_label_repel(aes(label=ifelse(drug2.low.deg$pvalue < cutoff.pvalue & abs(drug2.low.deg$diff.low2other) > 3*cutoff.diff, as.character(drug2.low.deg$name), "")), size = 2, color = ifelse(drug2.low.deg$diff.low2other > 0, "blue", "red"), segment.size=0.2) +
    labs(x="Sensitivity difference", y="-log10(p value)", title=paste0("Signature [", signature.name, "] Low vs. others"))
    p2 <- p2 + theme_classic() + rremove("legend")
    # p2

    p3 <- ggplot(data = drug2.high.deg, mapping = aes(x = diff.high2other, y = (-1)*log(pvalue, 10))) + 
    geom_point(size=1, color= ifelse(drug2.high.deg$pvalue < cutoff.pvalue & drug2.high.deg$diff.high2other > 3*cutoff.diff, "blue", ifelse(drug2.high.deg$pvalue < cutoff.pvalue & drug2.high.deg$diff.high2other < 3*cutoff.diff*(-1), "red", "grey60")))+ 
    # xlim(-1.5,1.5) + 
    # ylim(0,4) + 
    geom_hline(yintercept = (-1)*log(cutoff.pvalue, 10), linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = 3*cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = (-1)*3*cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    annotate("text", x=min(na.omit(drug2.high.deg$diff.high2other), na.omit(drug2.high.deg$diff.low2other)), y=max(na.omit((-1)*log(drug2.high.deg$pvalue, 10)), na.omit((-1)*log(drug2.low.deg$pvalue, 10)))*1.1, parse=FALSE, label = paste0("Signature.high cell lines: ", nrow(drug2.share.signature.high)), color = "red", hjust = 0) + 
    geom_label_repel(aes(label=ifelse(drug2.high.deg$pvalue < cutoff.pvalue & abs(drug2.high.deg$diff.high2other) > 3*cutoff.diff, as.character(drug2.high.deg$moa), "")), size = 2, color = ifelse(drug2.high.deg$diff.high2other > 0, "blue", "red"), segment.size=0.2) +
    labs(x="Sensitivity difference", y="-log10(p value)", title=paste0("Signature [", signature.name, "] High vs. others"))
    p3 <- p3 + theme_classic() + rremove("legend")
    # p3

    p4 <- ggplot(data = drug2.low.deg, mapping = aes(x = diff.low2other, y = (-1)*log(pvalue, 10))) + 
    geom_point(size=1, color= ifelse(drug2.low.deg$pvalue < cutoff.pvalue & drug2.low.deg$diff.low2other > 3*cutoff.diff, "blue", ifelse(drug2.low.deg$pvalue < cutoff.pvalue & drug2.low.deg$diff.low2other < 3*cutoff.diff*(-1), "red", "grey60")))+ 
    # xlim(-1.5,1.5) + 
    # ylim(0,4) + 
    geom_hline(yintercept = (-1)*log(cutoff.pvalue, 10), linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = 3*cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    geom_vline(xintercept = (-1)*3*cutoff.diff, linetype="dashed", colour="grey30", size=0.2) + 
    annotate("text", x=min(na.omit(drug2.high.deg$diff.high2other), na.omit(drug2.high.deg$diff.low2other)), y=max(na.omit((-1)*log(drug2.high.deg$pvalue, 10)), na.omit((-1)*log(drug2.low.deg$pvalue, 10)))*1.1, parse=FALSE, label = paste0("Signature.low cell lines: ", nrow(drug2.share.signature.low)), color = "red", hjust = 0) + 
    geom_label_repel(aes(label=ifelse(drug2.low.deg$pvalue < cutoff.pvalue & abs(drug2.low.deg$diff.low2other) > 3*cutoff.diff, as.character(drug2.low.deg$moa), "")), size = 2, color = ifelse(drug2.low.deg$diff.low2other > 0, "blue", "red"), segment.size=0.2) +
    labs(x="Sensitivity difference", y="-log10(p value)", title=paste0("Signature [", signature.name, "] Low vs. others"))
    p4 <- p4 + theme_classic() + rremove("legend")
    # p4

    # Arranging the plot using cowplot
    p = suppressWarnings(plot_grid(p1, p2, p3, p4, ncol = 2, align = "hv", rel_widths = c(1,1), rel_heights = c(1,1)))
    p <- ggplotGrob(p)
    ggsave(paste0("drug.PRISM/dep_", signature.name, "_hml", cutoff.percentile, "_drug.high.low.deg.pdf"), p, width=8.5, height=8.5, limitsize=FALSE, device="pdf")

    message("[06/14] Analysis-drug.PRISM: done!")
    #############################################

    # Dependency
    dir.create(file.path(outputDir, signature.name, "dependency"), showWarnings = FALSE)
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
    write.csv(dependency.high.deg, paste0("dependency/dep_", signature.name, "_hml", cutoff.percentile, "_dependency.high.deg.csv"), quote=TRUE)

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
    write.csv(dependency.low.deg, paste0("dependency/dep_", signature.name, "_hml", cutoff.percentile, "_dependency.low.deg.csv"), quote=TRUE)

    # cutoff.qvalue = 0.1
    # cutoff.diff = 0.1
    set.seed(42)
    p1 <- ggplot(data = dependency.high.deg, mapping = aes(x = diff.high2other, y = (-1)*log(qvalue, 10))) + 
    geom_point(size=1, color= ifelse(dependency.high.deg$qvalue < cutoff.qvalue & dependency.high.deg$diff.high2other > cutoff.diff, "blue", ifelse(dependency.high.deg$qvalue < cutoff.qvalue & dependency.high.deg$diff.high2other < cutoff.diff*(-1), "red", "grey60")))+ 
    xlim(-0.5,0.5) + 
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
    xlim(-0.5,0.5) + 
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
    p <- ggplotGrob(p)
    ggsave(paste0("dependency/dep_", signature.name, "_hml", cutoff.percentile, "_dependency.high.low.deg.q", cutoff.qvalue, ".pdf"), p, width=8.5, height=4.5, limitsize=FALSE, device="pdf")

    message("[07/14] Analysis-dependency: done!")
    #############################################
 
    # Expression
    dir.create(file.path(outputDir, signature.name, "expression"), showWarnings = FALSE)
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
    write.csv(expression.high.deg, paste0("expression/dep_", signature.name, "_hml", cutoff.percentile, "_expression.high.deg.csv"), quote=TRUE)

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
    write.csv(expression.low.deg, paste0("expression/dep_", signature.name, "_hml", cutoff.percentile, "_expression.low.deg.csv"), quote=TRUE)

    expression.high.deg = read.csv(paste0("expression/dep_", signature.name, "_hml", cutoff.percentile, "_expression.high.deg.csv"), header=TRUE, row.names=1)
    expression.low.deg = read.csv(paste0("expression/dep_", signature.name, "_hml", cutoff.percentile, "_expression.low.deg.csv"), header=TRUE, row.names=1)
    cutoff.qvalue = 0.1
    cutoff.fc = 2
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
    p <- ggplotGrob(p)
    ggsave(paste0("expression/dep_", signature.name, "_hml", cutoff.percentile, "_expression.high.low.deg.q", cutoff.qvalue, ".pdf"), p, width=8.5, height=4.5, limitsize=FALSE, device="pdf")

    message("[08/14] Analysis-expression: done!")
    #############################################
 
    # Genome instability
    dir.create(file.path(outputDir, signature.name, "genome.instability"), showWarnings = FALSE)
    # Tumor Mutation burden
    dir.create(file.path(outputDir, signature.name, "genome.instability/TMB"), showWarnings = FALSE)

    TMB.share = TMB[rownames(TMB) %in% rownames(dep.t.signature),,drop=FALSE]
    head(TMB.share)
    dim(TMB.share)
    # 554   1

    TMB.share.signature.high = TMB.share[rownames(TMB.share) %in% rownames(dep.t.signature.high),,drop=FALSE]
    TMB.share.signature.high$signature = paste0(signature.name, ".dep.high")
    head(TMB.share.signature.high)
    dim(TMB.share.signature.high)
    # 56  2
    TMB.share.signature.low = TMB.share[rownames(TMB.share) %in% rownames(dep.t.signature.low),,drop=FALSE]
    TMB.share.signature.low$signature = paste0(signature.name, ".dep.low")
    head(TMB.share.signature.low)
    dim(TMB.share.signature.low)
    # 56  2
    TMB.share.signature.mid = TMB.share[!rownames(TMB.share) %in% c(rownames(dep.t.signature.high), rownames(dep.t.signature.low)),,drop=FALSE]
    TMB.share.signature.mid$signature = paste0(signature.name, ".dep.mid")
    head(TMB.share.signature.mid)
    dim(TMB.share.signature.mid)
    # 442   2

    TMB.share.signature = rbind(TMB.share.signature.high, TMB.share.signature.mid, TMB.share.signature.low)
    head(TMB.share.signature)
    dim(TMB.share.signature)
    # 554   2
    write.csv(TMB.share.signature, paste0("genome.instability/TMB/dep_", signature.name, "_hml", cutoff.percentile, "_TMB.csv"), quote=TRUE)
    # write.csv(TMB.share.signature, paste0("genome.instability/TMB/dep_", signature.name, "_hml", cutoff.percentile, "_TMB.sig.csv"), quote=TRUE)

    group.high = TMB.share.signature[TMB.share.signature$signature %like% "high",,drop=FALSE]
    group.mid = TMB.share.signature[TMB.share.signature$signature %like% "mid",,drop=FALSE]
    group.low = TMB.share.signature[TMB.share.signature$signature %like% "low",,drop=FALSE]
    # t.test(group1[,1], group2[,1])$p.value
    # 0.001361099
    p.hm = paste("p.hm = ", signif(t.test(group.high$TMB, group.mid$TMB)$p.value, 4))
    p.hl = paste("p.hl = ", signif(t.test(group.high$TMB, group.low$TMB)$p.value, 4))
    p.ml = paste("p.ml = ", signif(t.test(group.mid$TMB, group.low$TMB)$p.value, 4))
    TMB.share.signature$signature = gsub(paste0(signature.name, "\\.dep\\.high"), "1.high", TMB.share.signature$signature)
    TMB.share.signature$signature = gsub(paste0(signature.name, "\\.dep\\.mid"), "2.mid", TMB.share.signature$signature)
    TMB.share.signature$signature = gsub(paste0(signature.name, "\\.dep\\.low"), "3.low", TMB.share.signature$signature)

    p1 <- ggplot(TMB.share.signature, aes(x=signature, y= log10(TMB), fill=signature)) + 
      geom_violin(trim=FALSE, linetype="blank", na.rm=TRUE)+
      geom_boxplot(width=0.05, fill="white", outlier.size=0.1, na.rm=TRUE)+
      annotate("text", x=0.5, y=log10(max(na.omit(TMB.share.signature$TMB))*1.5), parse=FALSE, hjust=0, label = p.hm)+
      annotate("text", x=0.5, y=log10(max(na.omit(TMB.share.signature$TMB))*1.1), parse=FALSE, hjust=0, label = p.hl)+
      annotate("text", x=0.5, y=log10(max(na.omit(TMB.share.signature$TMB))*0.8), parse=FALSE, hjust=0, label = p.ml)+
      labs(title=paste0("TMB of signature [", signature.name, "]"), x=paste0(signature.name, " score"), y = "TMB (log10)") +  theme_classic() + rremove("legend")
    ggsave(paste0("genome.instability/TMB/dep_", signature.name, "_hml", cutoff.percentile, "_TMB.pdf"), p1, width=3, height=4)

    message("[09/14] Analysis-TMB: done!")
    #############################################

    # CNV
    dir.create(file.path(outputDir, signature.name, "genome.instability/CNV"), showWarnings = FALSE)

    CNV.share = CNV[rownames(CNV) %in% rownames(dep.t.signature),,drop=FALSE]
    head(CNV.share)
    dim(CNV.share)
    # 554   3

    CNV.share.signature.high = CNV.share[rownames(CNV.share) %in% rownames(dep.t.signature.high),,drop=FALSE]
    CNV.share.signature.high$signature = paste0(signature.name, ".dep.high")
    head(CNV.share.signature.high)
    dim(CNV.share.signature.high)
    # 56  2
    CNV.share.signature.low = CNV.share[rownames(CNV.share) %in% rownames(dep.t.signature.low),,drop=FALSE]
    CNV.share.signature.low$signature = paste0(signature.name, ".dep.low")
    head(CNV.share.signature.low)
    dim(CNV.share.signature.low)
    # 56  2
    CNV.share.signature.mid = CNV.share[!rownames(CNV.share) %in% c(rownames(dep.t.signature.high), rownames(dep.t.signature.low)),,drop=FALSE]
    CNV.share.signature.mid$signature = paste0(signature.name, ".dep.mid")
    head(CNV.share.signature.mid)
    dim(CNV.share.signature.mid)
    # 442   2

    CNV.share.signature = rbind(CNV.share.signature.high, CNV.share.signature.mid, CNV.share.signature.low)
    head(CNV.share.signature)
    dim(CNV.share.signature)
    # 554   2
    write.csv(CNV.share.signature, paste0("genome.instability/CNV/dep_", signature.name, "_hml", cutoff.percentile, "_CNV.csv"), quote=TRUE)

    group.high = CNV.share.signature[CNV.share.signature$signature %like% "high",,drop=FALSE]
    group.mid = CNV.share.signature[CNV.share.signature$signature %like% "mid",,drop=FALSE]
    group.low = CNV.share.signature[CNV.share.signature$signature %like% "low",,drop=FALSE]
    # t.test(group1[,1], group2[,1])$p.value
    # 0.001361099
    p.hm = paste("p.hm = ", signif(t.test(group.high$CNV, group.mid$CNV)$p.value, 4))
    p.hl = paste("p.hl = ", signif(t.test(group.high$CNV, group.low$CNV)$p.value, 4))
    p.ml = paste("p.ml = ", signif(t.test(group.mid$CNV, group.low$CNV)$p.value, 4))
    CNV.share.signature$signature = gsub(paste0(signature.name, "\\.dep\\.high"), "1.high", CNV.share.signature$signature)
    CNV.share.signature$signature = gsub(paste0(signature.name, "\\.dep\\.mid"), "2.mid", CNV.share.signature$signature)
    CNV.share.signature$signature = gsub(paste0(signature.name, "\\.dep\\.low"), "3.low", CNV.share.signature$signature)

    p1 <- ggplot(CNV.share.signature, aes(x=signature, y= log10(CNV), fill=signature)) + 
      geom_violin(trim=FALSE, linetype="blank", na.rm=TRUE)+
      geom_boxplot(width=0.05, fill="white", outlier.size=0.1, na.rm=TRUE)+
      annotate("text", x=0.5, y=log10(max(na.omit(CNV.share.signature$CNV))*1.5), parse=FALSE, hjust=0, label = p.hm)+
      annotate("text", x=0.5, y=log10(max(na.omit(CNV.share.signature$CNV))*1.1), parse=FALSE, hjust=0, label = p.hl)+
      annotate("text", x=0.5, y=log10(max(na.omit(CNV.share.signature$CNV))*0.8), parse=FALSE, hjust=0, label = p.ml)+
      labs(title=paste0("CNV of signature [", signature.name, "]"), x=paste0(signature.name, " score"), y = "CNV (log10)") +  theme_classic() + rremove("legend")
    ggsave(paste0("genome.instability/CNV/dep_", signature.name, "_hml", cutoff.percentile, "_CNV.pdf"), p1, width=3, height=4)

    message("[10/14] Analysis-CNV: done!")
    #############################################

    # Microsatellite instability (MSI)
    dir.create(file.path(outputDir, signature.name, "genome.instability/MSI"), showWarnings = FALSE)

    MSI.share = MSI[rownames(MSI) %in% rownames(dep.t.signature),,drop=FALSE]
    head(MSI.share)
    dim(MSI.share)
    # 474  11

    MSI.share.signature.high = MSI.share[rownames(MSI.share) %in% rownames(dep.t.signature.high),,drop=FALSE]
    MSI.share.signature.high$signature = paste0(signature.name, ".dep.high")
    head(MSI.share.signature.high)
    dim(MSI.share.signature.high)
    # 56  2
    MSI.share.signature.low = MSI.share[rownames(MSI.share) %in% rownames(dep.t.signature.low),,drop=FALSE]
    MSI.share.signature.low$signature = paste0(signature.name, ".dep.low")
    head(MSI.share.signature.low)
    dim(MSI.share.signature.low)
    # 56  2
    MSI.share.signature.mid = MSI.share[!rownames(MSI.share) %in% c(rownames(dep.t.signature.high), rownames(dep.t.signature.low)),,drop=FALSE]
    MSI.share.signature.mid$signature = paste0(signature.name, ".dep.mid")
    head(MSI.share.signature.mid)
    dim(MSI.share.signature.mid)
    # 442   2

    MSI.share.signature = rbind(MSI.share.signature.high, MSI.share.signature.mid, MSI.share.signature.low)
    head(MSI.share.signature)
    dim(MSI.share.signature)
    # 554   2
    write.csv(MSI.share.signature, paste0("genome.instability/MSI/dep_", signature.name, "_hml", cutoff.percentile, "_MSI.csv"), quote=TRUE)
    # write.csv(MSI.share.signature, paste0("genome.instability/MSI/dep_", signature.name, "_hml", cutoff.percentile, "_MSI.GDSC.csv"), quote=TRUE)

    group.high = MSI.share.signature[MSI.share.signature$signature %like% "high",,drop=FALSE]
    group.mid = MSI.share.signature[MSI.share.signature$signature %like% "mid",,drop=FALSE]
    group.low = MSI.share.signature[MSI.share.signature$signature %like% "low",,drop=FALSE]
    # t.test(group1[,1], group2[,1])$p.value
    # 0.001361099
    p.hm = paste("p.hm = ", signif(t.test(group.high$MSI, group.mid$MSI)$p.value, 4))
    p.hl = paste("p.hl = ", signif(t.test(group.high$MSI, group.low$MSI)$p.value, 4))
    p.ml = paste("p.ml = ", signif(t.test(group.mid$MSI, group.low$MSI)$p.value, 4))
    MSI.share.signature$signature = gsub(paste0(signature.name, "\\.dep\\.high"), "1.high", MSI.share.signature$signature)
    MSI.share.signature$signature = gsub(paste0(signature.name, "\\.dep\\.mid"), "2.mid", MSI.share.signature$signature)
    MSI.share.signature$signature = gsub(paste0(signature.name, "\\.dep\\.low"), "3.low", MSI.share.signature$signature)

    p1 <- ggplot(MSI.share.signature, aes(x=signature, y= log10(MSI), fill=signature)) + 
      geom_violin(trim=FALSE, linetype="blank", na.rm=TRUE)+
      geom_boxplot(width=0.05, fill="white", outlier.size=0.1, na.rm=TRUE)+
      annotate("text", x=0.5, y=log10(max(na.omit(MSI.share.signature$MSI))*1.5), parse=FALSE, hjust=0, label = p.hm)+
      annotate("text", x=0.5, y=log10(max(na.omit(MSI.share.signature$MSI))*1.1), parse=FALSE, hjust=0, label = p.hl)+
      annotate("text", x=0.5, y=log10(max(na.omit(MSI.share.signature$MSI))*0.8), parse=FALSE, hjust=0, label = p.ml)+
      labs(title=paste0("MSI of signature [", signature.name, "]"), x=paste0(signature.name, " score"), y = "MSI (log10)") +  theme_classic() + rremove("legend")
    # p1
    ggsave(paste0("genome.instability/MSI/dep_", signature.name, "_hml", cutoff.percentile, "_MSI.pdf"), p1, width=3, height=4)
    # ggsave(paste0("genome.instability/MSI/dep_", signature.name, "_hml", cutoff.percentile, "_MSI.GDSC.pdf"), p1, width=3, height=4)

    message("[11/14] Analysis-MSI: done!")
    #############################################
 
    # ISG
    dir.create(file.path(outputDir, signature.name, "ISG"), showWarnings = FALSE)

    ISG.share = ISG[rownames(ISG) %in% rownames(dep.t.signature),,drop=FALSE]
    head(ISG.share)
    dim(ISG.share)
    # 558   9

    ISG.share.signature.high = ISG.share[rownames(ISG.share) %in% rownames(dep.t.signature.high),,drop=FALSE]
    ISG.share.signature.high$signature = paste0(signature.name, ".dep.high")
    head(ISG.share.signature.high)
    dim(ISG.share.signature.high)
    # 56  2
    ISG.share.signature.low = ISG.share[rownames(ISG.share) %in% rownames(dep.t.signature.low),,drop=FALSE]
    ISG.share.signature.low$signature = paste0(signature.name, ".dep.low")
    head(ISG.share.signature.low)
    dim(ISG.share.signature.low)
    # 56  2
    ISG.share.signature.mid = ISG.share[!rownames(ISG.share) %in% c(rownames(dep.t.signature.high), rownames(dep.t.signature.low)),,drop=FALSE]
    ISG.share.signature.mid$signature = paste0(signature.name, ".dep.mid")
    head(ISG.share.signature.mid)
    dim(ISG.share.signature.mid)
    # 442   2

    ISG.share.signature = rbind(ISG.share.signature.high, ISG.share.signature.mid, ISG.share.signature.low)
    head(ISG.share.signature)
    dim(ISG.share.signature)
    # 554   2
    write.csv(ISG.share.signature, paste0("ISG/dep_", signature.name, "_hml", cutoff.percentile, "_ISG.csv"), quote=TRUE)

    group.high = ISG.share.signature[ISG.share.signature$signature %like% "high",,drop=FALSE]
    group.mid = ISG.share.signature[ISG.share.signature$signature %like% "mid",,drop=FALSE]
    group.low = ISG.share.signature[ISG.share.signature$signature %like% "low",,drop=FALSE]
    # t.test(group1[,1], group2[,1])$p.value
    # 0.001361099
    p.hm = paste("p.hm = ", signif(t.test(group.high$ISG, group.mid$ISG)$p.value, 4))
    p.hl = paste("p.hl = ", signif(t.test(group.high$ISG, group.low$ISG)$p.value, 4))
    p.ml = paste("p.ml = ", signif(t.test(group.mid$ISG, group.low$ISG)$p.value, 4))
    ISG.share.signature$signature = gsub(paste0(signature.name, "\\.dep\\.high"), "1.high", ISG.share.signature$signature)
    ISG.share.signature$signature = gsub(paste0(signature.name, "\\.dep\\.mid"), "2.mid", ISG.share.signature$signature)
    ISG.share.signature$signature = gsub(paste0(signature.name, "\\.dep\\.low"), "3.low", ISG.share.signature$signature)

    p1 <- ggplot(ISG.share.signature, aes(x=signature, y= ISG, fill=signature)) + 
      geom_violin(trim=FALSE, linetype="blank", na.rm=TRUE)+
      geom_boxplot(width=0.05, fill="white", outlier.size=0.1, na.rm=TRUE)+
      annotate("text", x=0.5, y=max(na.omit(ISG.share.signature$ISG))*1.2, parse=FALSE, hjust=0, label = p.hm)+
      annotate("text", x=0.5, y=max(na.omit(ISG.share.signature$ISG))*1.1, parse=FALSE, hjust=0, label = p.hl)+
      annotate("text", x=0.5, y=max(na.omit(ISG.share.signature$ISG))*1.0, parse=FALSE, hjust=0, label = p.ml)+
      labs(title=paste0("ISG of signature [", signature.name, "]"), x=paste0(signature.name, " score"), y = "ISG") +  theme_classic() + rremove("legend")
    # p1
    ggsave(paste0("ISG/dep_", signature.name, "_hml", cutoff.percentile, "_ISG.pdf"), p1, width=3, height=4)

    message("[12/14] Analysis-ISG: done!")
    #############################################
 
    # mRNAsi
    dir.create(file.path(outputDir, signature.name, "mRNAsi"), showWarnings = FALSE)

    mRNAsi.share = mRNAsi[rownames(mRNAsi) %in% rownames(dep.t.signature),,drop=FALSE]
    head(mRNAsi.share)
    dim(mRNAsi.share)
    # 554   1

    mRNAsi.share.signature.high = mRNAsi.share[rownames(mRNAsi.share) %in% rownames(dep.t.signature.high),,drop=FALSE]
    mRNAsi.share.signature.high$signature = paste0(signature.name, ".dep.high")
    head(mRNAsi.share.signature.high)
    dim(mRNAsi.share.signature.high)
    # 56  2
    mRNAsi.share.signature.low = mRNAsi.share[rownames(mRNAsi.share) %in% rownames(dep.t.signature.low),,drop=FALSE]
    mRNAsi.share.signature.low$signature = paste0(signature.name, ".dep.low")
    head(mRNAsi.share.signature.low)
    dim(mRNAsi.share.signature.low)
    # 56  2
    mRNAsi.share.signature.mid = mRNAsi.share[!rownames(mRNAsi.share) %in% c(rownames(dep.t.signature.high), rownames(dep.t.signature.low)),,drop=FALSE]
    mRNAsi.share.signature.mid$signature = paste0(signature.name, ".dep.mid")
    head(mRNAsi.share.signature.mid)
    dim(mRNAsi.share.signature.mid)
    # 442   2

    mRNAsi.share.signature = rbind(mRNAsi.share.signature.high, mRNAsi.share.signature.mid, mRNAsi.share.signature.low)
    head(mRNAsi.share.signature)
    dim(mRNAsi.share.signature)
    # 554   2
    write.csv(mRNAsi.share.signature, paste0("mRNAsi/dep_", signature.name, "_hml", cutoff.percentile, "_mRNAsi.csv"), quote=TRUE)

    group.high = mRNAsi.share.signature[mRNAsi.share.signature$signature %like% "high",,drop=FALSE]
    group.mid = mRNAsi.share.signature[mRNAsi.share.signature$signature %like% "mid",,drop=FALSE]
    group.low = mRNAsi.share.signature[mRNAsi.share.signature$signature %like% "low",,drop=FALSE]
    # t.test(group1[,1], group2[,1])$p.value
    # 0.001361099
    p.hm = paste("p.hm = ", signif(t.test(group.high$mRNAsi, group.mid$mRNAsi)$p.value, 4))
    p.hl = paste("p.hl = ", signif(t.test(group.high$mRNAsi, group.low$mRNAsi)$p.value, 4))
    p.ml = paste("p.ml = ", signif(t.test(group.mid$mRNAsi, group.low$mRNAsi)$p.value, 4))
    mRNAsi.share.signature$signature = gsub(paste0(signature.name, "\\.dep\\.high"), "1.high", mRNAsi.share.signature$signature)
    mRNAsi.share.signature$signature = gsub(paste0(signature.name, "\\.dep\\.mid"), "2.mid", mRNAsi.share.signature$signature)
    mRNAsi.share.signature$signature = gsub(paste0(signature.name, "\\.dep\\.low"), "3.low", mRNAsi.share.signature$signature)

    p1 <- ggplot(mRNAsi.share.signature, aes(x=signature, y= mRNAsi, fill=signature)) + 
      geom_violin(trim=FALSE, linetype="blank", na.rm=TRUE)+
      geom_boxplot(width=0.05, fill="white", outlier.size=0.1, na.rm=TRUE)+
      annotate("text", x=0.5, y=max(na.omit(mRNAsi.share.signature$mRNAsi))*1.2, parse=FALSE, hjust=0, label = p.hm)+
      annotate("text", x=0.5, y=max(na.omit(mRNAsi.share.signature$mRNAsi))*1.1, parse=FALSE, hjust=0, label = p.hl)+
      annotate("text", x=0.5, y=max(na.omit(mRNAsi.share.signature$mRNAsi))*1.0, parse=FALSE, hjust=0, label = p.ml)+
      labs(title=paste0("mRNAsi of signature [", signature.name, "]"), x=paste0(signature.name, " score"), y = "mRNAsi") +  theme_classic() + rremove("legend")
    # p1
    ggsave(paste0("mRNAsi/dep_", signature.name, "_hml", cutoff.percentile, "_mRNAsi.pdf"), p1, width=3, height=4)

    message("[13/14] Analysis-mRNAsi: done!")
    #############################################

    # EMT
    dir.create(file.path(outputDir, signature.name, "EMT"), showWarnings = FALSE)

    EMT.share = EMT[rownames(EMT) %in% rownames(dep.t.signature),,drop=FALSE]
    head(EMT.share)
    dim(EMT.share)
    # 250   3

    EMT.share.signature.high = EMT.share[rownames(EMT.share) %in% rownames(dep.t.signature.high),,drop=FALSE]
    EMT.share.signature.high$signature = paste0(signature.name, ".dep.high")
    head(EMT.share.signature.high)
    dim(EMT.share.signature.high)
    # 56  2
    EMT.share.signature.low = EMT.share[rownames(EMT.share) %in% rownames(dep.t.signature.low),,drop=FALSE]
    EMT.share.signature.low$signature = paste0(signature.name, ".dep.low")
    head(EMT.share.signature.low)
    dim(EMT.share.signature.low)
    # 56  2
    EMT.share.signature.mid = EMT.share[!rownames(EMT.share) %in% c(rownames(dep.t.signature.high), rownames(dep.t.signature.low)),,drop=FALSE]
    EMT.share.signature.mid$signature = paste0(signature.name, ".dep.mid")
    head(EMT.share.signature.mid)
    dim(EMT.share.signature.mid)
    # 442   2

    EMT.share.signature = rbind(EMT.share.signature.high, EMT.share.signature.mid, EMT.share.signature.low)
    head(EMT.share.signature)
    dim(EMT.share.signature)
    # 554   2
    write.csv(EMT.share.signature, paste0("EMT/dep_", signature.name, "_hml", cutoff.percentile, "_EMT.csv"), quote=TRUE)

    group.high = EMT.share.signature[EMT.share.signature$signature %like% "high",,drop=FALSE]
    group.mid = EMT.share.signature[EMT.share.signature$signature %like% "mid",,drop=FALSE]
    group.low = EMT.share.signature[EMT.share.signature$signature %like% "low",,drop=FALSE]
    # t.test(group1[,1], group2[,1])$p.value
    # 0.001361099
    p.hm = paste("p.hm = ", signif(t.test(group.high$EMT, group.mid$EMT)$p.value, 4))
    p.hl = paste("p.hl = ", signif(t.test(group.high$EMT, group.low$EMT)$p.value, 4))
    p.ml = paste("p.ml = ", signif(t.test(group.mid$EMT, group.low$EMT)$p.value, 4))
    EMT.share.signature$signature = gsub(paste0(signature.name, "\\.dep\\.high"), "1.high", EMT.share.signature$signature)
    EMT.share.signature$signature = gsub(paste0(signature.name, "\\.dep\\.mid"), "2.mid", EMT.share.signature$signature)
    EMT.share.signature$signature = gsub(paste0(signature.name, "\\.dep\\.low"), "3.low", EMT.share.signature$signature)

    p1 <- ggplot(EMT.share.signature, aes(x=signature, y= EMT, fill=signature)) + 
      geom_violin(trim=FALSE, linetype="blank", na.rm=TRUE)+
      geom_boxplot(width=0.05, fill="white", outlier.size=0.1, na.rm=TRUE)+
      annotate("text", x=0.5, y=max(na.omit(EMT.share.signature$EMT))*1.5, parse=FALSE, hjust=0, label = p.hm)+
      annotate("text", x=0.5, y=max(na.omit(EMT.share.signature$EMT))*1.3, parse=FALSE, hjust=0, label = p.hl)+
      annotate("text", x=0.5, y=max(na.omit(EMT.share.signature$EMT))*1.1, parse=FALSE, hjust=0, label = p.ml)+
      labs(title=paste0("EMT of signature [", signature.name, "]"), x=paste0(signature.name, " score"), y = "EMT") +  theme_classic() + rremove("legend")
    # p1
    ggsave(paste0("EMT/dep_", signature.name, "_hml", cutoff.percentile, "_EMT.pdf"), p1, width=3, height=4)

    message("[14/14] Analysis-EMT: done!")
    #############################################
    
    message("*** deplink finished successfully!")
    message(paste0("*** Please check results in output folder: ", getwd()))
}

