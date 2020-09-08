##' 'load_libs' loads libraries required by "deplink" 
##'
##'
##' @title load_libs
##' @return 
##' @author Xiao Chen
##' @references 
##' @examples
##' outputDir = "~/"


## Main
# load_libs <- function {
	######################################################
	## Loading data
	### meta ###
	meta <- read.csv(system.file("extdata", "avana_ceres_celllines_2019q4.csv", package = "deplink"), header=TRUE, row.names=1)
	head(meta)
	dim(meta)
	# 1677    8
    message("[01/15] Library-meta loaded...")

	# Convert to TCGA Study Abbreviations
	# https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations
	# DepMapID to TCGA tumor type: C:\workspace\OneDrive - cumc.columbia.edu\Box Sync\CX\20190122_cell_vulnerability\data\cellline\meta\Cell_lines_annotations_20181226.csv
	TCGA.tumor <- read.csv(system.file("extdata", "Cell_lines_annotations_20181226.csv", package = "deplink"), header=TRUE)
	TCGA.tumor = TCGA.tumor[!duplicated(TCGA.tumor$depMapID) & !is.na(TCGA.tumor$depMapID),,drop=FALSE]
	rownames(TCGA.tumor) = TCGA.tumor$depMapID
	TCGA.tumor = TCGA.tumor[,c("tcga_code"),drop=FALSE]
	head(TCGA.tumor)
	dim(TCGA.tumor)
	# 1457  1
    message("[02/15] Library-TCGA.tumor loaded...")

	### Dependency ### -1 as essential, 0 as non-essential
	dep <- read.csv(unz(system.file("extdata", "avana_ceres_gene_effects_2019q4_noNA_3digits.zip", package = "deplink"), "avana_ceres_gene_effects_2019q4_noNA_3digits.csv"), header=TRUE, row.names=1)
	head(dep)
	dim(dep)
	# 18333   689
	dep.t = as.data.frame(t(dep))
	rownames(dep.t) = gsub("\\.", "-", rownames(dep.t))
	head(dep.t)
	dim(dep.t)
	# 689 18333

	dep.t.meta = merge(dep.t, meta, by="row.names", all=FALSE)
	rownames(dep.t.meta) = dep.t.meta[,1]
	dep.t.meta = dep.t.meta[,-1]
	head(dep.t.meta)
	dim(dep.t.meta)
	# 667 18341
	dep.t.meta = merge(dep.t.meta, TCGA.tumor, by="row.names", all=FALSE)
	rownames(dep.t.meta) = dep.t.meta[,1]
	dep.t.meta = dep.t.meta[,-1]
	head(dep.t.meta)
	dim(dep.t.meta)
	# 580 18342
    message("[03/15] Library-dependency loaded...")

	### Mutation ###
	mutation <- read.csv(unz(system.file("extdata", "depmap_19q2_mutation_calls.damage.hotspot_matrix.zip", package = "deplink"), "depmap_19q2_mutation_calls.damage.hotspot_matrix.csv"), header=TRUE, row.names=1)
	# mutation = as.data.frame(mutation[rownames(mutation) %in% rownames(dep.t),,drop=FALSE])
	mutation = mutation[,colSums(mutation)>1]
	head(mutation)
	dim(mutation)
	# 1630 16491
    message("[04/15] Library-mutation loaded...")

	### COSMIC ###
	cosmic <- read.csv(system.file("extdata", "CCLE_COSMIC_nodup.csv", package = "deplink"), header=TRUE, row.names=2)
	head(cosmic)
	dim(cosmic)
	# 922  32
	cosmic.share = cosmic[rownames(cosmic) %in% rownames(dep.t), -c(1,2),drop=FALSE]
	head(cosmic.share)
	dim(cosmic.share)
	# 396  30
    message("[05/15] Library-cosmic loaded...")

	### Expression ### 
	exp.TPM <- read.csv(unz(system.file("extdata", "CCLE_depMap_19q4_TPM_2digits.zip", package = "deplink"), "CCLE_depMap_19q4_TPM_2digits.csv"), header=TRUE, row.names=1)
	head(exp.TPM)
	dim(exp.TPM)
	# 16950   554
	exp.TPM.t = as.data.frame(t(exp.TPM))
	rownames(exp.TPM.t) = gsub("\\.", "-", rownames(exp.TPM.t))
	head(exp.TPM.t)
	dim(exp.TPM.t)
	# 554 16950
    message("[06/15] Library-expression loaded...")

	### Chromatin modification ###
	chromatin <- read.csv(system.file("extdata", "CCLE_GlobalChromatinProfiling_20181130_single.csv", package = "deplink"), header=TRUE, row.names=1)
	# all except "H3K18ac0", "H3K23ub1", "H3K56me1", "H3K4ac1"
	chromatin = chromatin[order(rownames(chromatin)), !colnames(chromatin) %in% c("H3K18ac0", "H3K23ub1", "H3K56me1", "H3K4ac1")]
	head(chromatin)
	dim(chromatin)
	# 897  31
    message("[07/15] Library-chromatin.modification loaded...")

	### Drug sensitivity - GDSC ### -2 as sensitive, 2 as resistant
	drug.meta <- read.csv(system.file("extdata", "GDSC_Screened_Compounds.csv", package = "deplink"), header=TRUE, row.names=1)
	head(drug.meta)
	dim(drug.meta)
	# 267   4
	drug <- read.csv(system.file("extdata", "GDSC_v17.3_fitted_dose_response.zscore_DepMap.csv", package = "deplink"), header=TRUE, row.names=1)
	colnames(drug) = gsub("\\.", ":", colnames(drug))
	head(drug)
	dim(drug)
	# 969 266
    message("[08/15] Library-drug.GDSC loaded...")

	### Drug sensitivity - PRISM ### -2 as sensitive, 2 as resistant
	drug2.meta <- read.csv(system.file("extdata", "PRISM_primary_replicate_collapsed_treatment_info.csv", package = "deplink"), header=TRUE, row.names=1)
	rownames(drug2.meta) = gsub(":", ".", rownames(drug2.meta))
	rownames(drug2.meta) = gsub("-", ".", rownames(drug2.meta))
	head(drug2.meta)
	dim(drug2.meta)
	# 4686   10
	drug2 <- read.csv(unz(system.file("extdata", "PRISM_primary_replicate_collapsed_logfold_change.zip", package = "deplink"), "PRISM_primary_replicate_collapsed_logfold_change.csv"), header=TRUE, row.names=1)
	colnames(drug2) = gsub(":", ".", colnames(drug2))
	head(drug2)
	dim(drug2)
	# 578 4686
    message("[09/15] Library-drug.PRISM loaded...")
    
	TMB <- read.csv(system.file("extdata", "depmap_19q1_mutation_calls.TMB.csv", package = "deplink"), header=TRUE, row.names=1)
	head(TMB)
	dim(TMB)
	# 1601    1
    message("[10/15] Library-TMB loaded...")
	
	CNV <- read.csv(system.file("extdata", "public_19q1_gene_cn_nodup.share.sig.cutoff1.csv", package = "deplink"), header=TRUE, row.names=1)
	colnames(CNV) = gsub("burden", "CNV", colnames(CNV))
	head(CNV)
	dim(CNV)
	# 554   3
    message("[11/15] Library-CNV loaded...")

	MSI <- read.csv(system.file("extdata", "Ghandi2019_MSI.csv", package = "deplink"), header=TRUE, row.names=2)
	colnames(MSI) = gsub("CCLE.hc.msi_del", "MSI", colnames(MSI))
	head(MSI)
	dim(MSI)
	# 1330   11
    message("[12/15] Library-MSI loaded...")

	ISG <- read.csv(system.file("extdata", "avana_ceres_celllines_2019q4_meta_ISG.score.csv", package = "deplink"), header=TRUE, row.names=1)
	colnames(ISG) = gsub("ISG.score", "ISG", colnames(ISG))
	head(ISG)
	dim(ISG)
	# 1677    9
    message("[13/15] Library-ISG loaded...")

	mRNAsi <- read.csv(system.file("extdata", "CCLE_depMap_19q1_TPM_3digits_nodup_mRNAsi.csv", package = "deplink"), header=TRUE, row.names=1)
	colnames(mRNAsi) = gsub("mRNAsi", "mRNAsi", colnames(mRNAsi))
	head(mRNAsi)
	dim(mRNAsi)
	# 1165    1
    message("[14/15] Library-ISG loaded...")

	EMT <- read.csv(system.file("extdata", "CCLE_depMap_19q1_zscore_3digits_nodup.share.EMT.score.csv", package = "deplink"), header=TRUE, row.names=1)
	colnames(EMT) = gsub("EMT.score", "EMT", colnames(EMT))
	head(EMT)
	dim(EMT)
	# 1677    9
    message("[15/15] Library-ISG loaded...")
	
    message("All libraries are loaded successfully!")
# }

