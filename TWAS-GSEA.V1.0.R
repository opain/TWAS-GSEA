#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at Cardiff University under the supervision of Richard Anney and Andrew Pocklington.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
make_option("--twas_results", action="store", default=NA, type='character',
	help="File containing TWAS results [required]"),
make_option("--gmt_file", action="store", default=NA, type='character',
	help="File containing genes-sets in gmt format [optional]"),
make_option("--prop_file", action="store", default=NA, type='character',
	help="File containing gene property values [optional]"),
make_option("--use_twas_id", action="store", default=F, type='logical',
	help="Specify as T if gene property file contains TWAS IDs instead of entrez IDs [optional]"),
make_option("--expression_ref", action="store", default=NA, type='character',
	help="File containing predicted expression [optional]"),
make_option("--n_cores", action="store", default=1, type='numeric',
	help="Number of cores for permutation testing [optional]"),
make_option("--cor_window", action="store", default=5e6, type='numeric',
	help="Number of cores for permutation testing [optional]"),
make_option("--covar", action="store", default=c('GeneLength','NSNP','MODELCV.R2'), type='character',
	help="Covariates you would like to include. Specify 'none' if you don't want any. [optional]"),
make_option("--min_Ngenes", action="store", default=5, type='numeric',
	help="Minimum number of available genes required in gene set for analysis [optional]"),
make_option("--qqplot", action="store", default=T, type='logical',
	help="Specify as F if you do not want a qqplot [optional]"),
make_option("--probit_P_as_Z", action="store", default=T, type='logical',
	help="Specify as F if you want to used abs(TWAS.Z) as the outcome [optional]"),
make_option("--outlier_threshold", action="store", default=3, type='numeric',
	help="Threshold for truncating Z scores [optional]"),
make_option("--save_CorMat", action="store", default=F, type='logical',
	help="Specify T if you would like to save the correlation matrix [optional]"),
make_option("--input_CorMat", action="store", default=NA, type='character',
	help="RDS file containing previously made correlation matrix [optional]"),
make_option("--allow_duplicate_ID", action="store", default=F, type='logical',
	help="Specify T if you would like to include multiple copies of the same gene ID [optional]"),
make_option("--self_contained", action="store", default=F, type='logical',
	help="Specify T if you would like to perform self contained analysis [optional]"),
make_option("--competitive", action="store", default=T, type='logical',
	help="Specify F if you do not want to perform competitive analysis [optional]"),
make_option("--max_r2", action="store", default=0.9, type='numeric',
	help="Specify the R2 threshold between genes for pruning [optional]"),
make_option("--min_r2", action="store", default=0.0001, type='numeric',
	help="Specify the R2 threshold between genes assuming independence [optional]"),
make_option("--output", action="store", default=NA, type='character',
	help="Output file for results [required]")
)

# twas_results should be the output of FUSION.assoc_test.R or a file containing a column called 'FILE' containing the name of the feature (i.e. gene), 'ID' containing the ncbi external gene name, and a column called 'TWAS.P'.
# gmt_file should be a standard gmt file which is a tab delimited file with the name of the gene set in the first column, a second column for notes which can be NA, and then entrez IDs within the gene set in each column after wards.
# expression_ref should be a table readable by fread in the data.table package. Columns are genes, rows are individuals, the first two columns are assumed to be ID variables and will be ignored.

opt = parse_args(OptionParser(option_list=option_list))

sink(file = paste(opt$output,'.log',sep=''), append = F)
if(is.na(opt$twas_results) == T){
	cat('Either expression_ref or input_CorMat must be specified\n.')
	q()
}

if(is.na(opt$output) == T){
	cat('Either expression_ref or input_CorMat must be specified\n.')
	q()
}

if(is.na(opt$output) == T){
	cat('Either expression_ref or input_CorMat must be specified\n.')
	q()
}

if(is.na(opt$expression_ref) == T & is.na(opt$input_CorMat) == T){
	cat('Either expression_ref or input_CorMat must be specified\n.')
	q()
}

if(is.na(opt$gmt_file) == T & is.na(opt$prop_file) == T){
	cat('Either gmt_file or prop_file must be specified\n.')
	q()
}

if(is.na(opt$gmt_file) == F & is.na(opt$prop_file) == F){
	cat('Gene set and gene property analysis must be performed separately.\n.')
	q()
}

if(opt$self_contained == F & opt$competitive == F){
	cat('Both competitive and self_contained have been set to false\n.')
	q()
}

sink()

suppressMessages(library(data.table))
suppressMessages(library(GWASTools))
sink("/dev/null")
suppressMessages(library(WGCNA, quietly=T))
sink()
suppressMessages(library(Matrix))
suppressMessages(library(VGAM))
suppressMessages(library(biomaRt))
suppressMessages(library(qusage))
suppressMessages(library(lme4qtl))
suppressMessages(library(lme4))
suppressMessages(library(matrixcalc))
suppressMessages(library(pbkrtest))
suppressMessages(library(foreach))
suppressMessages(library(doMC))
registerDoMC(opt$n_cores)

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat(
'#################################################################
# TWAS-based Gene Set Enrichment Analysis
# V1.0 02/07/2018
# Remember to cite the lme4qtl package!
#################################################################

Options are:\n')
print(opt)
cat('Analysis started at',as.character(start.time),'\n')

# Read in TWAS results for the PFC.
TWAS<-data.frame(fread(opt$twas_results))
TWAS<-TWAS[!is.na(TWAS$TWAS.P),]
TWAS<-TWAS[!duplicated(TWAS$FILE),]

cat('TWAS contains',dim(TWAS)[1],'unique features with non-missing TWAS.P values.\n')

# Convert TWAS FILE column to match the gene expression column names in expression_ref
TWAS$FILE<-sub(".*/", "", TWAS$FILE)
TWAS$FILE<-sub(".wgt.RDat", "", TWAS$FILE)

# Remove features with duplicate IDs, retaining the feature with the best R2.
if(opt$allow_duplicate_ID == F){
	TWAS<-TWAS[order(TWAS$MODELCV.R2),]
	TWAS<-TWAS[!duplicated(TWAS$ID),]
	cat('Duplicate IDs removed, leaving',dim(TWAS)[1],'unique features.\n')
}

if(opt$probit_P_as_Z == T){
	# Create a normally distributed absolute TWAS.Z value using a probit transformation
	TWAS$TWAS.Z.NORM<-probit(1-TWAS$TWAS.P)
	TWAS$TWAS.Z.NORM[TWAS$TWAS.Z.NORM > opt$outlier_threshold] <- opt$outlier_threshold
	TWAS$TWAS.Z.NORM[TWAS$TWAS.Z.NORM < -opt$outlier_threshold] <- -opt$outlier_threshold
} else {
	# Remove the sign from TWAS.Z values.
	TWAS$TWAS.Z.NORM<-abs(TWAS$TWAS.Z)
	TWAS$TWAS.Z.NORM[TWAS$TWAS.Z.NORM > opt$outlier_threshold] <- opt$outlier_threshold
}

# Convert the TWAS gene names into the human Entrez IDs that are in the gene sets.
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
Genes<-getBM(attributes=c('external_gene_name','entrezgene','start_position','end_position'), mart = ensembl)

if(opt$use_twas_id == F){
	Genes<-Genes[!is.na(Genes$entrezgene),]
	Genes<-Genes[!duplicated(Genes$entrezgene),]
}

Genes<-Genes[!is.na(Genes$external_gene_name),]
Genes<-Genes[!duplicated(Genes$external_gene_name),]

TWAS_withEntrez<-merge(TWAS, Genes, by.x='ID', by.y='external_gene_name')

# Add a .5Mb window to the gene coordinates as this is the window for including SNPs as predictors.
TWAS_withEntrez$start_position<-TWAS_withEntrez$start_position-5e5
TWAS_withEntrez$start_position[TWAS_withEntrez$start_position < 0]<-0
TWAS_withEntrez$end_position<-TWAS_withEntrez$end_position+5e5

if(opt$use_twas_id == F){
	cat(dim(TWAS_withEntrez)[1],'genes have entrez IDs.\n')
}

TWAS_GS<-TWAS_withEntrez[c('FILE','entrezgene','ID','TWAS.Z.NORM','MODELCV.R2','CHR','NSNP','start_position','end_position')]
TWAS_GS$GeneLength<-TWAS_GS$end_position-TWAS_GS$start_position
TWAS_GS<-TWAS_GS[!is.na(TWAS_GS$entrezgene),]

if(is.na(opt$gmt_file) == F){
	# Read in gene sets of interest.
	gene_sets<-read.gmt(opt$gmt_file)

	cat('Gene set file contained', length(gene_sets),'gene sets.\n')

	# Create column for each gene set, indicating whether each gene is a member.
	TWAS_GS_Mem<-data.frame(TWAS_GS, foreach(i=1:length(gene_sets), .combine=cbind) %dopar% {
		temp<-data.frame(TWAS_GS$entrezgene %in% as.character(unlist(gene_sets[i])))
		names(temp)<-names(gene_sets[i])
		temp
	})

	# Remove gene sets which have fewer than opt$min_Ngenes genes available in the TWAS
	TWAS_GS_Mem_only<-TWAS_GS_Mem[-1:-10]
	TWAS_GS_Mem_only_min5<-names(TWAS_GS_Mem_only)[colSums(TWAS_GS_Mem_only) >= opt$min_Ngenes]
	TWAS_GS_Mem_min5<-cbind(TWAS_GS_Mem[1:10], TWAS_GS_Mem_only[TWAS_GS_Mem_only_min5])

	gene_sets_min5<-names(gene_sets[TWAS_GS_Mem_only_min5])

	cat(length(gene_sets_min5),'gene sets have a sufficient number of genes available in the TWAS.\n')
}

if(is.na(opt$prop_file) == F){
	# Read in gene property file
	gene_prop<-data.frame(fread(opt$prop_file))

	cat('Gene property file contained', dim(gene_prop)[2]-1,'properties.\n')

	# Merge with the TWAS information.
	if(opt$use_twas_id == F){
		TWAS_GS_Prop<-merge(TWAS_GS, gene_prop, by.x='entrezgene', by.y='ID')
	} else {
		TWAS_GS_Prop<-merge(TWAS_GS, gene_prop, by.x='ID', by.y='ID')
	}

	TWAS_GS_Mem_min5<-TWAS_GS_Prop

	cat(dim(TWAS_GS_Mem_min5)[1],'genes will be included in the gene property analysis.\n')

	gene_sets_min5<-names(gene_prop[-1])
}

# Sort the results by location
TWAS_GS_Mem_min5<-TWAS_GS_Mem_min5[order(TWAS_GS_Mem_min5$CHR,TWAS_GS_Mem_min5$start_position,TWAS_GS_Mem_min5$end_position),]

if(is.na(opt$input_CorMat) == T){
	# Read in predicted gene expression values for this set of tissue weights
	GeneX_all<-data.frame(fread(opt$expression_ref))
	GeneX_all<-GeneX_all[-1:-2]
	GeneX_all<-GeneX_all[,apply(GeneX_all,2,function(x) !all(x==0))] 

	cat('Gene expression values contain', dim(GeneX_all)[2]-2,'non-zero variance features and',dim(GeneX_all)[1],'individuals.\n')

	# Extract genes available in TWAS and correlation matrix
	genes_overlap<-intersect(TWAS_GS_Mem_min5$FILE, names(GeneX_all))
	TWAS_GS_Mem_min5<-TWAS_GS_Mem_min5[(TWAS_GS_Mem_min5$FILE %in% genes_overlap),]
	GeneX_all<-GeneX_all[(names(GeneX_all) %in% genes_overlap)]
	GeneX_all<-GeneX_all[match(TWAS_GS_Mem_min5$FILE, names(GeneX_all))]

	cat(dim(TWAS_GS_Mem_min5)[1],'features are available in both TWAS and gene expression data.\n')

	##########
	# Create block wise correlation matrix for all genes in TWAS
	##########

	# Determine gene blocks.
	TWAS_GS_Mem_min5$Block<-NA
	for(i in 1:dim(TWAS_GS_Mem_min5)[1]){
		if(i == 1){
			TWAS_GS_Mem_min5$Block[i]<-1
		} else {
			if(i > 1 & TWAS_GS_Mem_min5$CHR[i] == TWAS_GS_Mem_min5$CHR[i-1] & TWAS_GS_Mem_min5$end_position[i] > (TWAS_GS_Mem_min5$start_position[i-1] - opt$cor_window) & TWAS_GS_Mem_min5$start_position[i] < (TWAS_GS_Mem_min5$end_position[i-1] + opt$cor_window)){
				TWAS_GS_Mem_min5$Block[i]<-TWAS_GS_Mem_min5$Block[i-1]
			}
			if(!(i > 1 & TWAS_GS_Mem_min5$CHR[i] == TWAS_GS_Mem_min5$CHR[i-1] & TWAS_GS_Mem_min5$end_position[i] > (TWAS_GS_Mem_min5$start_position[i-1] - opt$cor_window) & TWAS_GS_Mem_min5$start_position[i] < (TWAS_GS_Mem_min5$end_position[i-1] + opt$cor_window))){
				TWAS_GS_Mem_min5$Block[i]<-TWAS_GS_Mem_min5$Block[i-1]+1
			}
		}
	}

	cat('The genes could be separated into',length(unique(TWAS_GS_Mem_min5$Block)),'blocks.\n')

	# Calculate correlation matrix for each block, remove values for genes more than 5Mbs apart, and make it positive definite
	cor_block_all<-Matrix(0, nrow = 0, ncol = 0, sparse = TRUE)
	TWAS_GS_Mem_min5_Block_all<-NULL
	for(i in unique(TWAS_GS_Mem_min5$Block)){
		if(sum(TWAS_GS_Mem_min5$Block == i) == 1){
			single_value<-Matrix(1, nrow = 1, ncol = 1, sparse = TRUE)
			cor_block_all<-bdiag(cor_block_all,single_value)
			TWAS_GS_Mem_min5_Block<-TWAS_GS_Mem_min5[TWAS_GS_Mem_min5$Block == i,]
		} else {
			cor_block<-WGCNA::cor(as.matrix(GeneX_all[(names(GeneX_all) %in% TWAS_GS_Mem_min5$FILE[TWAS_GS_Mem_min5$Block == i])]), method='pearson')
			# Remove genes that have a correlation above 0.98 and the same ID
			tmp<-cor_block
			tmp[!lower.tri(tmp)] <- 0
			keep <- colnames(cor_block)[!apply(tmp,2,function(x) any(x > sqrt(opt$max_r2)))]
			cor_block_2 <- cor_block[(colnames(cor_block) %in% keep),(colnames(cor_block) %in% keep)]
			
			TWAS_GS_Mem_min5_Block<-TWAS_GS_Mem_min5[which(TWAS_GS_Mem_min5$Block == i),]
			TWAS_GS_Mem_min5_Block<-TWAS_GS_Mem_min5_Block[(TWAS_GS_Mem_min5_Block$FILE %in% keep),]
			
			if(length(cor_block_2) == 1){
				single_value<-Matrix(1, nrow = 1, ncol = 1, sparse = TRUE)
				cor_block_all<-bdiag(cor_block_all,single_value)
			} else {
				cor_block_2<-cor_block_2[match(TWAS_GS_Mem_min5_Block$FILE, colnames(cor_block_2)),match(TWAS_GS_Mem_min5_Block$FILE, rownames(cor_block_2))]
			
				if(is.positive.definite(as.matrix(cor_block)) == F){
				cor_block_2<-nearPD(cor_block_2,corr=T)$mat
				}
				sparse_struc<-Matrix(0, nrow = dim(TWAS_GS_Mem_min5_Block)[1], ncol = dim(TWAS_GS_Mem_min5_Block)[1], sparse = TRUE)
				for(j in 1:dim(TWAS_GS_Mem_min5_Block)[1]){
					temp<-(TWAS_GS_Mem_min5_Block$CHR == TWAS_GS_Mem_min5_Block$CHR[j] & TWAS_GS_Mem_min5_Block$end_position > (TWAS_GS_Mem_min5_Block$start_position[j] - opt$cor_window) & TWAS_GS_Mem_min5_Block$start_position < (TWAS_GS_Mem_min5_Block$end_position[j] + opt$cor_window))
					sparse_struc[,j][temp]<-1
				}
				cor_block_2[(sparse_struc[,] != 1)@x]<-0
				cor_block_2[abs(cor_block_2) < sqrt(opt$min_r2)]<-0
				is.positive.definite(as.matrix(cor_block_2))
				if(is.positive.definite(as.matrix(cor_block_2)) == F){
					cor_block_2<-nearPD(cor_block_2,corr=T)$mat
				}
				cor_block_2<-Matrix(cor_block_2, sparse=T)
				cor_block_all<-bdiag(cor_block_all,cor_block_2)
			}
		}
		TWAS_GS_Mem_min5_Block_all<-rbind(TWAS_GS_Mem_min5_Block_all, TWAS_GS_Mem_min5_Block)
		cat('Block',i,'complete.\n')
	}

	# Set dimnames
	cor_block_all@Dimnames<-list(TWAS_GS_Mem_min5_Block_all$FILE,TWAS_GS_Mem_min5_Block_all$FILE)

	# Calculate the proportion of sparse values
	prop_sparse<-sum(cor_block_all == 0)/(dim(cor_block_all)[1]*dim(cor_block_all)[2])

	cat('The correlation matrix of gene expression is ',prop_sparse*100,'% sparse.\n',sep='')
	cat('After pruning',dim(TWAS_GS_Mem_min5_Block_all)[1],'features remain.\n')

	TWAS_GS_Mem_min5<-TWAS_GS_Mem_min5_Block_all

	if(opt$save_CorMat ==T){
		saveRDS(cor_block_all,paste(opt$output,'.CorMat.RDS',sep=''))
	}
}

if(is.na(opt$input_CorMat) == F){
	cor_block_all<-readRDS(opt$input_CorMat)

	cat('Precomputed correlation matrix contains', dim(cor_block_all)[2],'features.\n')

	genes_overlap<-intersect(TWAS_GS_Mem_min5$FILE, colnames(cor_block_all))
	TWAS_GS_Mem_min5<-TWAS_GS_Mem_min5[(TWAS_GS_Mem_min5$FILE %in% genes_overlap),]
	cor_block_all<-cor_block_all[(colnames(cor_block_all) %in% genes_overlap),(colnames(cor_block_all) %in% genes_overlap)]
	cor_block_all<-cor_block_all[match(TWAS_GS_Mem_min5$FILE, colnames(cor_block_all)),match(TWAS_GS_Mem_min5$FILE, colnames(cor_block_all))]

	cat(dim(TWAS_GS_Mem_min5)[1],'features are available in both TWAS and gene expression data.\n')
}

sink()

if(opt$competitive == T){
	if(opt$covar != 'none'){
		# Run regression for each gene set
		pb <- txtProgressBar(min = 1, max = length(gene_sets_min5), style = 3)
		Results<-foreach(i=1:length(gene_sets_min5), .combine=rbind) %dopar% {
			mod <- relmatLmer(as.formula(paste('TWAS.Z.NORM ~ TWAS_GS_Mem_min5[,gene_sets_min5[i]]', paste(opt$covar,collapse=' + '), '(1|FILE)', sep=' + ')), TWAS_GS_Mem_min5, relmat = list(FILE = cor_block_all))
			coefs<-data.frame(coef(summary(mod)))
			setTxtProgressBar(pb, i)
			if(is.na(opt$prop_file) == T){
			data.frame(	GeneSet=gene_sets_min5[i],
						NGenesAvail=paste(sum(TWAS_GS_Mem_min5[,names(TWAS_GS_Mem_min5) == gene_sets_min5[i]]), '/', length(gene_sets_min5[[i]]), sep=''),
						Estimate=coefs$Estimate[2],
						SE=coefs$Std..Error[2],
						T=coefs$t.value[2],
						P=(1 - pnorm(coefs$t.value[2])),
						row.names=paste(i))
			} else {
				data.frame(	GeneSet=gene_sets_min5[i],
						Estimate=coefs$Estimate[2],
						SE=coefs$Std..Error[2],
						T=coefs$t.value[2],
						P=(1 - pnorm(coefs$t.value[2])),
						row.names=paste(i))
			}
		}
	} else {
		# Run regression for each gene set
		pb <- txtProgressBar(min = 1, max = length(gene_sets_min5), style = 3)
		Results<-foreach(i=1:length(gene_sets_min5), .combine=rbind) %dopar% {
			mod <- relmatLmer(as.formula(paste('TWAS.Z.NORM ~ TWAS_GS_Mem_min5[gene_sets_min5[i]]', '(1|FILE)', sep=' + ')), TWAS_GS_Mem_min5, relmat = list(FILE = cor_block_all))
			coefs<-data.frame(coef(summary(mod)))
			setTxtProgressBar(pb, i)
			if(is.na(opt$prop_file) == T){
				data.frame(	GeneSet=gene_sets_min5[i],
							NGenesAvail=paste(sum(TWAS_GS_Mem_min5[,names(TWAS_GS_Mem_min5) == gene_sets_min5[i]]), '/', length(gene_sets_min5[[i]]), sep=''),
							Estimate=coefs$Estimate[2],
							SE=coefs$Std..Error[2],
							T=coefs$t.value[2],
							P=(1 - pnorm(coefs$t.value[2])),
							row.names=paste(i))
			} else {
				data.frame(	GeneSet=gene_sets_min5[i],
							Estimate=coefs$Estimate[2],
							SE=coefs$Std..Error[2],
							T=coefs$t.value[2],
							P=(1 - pnorm(coefs$t.value[2])),
							row.names=paste(i))
			}
		}
	}
}

if(opt$self_contained == T){
	# Self contained analysis
	pb <- txtProgressBar(min = 1, max = length(gene_sets_min5), style = 3)
	Results_SelfCont<-foreach(i=1:length(gene_sets_min5), .combine=rbind) %dopar% {
		TWAS_GS_Mem_min5_selfCont<-TWAS_GS_Mem_min5[TWAS_GS_Mem_min5[[gene_sets_min5[i]]] == T,] 
		mod <- relmatLmer(TWAS.Z.NORM ~ (1|FILE), TWAS_GS_Mem_min5_selfCont, relmat = list(FILE = cor_block_all[(colnames(cor_block_all) %in% TWAS_GS_Mem_min5_selfCont$FILE),(colnames(cor_block_all) %in% TWAS_GS_Mem_min5_selfCont$FILE)]))
			coefs<-data.frame(coef(summary(mod)))
			df.KR<-get_ddf_Lb(mod, fixef(mod))
			setTxtProgressBar(pb, i)
		if(is.na(opt$prop_file) == T){
			data.frame(	GeneSet=gene_sets_min5[i],
						NGenesAvail=paste(sum(TWAS_GS_Mem_min5[,names(TWAS_GS_Mem_min5) == gene_sets_min5[i]]), '/', length(gene_sets_min5[[i]]), sep=''),
						Estimate=coefs$Estimate[1],
						SE=coefs$Std..Error[1],
						T=coefs$t.value[1],
						P=(1 - pt(coefs$t.value[1], df.KR,lower=T)),
						row.names=paste(i))
		} else {
			data.frame(	GeneSet=gene_sets_min5[i],
					Estimate=coefs$Estimate[1],
					SE=coefs$Std..Error[1],
					T=coefs$t.value[1],
					P=(1 - pt(coefs$t.value[1], df.KR,lower=T)),
					row.names=paste(i))
		}
	}
}

# Write out results
if(opt$competitive == T){
	write.table(Results, paste(opt$output,'.competitive.txt',sep=''), col.names=T, row.names=F, quote=F)
}
if(opt$self_contained == T){
	write.table(Results_SelfCont, paste(opt$output,'.self_contained.txt',sep=''), col.names=T, row.names=F, quote=F)
}

if(opt$qqplot == T){
	# Create qqplot
	if(opt$competitive == T){
		png(paste(opt$output,'.competitive.png',sep=''))
		qqPlot(Results$P)
		dev.off()
	}
	if(opt$self_contained == T){
		png(paste(opt$output,'.self_contained.png',sep=''))
		qqPlot(Results_SelfCont$P)
		dev.off()
	}
}

end.time <- Sys.time()
time.taken <- end.time - start.time

sink(file = paste(opt$output,'.log',sep=''), append = TRUE)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),sep=,'\n')
sink()

