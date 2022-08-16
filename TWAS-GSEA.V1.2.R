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
make_option("--expression_ref", action="store", default=NA, type='character',
	help="File containing predicted expression [optional]"),
make_option("--n_cores", action="store", default=1, type='numeric',
	help="Number of cores for permutation testing [optional]"),
make_option("--covar", action="store", default=c('none'), type='character',
	help="Covariates you would like to include. Specify 'none' if you don't want any. [optional]"),
make_option("--pos", action="store", default=NA, type='character',
	help="File containing positions of each feature [required]"),
make_option("--weights", action="store", default=NA, type='character',
	help="Variable used to weight observations [optional]"),
make_option("--use_alt_id", action="store", default=NA, type='character',
	help="Specify alternative column name to match IDs to gene set or property file instead of entrez IDs [optional]"),
make_option("--cor_window", action="store", default=5e6, type='numeric',
	help="Size of window for correlations between genes [optional]"),
make_option("--min_Ngenes", action="store", default=5, type='numeric',
	help="Minimum number of available genes required in gene set for analysis [optional]"),
make_option("--qqplot", action="store", default=T, type='logical',
	help="Specify as F if you do not want a qqplot [optional]"),
make_option("--probit_P_as_Z", action="store", default=T, type='logical',
	help="Specify as F if you want to used abs(TWAS.Z) as the outcome [optional]"),
make_option("--p_cor_method", action="store", default='fdr', type='character',
	help="Select method for correction of multiple testing. Options are the same as the method option for the p.adjust function [optional]"),
make_option("--directional", action="store", default=F, type='logical',
  help="Specify as T if you want to use TWAS.Z as the outcome (i.e. take into account the direction of TWAS association) [optional]"),
make_option("--outlier_threshold", action="store", default='-3,6', type='character',
	help="Thresholds for truncating Z scores [optional]"),
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
make_option("--max_r2", action="store", default=1, type='numeric',
	help="Specify the R2 threshold between genes for pruning [optional]"),
make_option("--min_r2", action="store", default=0.0001, type='numeric',
	help="Specify the R2 threshold between genes assuming independence [optional]"),
make_option("--linear_p_thresh", action="store", default=NA, type='numeric',
	help="Linear model p-value threshold for mixed model analysis [optional]"),
make_option("--output", action="store", default=NA, type='character',
	help="Output file for results [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

opt$covar<- as.character(unlist(strsplit(opt$covar,',')))
opt$outlier_threshold<- as.numeric(unlist(strsplit(opt$outlier_threshold,',')))

sink(file = paste(opt$output,'.log',sep=''), append = F)
if(is.na(opt$twas_results) == T){
	cat('Either expression_ref or input_CorMat must be specified.\n')
	q()
}

if(is.na(opt$output) == T){
	cat('Either expression_ref or input_CorMat must be specified\n.')
	q()
}

if(is.na(opt$pos) == T){
	cat('pos file for weights must be specified.\n')
	q()
}

if(is.na(opt$output) == T){
	cat('Either expression_ref or input_CorMat must be specified.\n')
	q()
}

if(is.na(opt$expression_ref) == T & is.na(opt$input_CorMat) == T){
	cat('Either expression_ref or input_CorMat must be specified for mixed models.\n')
}

if(opt$probit_P_as_Z == T & opt$directional == T){
  cat('Both probit_P_as_Z and directional equal TRUE. probit_P_as_Z will be set to FALSE.\n')
  opt$probit_P_as_Z<-F
}

if(is.na(opt$gmt_file) == T & is.na(opt$prop_file) == T){
	cat('Either gmt_file or prop_file must be specified.\n')
	q()
}

if(is.na(opt$gmt_file) == F & is.na(opt$prop_file) == F){
	cat('Gene set and gene property analysis must be performed separately.\n')
	q()
}

if(opt$self_contained == F & opt$competitive == F){
	cat('Both competitive and self_contained have been set to false.\n')
}

sink()

suppressMessages(library(data.table))
suppressMessages(library(GWASTools))
sink("/dev/null")
suppressMessages(library(WGCNA, quietly=T))
sink()
suppressMessages(library(Matrix))
suppressMessages(library(VGAM))
suppressMessages(library(speedglm))
suppressMessages(library(biomaRt))
suppressMessages(library(qusage))
suppressMessages(library(gdata))
suppressMessages(library(lme4qtl))
suppressMessages(library(lme4))
suppressMessages(library(matrixcalc))
suppressMessages(library(pbkrtest))
suppressMessages(library(foreach))
suppressMessages(library(doMC))
registerDoMC(opt$n_cores)

# Make bdiag function that retains column names
bdiag_withNames<-function(x,y){
	tmp<-bdiag(x,y)
	colnames(tmp)<-c(dimnames(x)[[1]],dimnames(y)[[1]])
	rownames(tmp)<-c(dimnames(x)[[1]],dimnames(y)[[1]])
	return(tmp)
}

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat(
'#################################################################
# TWAS-based Gene Set Enrichment Analysis
# V1.2 26/09/2018
#################################################################

Options are:\n')
print(opt)
cat('Analysis started at',as.character(start.time),'\n')

# Read in TWAS results, removing duplicates and rows with missing TWAS.P.
sink()
TWAS<-data.frame(fread(opt$twas_results))
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('TWAS results file contains',dim(TWAS)[1],'rows.\n')

# Update FILE column to match pos file
file_list<-strsplit(as.character(TWAS$FILE),'/')
file_tab<-lapply(file_list, function(x) x[(length(x)-1):length(x)])
tmp<-data.frame(do.call(rbind, file_tab))
TWAS$FILE<-do.call(paste, c(tmp[,(dim(tmp)[2]-1):dim(tmp)[2]], sep="/"))
rm(file_list,file_tab,tmp)

# Read in .pos file and update P0 and P1 values (This is abug in FUSION).
pos<-data.frame(fread(opt$pos))
TWAS$P0<-NULL
TWAS$P1<-NULL
TWAS<-merge(TWAS,pos[c('WGT','P0','P1')],by.x='FILE',by.y='WGT')
cat('Positional information is available for',dim(TWAS)[1],'TWAS features.\n')

# Add a .5Mb window to the gene coordinates as this is the window for including SNPs as predictors
TWAS$P0<-TWAS$P0-5e5
TWAS$P0[TWAS$P0 < 0]<-0
TWAS$P1<-TWAS$P1+5e5

# Create gene length variable.
TWAS$GeneLength<-TWAS$P1-TWAS$P0

# Remove duplicate features or features with missing TWAS.P
TWAS<-TWAS[!is.na(TWAS$TWAS.P),]
TWAS<-TWAS[!duplicated(TWAS$FILE),]
cat('TWAS contains',dim(TWAS)[1],'unique features with non-missing TWAS.P values.\n')

if(!is.na(opt$use_alt_id)){
	if(opt$use_alt_id == 'ID'){
		# Create an alternate ID column (just for code simplicity later on)
		TWAS$Alt_ID<-TWAS$ID
	} else {
		# Rename the alternate ID column 'Alt_ID'
		names(TWAS)[grepl(opt$use_alt_id, names(TWAS))]<-'Alt_ID'
	}
}

# Convert TWAS 'FILE' column to match the gene expression column names in expression_ref
TWAS$FILE<-sub(".*/", "", TWAS$FILE)
TWAS$FILE<-sub(".wgt.RDat", "", TWAS$FILE)

if(opt$allow_duplicate_ID == F){
	# Remove features with duplicate IDs, retaining the feature with the best R2
	TWAS<-TWAS[order(TWAS$MODELCV.R2),]
	if(is.na(opt$use_alt_id)){
		TWAS<-TWAS[!duplicated(TWAS$ID),]
	} else {
		TWAS<-TWAS[!duplicated(TWAS$Alt_ID),]
	}
	cat('Duplicate IDs removed, leaving',dim(TWAS)[1],'features.\n')
}

if(opt$probit_P_as_Z == T){
  # Create a normally distributed absolute TWAS.Z value using a probit transformation
  TWAS$ZSCORE<-probitlink(1-TWAS$TWAS.P)
  TWAS$ZSCORE[TWAS$ZSCORE < opt$outlier_threshold[1]] <- opt$outlier_threshold[1]
  TWAS$ZSCORE[TWAS$ZSCORE > opt$outlier_threshold[2]] <- opt$outlier_threshold[2]
} 

if(opt$directional == T){
  # Use TWAS.Z as is
  TWAS$ZSCORE<-TWAS$TWAS.Z
  TWAS$ZSCORE[TWAS$ZSCORE > opt$outlier_threshold[2]] <- opt$outlier_threshold[2]
}

if(opt$probit_P_as_Z == F & opt$directional == F){
  # Remove the sign from TWAS.Z values
  TWAS$ZSCORE<-abs(TWAS$TWAS.Z)
  TWAS$ZSCORE[TWAS$ZSCORE > opt$outlier_threshold[2]] <- opt$outlier_threshold[2]
}

if(is.na(opt$use_alt_id)){
	# Merge TWAS data with reference to retrieve entrez IDs
  biomartCacheClear()
	ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
	Genes<-getBM(attributes=c('external_gene_name','entrezgene_id'), mart = ensembl)

	# Remove genes from ensembl info that have duplicate IDs
	Genes<-Genes[!is.na(Genes$entrezgene_id),]
	Genes<-Genes[!duplicated(Genes$entrezgene_id),]
	Genes<-Genes[!is.na(Genes$external_gene_name),]
	Genes<-Genes[!duplicated(Genes$external_gene_name),]

	# Merge TWAS with ensembl info
	TWAS<-merge(TWAS, Genes, by.x='ID', by.y='external_gene_name')
	cat(dim(TWAS)[1],'features have entrez IDs.\n')
}

if(is.na(opt$gmt_file) == F){
	# Read in gene sets of interest
	gene_sets<-read.gmt(opt$gmt_file)
	names(gene_sets)<-gsub("[[:punct:]]", ".", names(gene_sets))
	
	cat('Gene set file contained', length(gene_sets),'gene sets.\n')

	# Create column for each gene set, indicating whether each gene is a member
	TWAS_GS_Mem<-data.frame(TWAS, foreach(i=1:length(gene_sets), .combine=cbind) %dopar% {
		if(is.na(opt$use_alt_id)){
			temp<-data.frame(TWAS$entrezgene_id %in% as.character(unlist(gene_sets[i])))
		} else {
			temp<-data.frame(TWAS$Alt_ID %in% as.character(unlist(gene_sets[i])))
		}
		names(temp)<-names(gene_sets[i])
		temp
	})

	# Remove gene sets which have fewer than opt$min_Ngenes genes available in the TWAS
	TWAS_GS_Mem_only<-TWAS_GS_Mem[(names(TWAS_GS_Mem) %in% names(gene_sets))]
	TWAS_GS_Mem_only_clean<-names(TWAS_GS_Mem_only)[colSums(TWAS_GS_Mem_only) >= opt$min_Ngenes]
	TWAS_GS_Mem_clean<-cbind(TWAS_GS_Mem[!(names(TWAS_GS_Mem) %in% names(gene_sets))], TWAS_GS_Mem_only[TWAS_GS_Mem_only_clean])

	gene_sets_clean<-names(gene_sets[TWAS_GS_Mem_only_clean])

	cat(length(gene_sets_clean),'gene sets have a sufficient number of genes available in the TWAS.\n')
}

if(is.na(opt$prop_file) == F){
	# Read in gene property file
	sink()
	gene_prop<-data.frame(fread(opt$prop_file))
	sink(file = paste(opt$output,'.log',sep=''), append = T)

	cat('Gene property file contained', dim(gene_prop)[2]-1,'properties.\n')

	# Merge with the TWAS data
	if(is.na(opt$use_alt_id)){
		TWAS_GS_Prop<-merge(TWAS, gene_prop, by.x='entrezgene_id', by.y='ID')
	} else {
		TWAS_GS_Prop<-merge(TWAS, gene_prop, by.x='Alt_ID', by.y='ID')
	}
	
	# Retain gene properties with >= opt$min_Ngenes which have non-zero values.
	TWAS_GS_Prop_only<-TWAS_GS_Prop[(names(TWAS_GS_Prop) %in% names(gene_prop)[-1])]
	TWAS_GS_Prop_only_clean<-names(TWAS_GS_Prop_only)[colSums(abs(TWAS_GS_Prop_only)) >= opt$min_Ngenes]
	TWAS_GS_Prop_clean<-cbind(TWAS_GS_Prop[!(names(TWAS_GS_Prop) %in% names(gene_prop)[-1])], TWAS_GS_Prop_only[TWAS_GS_Prop_only_clean])
	
	for(i in TWAS_GS_Prop_only_clean){
	  TWAS_GS_Prop_clean[[i]]<-as.numeric(scale(TWAS_GS_Prop_clean[[i]]))
	}

	TWAS_GS_Mem_clean<-TWAS_GS_Prop_clean

	cat(dim(TWAS_GS_Mem_clean)[1],'genes will be included in the gene property analysis.\n')

	gene_sets_clean<-TWAS_GS_Prop_only_clean
	cat(length(gene_sets_clean),'gene properties have a sufficient number of genes available with non-zero property in the TWAS.\n')
	
}

#########
# Perform standard linear regression without accounting for correlation between genes.
#########
				
cat('Performing competitive linear model... ')
Linear_Results<-foreach(i=1:length(gene_sets_clean), .combine=rbind) %dopar% {
	tryCatch({
		if(opt$covar != 'none'){
			if(is.na(opt$weights)){
				nest_mod<-speedlm.fit(y=TWAS_GS_Mem_clean$ZSCORE, X=cbind(1,as.matrix(TWAS_GS_Mem_clean[c(gene_sets_clean[i],opt$covar)])))
			} else {
				nest_mod<-speedlm.wfit(y=TWAS_GS_Mem_clean$ZSCORE, X=cbind(1,as.matrix(TWAS_GS_Mem_clean[c(gene_sets_clean[i],opt$covar)])),w=abs(TWAS_GS_Mem_clean[[opt$weights]]))
			}
		} else {
			if(is.na(opt$weights)){
				nest_mod<-speedlm.fit(y=TWAS_GS_Mem_clean$ZSCORE, X=cbind(1,as.matrix(TWAS_GS_Mem_clean[c(gene_sets_clean[i])])))
			} else {
				nest_mod<-speedlm.wfit(y=TWAS_GS_Mem_clean$ZSCORE, X=cbind(1,as.matrix(TWAS_GS_Mem_clean[c(gene_sets_clean[i])])),w=abs(TWAS_GS_Mem_clean[[opt$weights]]))
			}
		}
		
		sum<-summary(nest_mod)
		if(i == floor(length(gene_sets_clean)/100*10)){cat('10% ')}
		if(i == floor(length(gene_sets_clean)/100*20)){cat('20% ')}
		if(i == floor(length(gene_sets_clean)/100*30)){cat('30% ')}
		if(i == floor(length(gene_sets_clean)/100*40)){cat('40% ')}
		if(i == floor(length(gene_sets_clean)/100*50)){cat('50% ')}
		if(i == floor(length(gene_sets_clean)/100*60)){cat('60% ')}
		if(i == floor(length(gene_sets_clean)/100*70)){cat('70% ')}
		if(i == floor(length(gene_sets_clean)/100*80)){cat('80% ')}
		if(i == floor(length(gene_sets_clean)/100*90)){cat('90% ')}
		if(i == floor(length(gene_sets_clean)/100*100)){cat('100% ')}
		
		if(is.na(opt$prop_file)){
			data.frame(	GeneSet=gene_sets_clean[i],
						Est=coef(sum)[2, 1],
						SE=coef(sum)[2, 2],
						T=coef(sum)[2, 3],
						N_Mem_Avail=sum(TWAS_GS_Mem_clean[c(gene_sets_clean[i])]==T),
						N_Mem=length(gene_sets[[which(names(gene_sets) == gene_sets_clean[i])]]),
						P=pt(coef(sum)[2, 3], sum$df, lower=FALSE))
		} else {
			data.frame(	GeneSet=gene_sets_clean[i],
					Est=coef(sum)[2, 1],
					SE=coef(sum)[2, 2],
					T=coef(sum)[2, 3],
					P=pt(coef(sum)[2, 3], sum$df, lower=FALSE))
		}
	
	}, error=function(e) NULL)
}
cat('Done!\n')

# Extract gene sets/properties achieving opt$linear_p_thresh
if(is.na(opt$linear_p_thresh)){
	Linear_Results_temp<-Linear_Results
	Linear_Results_temp$P.CORR<-p.adjust(Linear_Results_temp$P, method=opt$p_cor_method)
	gene_sets_clean_forMLM<-as.character(Linear_Results_temp$GeneSet[Linear_Results_temp$P.CORR <= 0.1])
	cat('Using a', opt$p_cor_method,'corrected p-value threshold of 0.1 to select gene sets/properties for competitive mixed model analysis.\n')
} else {
	gene_sets_clean_forMLM<-as.character(Linear_Results$GeneSet[Linear_Results$P <= opt$linear_p_thresh])
	cat('Using a p-value threshold of', opt$linear_p_thresh,' to select gene sets/properties for competitive mixed model analysis.\n')
}

sink()

# Write out linear associations for all gene sets/properties
Linear_Results$P.CORR<-p.adjust(Linear_Results$P, method=opt$p_cor_method)
Linear_Results<-Linear_Results[order(Linear_Results$P),]
write.table(Linear_Results, paste(opt$output,'.linear.txt',sep=''), col.names=T, row.names=F, quote=F)

if(opt$qqplot == T){
	# Create qqplot
		png(paste(opt$output,'.linear.png',sep=''))
		qqPlot(Linear_Results$P)
		dev.off()
}

if(is.na(opt$gmt_file) == F){
	if(sum(Linear_Results$P.CORR <= 0.05) > 0){
		sink(file = paste(opt$output,'.linear.sig.txt',sep=''), append = F)
		Results_sig<-Linear_Results[Linear_Results$P.CORR <= 0.05,]
		Results_sig<-Results_sig[order(Results_sig$P.CORR),]
		for(i in 1:dim(Results_sig)[1]){
			TWAS_SigSet<-TWAS_GS_Mem_clean[TWAS_GS_Mem_clean[[as.character(Results_sig$GeneSet[i])]],]
			TWAS_SigSet<-TWAS_SigSet[c('FILE','ID','CHR','P0','P1','NSNP','NWGT','MODELCV.R2','TWAS.Z','TWAS.P','ZSCORE')]
			TWAS_SigSet$MODELCV.R2<-round(TWAS_SigSet$MODELCV.R2,3)
			TWAS_SigSet$TWAS.Z<-round(TWAS_SigSet$TWAS.Z,3)
			TWAS_SigSet$TWAS.P<-round(TWAS_SigSet$TWAS.P,3)
			TWAS_SigSet$ZSCORE<-round(TWAS_SigSet$ZSCORE,3)
			TWAS_SigSet<-TWAS_SigSet[order(TWAS_SigSet$CHR,TWAS_SigSet$P0),]
			cat('Set No.',i,': ',as.character(Results_sig$GeneSet[i]),' (P.CORR = ',Results_sig$P.CORR[i],')\n',sep='')
			TWAS_SigSet_header <- rbind(names(TWAS_SigSet) , TWAS_SigSet )
			write.fwf(TWAS_SigSet_header,sep='\t',append=T, colnames=F)
			cat('\n')
		}
		sink()
	}
}

sink(file = paste(opt$output,'.log',sep=''), append = T)
if((length(gene_sets_clean_forMLM) > 0 & opt$competitive == T) | opt$self_contained == T){
	cat('Mixed model competitive analysis will be performed for',length(gene_sets_clean_forMLM),'gene sets/properties.\n')
	
	# Sort the results by location
	TWAS_GS_Mem_clean<-TWAS_GS_Mem_clean[order(TWAS_GS_Mem_clean$CHR,TWAS_GS_Mem_clean$P0,TWAS_GS_Mem_clean$P1),]

	if(is.na(opt$input_CorMat) == T){
		# Read in predicted gene expression values for this set of tissue weights
		sink()
		if(substr(opt$expression_ref,(nchar(opt$expression_ref)+1)-3,nchar(opt$expression_ref)) == '.gz'){
			GeneX_all<-data.frame(fread(cmd=paste0('zcat ',opt$expression_ref)))
		} else {	
			GeneX_all<-data.frame(fread(opt$expression_ref))
		}
		
		sink(file = paste(opt$output,'.log',sep=''), append = T)
		GeneX_all<-GeneX_all[-1:-2]
		GeneX_all<-GeneX_all[,apply(GeneX_all,2,function(x) !(all(x==0) | all(is.na(x))))] 
		
		cat('Gene expression values contain', dim(GeneX_all)[2]-2,'non-zero variance features and',dim(GeneX_all)[1],'individuals.\n')
		
		# Extract genes available in TWAS and correlation matrix
		TWAS_GS_Mem_clean$FILE<-gsub(':','.',TWAS_GS_Mem_clean$FILE)
		TWAS_GS_Mem_clean$FILE<-gsub('-','.',TWAS_GS_Mem_clean$FILE)
		
		names(GeneX_all)<-gsub(':','.',names(GeneX_all))
		names(GeneX_all)<-gsub('-','.',names(GeneX_all))
		
		genes_overlap<-intersect(TWAS_GS_Mem_clean$FILE, names(GeneX_all))
		TWAS_GS_Mem_clean<-TWAS_GS_Mem_clean[(TWAS_GS_Mem_clean$FILE %in% genes_overlap),]
		GeneX_all<-GeneX_all[(names(GeneX_all) %in% genes_overlap)]
		GeneX_all<-GeneX_all[match(TWAS_GS_Mem_clean$FILE, names(GeneX_all))]
		
		cat(dim(TWAS_GS_Mem_clean)[1],'features are available in both TWAS and gene expression data.\n')
		
		##########
		# Create block wise correlation matrix for all genes in TWAS
		##########
		
		# Determine gene blocks.
		TWAS_GS_Mem_clean$Block<-NA
		for(i in 1:dim(TWAS_GS_Mem_clean)[1]){
			if(i == 1){
				TWAS_GS_Mem_clean$Block[i]<-1
			} else {
				if(i > 1 & TWAS_GS_Mem_clean$CHR[i] == TWAS_GS_Mem_clean$CHR[i-1] & TWAS_GS_Mem_clean$P1[i] > (TWAS_GS_Mem_clean$P0[i-1] - opt$cor_window) & TWAS_GS_Mem_clean$P0[i] < (TWAS_GS_Mem_clean$P1[i-1] + opt$cor_window)){
					TWAS_GS_Mem_clean$Block[i]<-TWAS_GS_Mem_clean$Block[i-1]
				}
				if(!(i > 1 & TWAS_GS_Mem_clean$CHR[i] == TWAS_GS_Mem_clean$CHR[i-1] & TWAS_GS_Mem_clean$P1[i] > (TWAS_GS_Mem_clean$P0[i-1] - opt$cor_window) & TWAS_GS_Mem_clean$P0[i] < (TWAS_GS_Mem_clean$P1[i-1] + opt$cor_window))){
					TWAS_GS_Mem_clean$Block[i]<-TWAS_GS_Mem_clean$Block[i-1]+1
				}
			}
		}
		
		cat('The genes could be separated into',length(unique(TWAS_GS_Mem_clean$Block)),'blocks.\n')
		
		# Calculate correlation matrix for each block, remove values for genes more than 5Mbs apart, and make it positive definite
		cat('Creating correlation matrix... ')
		cor_block_all<-foreach(i=unique(TWAS_GS_Mem_clean$Block), .combine=bdiag_withNames) %dopar% {
			if(sum(TWAS_GS_Mem_clean$Block == i) == 1){
				cor_block_2<-Matrix(1, nrow = 1, ncol = 1, sparse = TRUE)
				colnames(cor_block_2)<-TWAS_GS_Mem_clean$FILE[TWAS_GS_Mem_clean$Block == i]
				rownames(cor_block_2)<-TWAS_GS_Mem_clean$FILE[TWAS_GS_Mem_clean$Block == i]
			} else {
				cor_block<-abs(WGCNA::cor(as.matrix(GeneX_all[(names(GeneX_all) %in% TWAS_GS_Mem_clean$FILE[TWAS_GS_Mem_clean$Block == i])]), method='pearson'))
				# Remove genes with a correlation exceeding an r2 of opt$max_r2
				tmp<-cor_block
				tmp[!lower.tri(tmp)] <- 0
				keep <- colnames(cor_block)[!apply(tmp,2,function(x) any(abs(x) > sqrt(opt$max_r2)))]
				cor_block_2 <- cor_block[(colnames(cor_block) %in% keep),(colnames(cor_block) %in% keep)]
				TWAS_GS_Mem_clean_Block<-TWAS_GS_Mem_clean[which(TWAS_GS_Mem_clean$Block == i),]
				TWAS_GS_Mem_clean_Block<-TWAS_GS_Mem_clean_Block[(TWAS_GS_Mem_clean_Block$FILE %in% keep),]
				if(length(cor_block_2) == 1){
					cor_block_2<-Matrix(1, nrow = 1, ncol = 1, sparse = TRUE)
					colnames(cor_block_2)<-TWAS_GS_Mem_clean_Block$FILE
					rownames(cor_block_2)<-TWAS_GS_Mem_clean_Block$FILE
				} else {
					if(is.positive.definite(as.matrix(cor_block)) == F){
					cor_block_2<-nearPD(cor_block_2,corr=T)$mat
					}
					# If genes are on a different chromosome or more than opt$cor_window apart then set the correlation to 0
					sparse_struc<-Matrix(0, nrow = dim(TWAS_GS_Mem_clean_Block)[1], ncol = dim(TWAS_GS_Mem_clean_Block)[1], sparse = TRUE)
					for(j in 1:dim(TWAS_GS_Mem_clean_Block)[1]){
						temp<-(TWAS_GS_Mem_clean_Block$CHR == TWAS_GS_Mem_clean_Block$CHR[j] & TWAS_GS_Mem_clean_Block$P1 > (TWAS_GS_Mem_clean_Block$P0[j] - opt$cor_window) & TWAS_GS_Mem_clean_Block$P0 < (TWAS_GS_Mem_clean_Block$P1[j] + opt$cor_window))
						sparse_struc[,j][temp]<-1
					}
					cor_block_2[(sparse_struc[,] != 1)@x]<-0
					# Change correlations with an r2 less opt$min_r2 to 0
					cor_block_2[abs(cor_block_2) < sqrt(opt$min_r2)]<-0
					is.positive.definite(as.matrix(cor_block_2))
					if(is.positive.definite(as.matrix(cor_block_2)) == F){
						cor_block_2<-nearPD(cor_block_2,corr=T)$mat
					}
					cor_block_2<-Matrix(cor_block_2, sparse=T)
				}
			}
			if(i == floor(length(unique(TWAS_GS_Mem_clean$Block))/100*10)){cat('10% ')}
			if(i == floor(length(unique(TWAS_GS_Mem_clean$Block))/100*20)){cat('20% ')}
			if(i == floor(length(unique(TWAS_GS_Mem_clean$Block))/100*30)){cat('30% ')}
			if(i == floor(length(unique(TWAS_GS_Mem_clean$Block))/100*40)){cat('40% ')}
			if(i == floor(length(unique(TWAS_GS_Mem_clean$Block))/100*50)){cat('50% ')}
			if(i == floor(length(unique(TWAS_GS_Mem_clean$Block))/100*60)){cat('60% ')}
			if(i == floor(length(unique(TWAS_GS_Mem_clean$Block))/100*70)){cat('70% ')}
			if(i == floor(length(unique(TWAS_GS_Mem_clean$Block))/100*80)){cat('80% ')}
			if(i == floor(length(unique(TWAS_GS_Mem_clean$Block))/100*90)){cat('90% ')}
			if(i == floor(length(unique(TWAS_GS_Mem_clean$Block))/100*100)){cat('100% ')}
			
			cor_block_2
		}
		
		cat('Done!\n')
		
		# Calculate the proportion of sparse values
		prop_sparse<-sum(cor_block_all == 0)/(dim(cor_block_all)[1]*dim(cor_block_all)[2])
		
		cat('The correlation matrix of gene expression is ',prop_sparse*100,'% sparse.\n',sep='')
		cat('After pruning',dim(cor_block_all)[1],'features remain.\n')
		
		TWAS_GS_Mem_clean<-TWAS_GS_Mem_clean[(TWAS_GS_Mem_clean$FILE %in% colnames(cor_block_all)),]
		cor_block_all<-cor_block_all[match(TWAS_GS_Mem_clean$FILE, colnames(cor_block_all)),match(TWAS_GS_Mem_clean$FILE, colnames(cor_block_all))]
		
		if(opt$save_CorMat ==T){
			saveRDS(cor_block_all,paste(opt$output,'.CorMat.RDS',sep=''))
		}
	}

	if(is.na(opt$input_CorMat) == F){
		cor_block_all<-readRDS(opt$input_CorMat)
		
		cat('Precomputed correlation matrix contains', dim(cor_block_all)[2],'features.\n')
		
		TWAS_GS_Mem_clean$FILE<-gsub(':','.',TWAS_GS_Mem_clean$FILE)
		TWAS_GS_Mem_clean$FILE<-gsub('-','.',TWAS_GS_Mem_clean$FILE)
		
		genes_overlap<-intersect(TWAS_GS_Mem_clean$FILE, colnames(cor_block_all))
		TWAS_GS_Mem_clean<-TWAS_GS_Mem_clean[(TWAS_GS_Mem_clean$FILE %in% genes_overlap),]
		cor_block_all<-cor_block_all[(colnames(cor_block_all) %in% genes_overlap),(colnames(cor_block_all) %in% genes_overlap)]
		cor_block_all<-cor_block_all[match(TWAS_GS_Mem_clean$FILE, colnames(cor_block_all)),match(TWAS_GS_Mem_clean$FILE, colnames(cor_block_all))]
		
		cat(dim(TWAS_GS_Mem_clean)[1],'features are available in both TWAS and gene expression data.\n')
	}

}

if(length(gene_sets_clean_forMLM) != 0){
	if(opt$competitive == T){
		# Run regression without fixed effects
		cat('Modelling random effects for competitive analysis... ')
		if(is.na(opt$weights)){
			if(opt$covar != 'none'){
				mod <- relmatLmer(as.formula(paste('ZSCORE ~ ', paste(opt$covar,collapse=' + '), '(1|FILE)', sep=' + ')), TWAS_GS_Mem_clean, relmat = list(FILE = cor_block_all))
			} else {
				mod <- relmatLmer(as.formula(paste('ZSCORE ~ ', '(1|FILE)', sep=' + ')), TWAS_GS_Mem_clean, relmat = list(FILE = cor_block_all))
			}
		} else {
			if(opt$covar != 'none'){
				mod <- relmatLmer(as.formula(paste('ZSCORE ~ ', paste(opt$covar,collapse=' + '), '(1|FILE)', sep=' + ')), TWAS_GS_Mem_clean, relmat = list(FILE = cor_block_all), weights=abs(TWAS_GS_Mem_clean[[opt$weights]]))
			} else {
				mod <- relmatLmer(as.formula(paste('ZSCORE ~ ', '(1|FILE)', sep=' + ')), TWAS_GS_Mem_clean, relmat = list(FILE = cor_block_all), weights=abs(TWAS_GS_Mem_clean[[opt$weights]]))
			}
		}
		cat('Done!\n')

		# Refit model with genesets as fixed effect.
		cat('Modelling fixed effects for competitive analysis... ')
		Results_Comp<-foreach(i=1:length(gene_sets_clean_forMLM), .combine=rbind) %dopar% {
		  skip_to_next<-F
		  if(opt$covar != 'none'){
				mod_alt<-mod
				mod_X<-mod@pp$X
				mod_alt@pp <- merPredD(X=cbind(mod@pp$X[,1],TWAS_GS_Mem_clean[,gene_sets_clean_forMLM[i]],mod@pp$X[,2:(length(opt$covar)+1)]), Zt=mod@pp$Zt, Lambdat=mod@pp$Lambdat, Lind=mod@pp$Lind, theta=mod@pp$theta, n=nrow(mod@pp$X))
				tryCatch(mod2<-refit(mod_alt, TWAS_GS_Mem_clean$ZSCORE), error = function(e){skip_to_next <<- TRUE})
			} else {
				mod_alt<-mod
				mod_X<-mod@pp$X
				mod_alt@pp <- merPredD(X=cbind(mod@pp$X,TWAS_GS_Mem_clean[,gene_sets_clean_forMLM[i]]), Zt=mod@pp$Zt, Lambdat=mod@pp$Lambdat, Lind=mod@pp$Lind, theta=mod@pp$theta, n=nrow(mod@pp$X))
				tryCatch(mod2<-refit(mod_alt, TWAS_GS_Mem_clean$ZSCORE), error = function(e){skip_to_next <<- TRUE})
			}
			
		  if(skip_to_next == F){
		    
  			coefs<-data.frame(coef(summary(mod2)))
  			if(i == floor(length(gene_sets_clean_forMLM)/100*10)){cat('10% ')}
  			if(i == floor(length(gene_sets_clean_forMLM)/100*20)){cat('20% ')}
  			if(i == floor(length(gene_sets_clean_forMLM)/100*30)){cat('30% ')}
  			if(i == floor(length(gene_sets_clean_forMLM)/100*40)){cat('40% ')}
  			if(i == floor(length(gene_sets_clean_forMLM)/100*50)){cat('50% ')}
  			if(i == floor(length(gene_sets_clean_forMLM)/100*60)){cat('60% ')}
  			if(i == floor(length(gene_sets_clean_forMLM)/100*70)){cat('70% ')}
  			if(i == floor(length(gene_sets_clean_forMLM)/100*80)){cat('80% ')}
  			if(i == floor(length(gene_sets_clean_forMLM)/100*90)){cat('90% ')}
  			if(i == floor(length(gene_sets_clean_forMLM)/100*100)){cat('100% ')}
  			
  			if(is.na(opt$prop_file)){
  				data.frame(	GeneSet=gene_sets_clean_forMLM[i],
  							Estimate=coefs$Estimate[2],
  							SE=coefs$Std..Error[2],
  							T=coefs$t.value[2],
  							N_Mem_Avail=sum(TWAS_GS_Mem_clean[c(gene_sets_clean_forMLM[i])]==T),
  							N_Mem=length(gene_sets[[which(names(gene_sets) == gene_sets_clean_forMLM[i])]]),
  							P=(1 - pnorm(coefs$t.value[2])),
  							row.names=paste(i))
  			} else {
  				data.frame(	GeneSet=gene_sets_clean_forMLM[i],
  							Estimate=coefs$Estimate[2],
  							SE=coefs$Std..Error[2],
  							T=coefs$t.value[2],
  							P=(1 - pnorm(coefs$t.value[2])),
  							row.names=paste(i))
  			}
		  }
		}
		cat('Done!\n')
	}
}

if(opt$self_contained == T){
	# Self contained analysis. NOTE: This doesn't weight genes or allow for covariates.
	cat('Performing self-contained mixed model... ')
	Results_SelfCont<-foreach(i=1:length(gene_sets_clean), .combine=rbind) %dopar% {
		TWAS_GS_Mem_clean_selfCont<-TWAS_GS_Mem_clean[TWAS_GS_Mem_clean[[gene_sets_clean[i]]] == T,]
		if(dim(TWAS_GS_Mem_clean_selfCont)[1] > 1){
			mod <- relmatLmer(ZSCORE ~ (1|FILE), TWAS_GS_Mem_clean_selfCont, relmat = list(FILE = cor_block_all[(colnames(cor_block_all) %in% TWAS_GS_Mem_clean_selfCont$FILE),(colnames(cor_block_all) %in% TWAS_GS_Mem_clean_selfCont$FILE)]))
			coefs<-data.frame(coef(summary(mod)))
			df.KR<-get_ddf_Lb(mod, fixef(mod))
			if(i == floor(length(gene_sets_clean)/100*10)){cat('10% ')}
			if(i == floor(length(gene_sets_clean)/100*20)){cat('20% ')}
			if(i == floor(length(gene_sets_clean)/100*30)){cat('30% ')}
			if(i == floor(length(gene_sets_clean)/100*40)){cat('40% ')}
			if(i == floor(length(gene_sets_clean)/100*50)){cat('50% ')}
			if(i == floor(length(gene_sets_clean)/100*60)){cat('60% ')}
			if(i == floor(length(gene_sets_clean)/100*70)){cat('70% ')}
			if(i == floor(length(gene_sets_clean)/100*80)){cat('80% ')}
			if(i == floor(length(gene_sets_clean)/100*90)){cat('90% ')}
			if(i == floor(length(gene_sets_clean)/100*100)){cat('100% ')}
			if(is.na(opt$prop_file)){
				data.frame(	GeneSet=gene_sets_clean[i],
							Estimate=coefs$Estimate[1],
							SE=coefs$Std..Error[1],
							T=coefs$t.value[1],
							N_Mem_Avail=sum(TWAS_GS_Mem_clean[c(gene_sets_clean_forMLM[i])]==T),
							N_Mem=length(gene_sets[[which(names(gene_sets) == gene_sets_clean_forMLM[i])]]),
							P=(1 - pt(coefs$t.value[1], df.KR,lower=T)),
							row.names=paste(i))
			} else {
				data.frame(	GeneSet=gene_sets_clean[i],
							Estimate=coefs$Estimate[1],
							SE=coefs$Std..Error[1],
							T=coefs$t.value[1],
							P=(1 - pt(coefs$t.value[1], df.KR,lower=T)),
							row.names=paste(i))
			}
		}
	}
	cat('Done!\n')
}

sink()

# Write out results for all gene sets/properties
if(opt$competitive == T & length(gene_sets_clean_forMLM) != 0){
	Results_Comp$P.CORR<-p.adjust(Results_Comp$P, method=opt$p_cor_method, n=length(Linear_Results$P))
	Results_Comp<-Results_Comp[order(Results_Comp$P),]
	write.table(Results_Comp, paste(opt$output,'.competitive.txt',sep=''), col.names=T, row.names=F, quote=F)
}
if(opt$self_contained == T){
	Results_SelfCont$P.CORR<-p.adjust(Results_SelfCont$P, method=opt$p_cor_method, n=length(Linear_Results$P))
	Results_SelfCont<-Results_SelfCont[order(Results_SelfCont$P),]
	write.table(Results_SelfCont, paste(opt$output,'.self_contained.txt',sep=''), col.names=T, row.names=F, quote=F)
}

if(opt$qqplot == T){
	if(opt$competitive == T & length(gene_sets_clean_forMLM) != 0){
		png(paste(opt$output,'.competitive.png',sep=''))
		qqPlot(Results_Comp$P)
		dev.off()
	}
	if(opt$self_contained == T){
		png(paste(opt$output,'.self_contained.png',sep=''))
		qqPlot(Results_SelfCont$P)
		dev.off()
	}
}

# Write out gene-level results for significant gene sets.
if(is.na(opt$gmt_file) == F){	
	if(length(gene_sets_clean_forMLM) != 0 & opt$competitive == T){
		if(sum(Results_Comp$P.CORR <= 0.05) > 0){
			sink(file = paste(opt$output,'.competitive.sig.txt',sep=''), append = F)
			Results_sig<-Results_Comp[Results_Comp$P.CORR <= 0.05,]
			Results_sig<-Results_sig[order(Results_sig$P.CORR),]
			for(i in 1:dim(Results_sig)[1]){
				TWAS_SigSet<-TWAS_GS_Mem_clean[TWAS_GS_Mem_clean[[as.character(Results_sig$GeneSet[i])]],]
				TWAS_SigSet<-TWAS_SigSet[c('FILE','ID','CHR','P0','P1','NSNP','NWGT','MODELCV.R2','TWAS.Z','TWAS.P','ZSCORE')]
				TWAS_SigSet$MODELCV.R2<-round(TWAS_SigSet$MODELCV.R2,3)
				TWAS_SigSet$TWAS.Z<-round(TWAS_SigSet$TWAS.Z,3)
				TWAS_SigSet$TWAS.P<-round(TWAS_SigSet$TWAS.P,3)
				TWAS_SigSet$ZSCORE<-round(TWAS_SigSet$ZSCORE,3)
				TWAS_SigSet<-TWAS_SigSet[order(TWAS_SigSet$CHR,TWAS_SigSet$P0),]
				cat('Set No.',i,': ',as.character(Results_sig$GeneSet[i]),' (P.CORR = ',Results_sig$P.CORR[i],')\n',sep='')
				TWAS_SigSet_header <- rbind(names(TWAS_SigSet) , TWAS_SigSet )
				write.fwf(TWAS_SigSet_header,sep='\t',append=T, colnames=F)
				cat('\n')
			}
			sink()
		}
	}
	if(opt$self_contained == T){
		if(sum(Results_SelfCont$P.CORR <= 0.05) > 0){
			sink(file = paste(opt$output,'.self_contained.sig.txt',sep=''), append = F)
			Results_sig<-Results_SelfCont[Results_SelfCont$P.CORR <= 0.05,]
			Results_sig<-Results_sig[rev(order(Results_sig$Estimate)),]
			for(i in 1:dim(Results_sig)[1]){
				TWAS_SigSet<-TWAS_GS_Mem_clean[TWAS_GS_Mem_clean[[as.character(Results_sig$GeneSet[i])]],]
				TWAS_SigSet<-TWAS_SigSet[c('FILE','ID','CHR','P0','P1','NSNP','NWGT','MODELCV.R2','TWAS.Z','TWAS.P','ZSCORE')]
				TWAS_SigSet$MODELCV.R2<-round(TWAS_SigSet$MODELCV.R2,3)
				TWAS_SigSet$TWAS.Z<-round(TWAS_SigSet$TWAS.Z,3)
				TWAS_SigSet$TWAS.P<-round(TWAS_SigSet$TWAS.P,3)
				TWAS_SigSet$ZSCORE<-round(TWAS_SigSet$ZSCORE,3)
				TWAS_SigSet<-TWAS_SigSet[order(TWAS_SigSet$CHR,TWAS_SigSet$P0),]
				cat('Set No.',i,': ',as.character(Results_sig$GeneSet[i]),' (Mean = ',Results_sig$Estimate[i],', P.CORR = ',Results_sig$P.CORR[i],')\n',sep='')
				TWAS_SigSet_header <- rbind(names(TWAS_SigSet) , TWAS_SigSet )
				write.fwf(TWAS_SigSet_header,sep='\t',append=T, colnames=F)
				cat('\n')
			}
			sink()
		}
	}
}

end.time <- Sys.time()
time.taken <- end.time - start.time

sink(file = paste(opt$output,'.log',sep=''), append = TRUE)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),sep=,'\n')
sink()

