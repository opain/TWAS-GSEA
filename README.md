# TWAS-based Gene Set Enrichment Analysis (TWAS-GSEA)

TWAS-GSEA is a tool for performing gene set or gene property analysis based TWAS results. It uses a similar method to MAGMA, in that it uses a mixed model to test for enrichment, specifying a gene-gene correlation matrix as a random effect to avoid bias due to non-independent observations. It uses the fantastic package lme4qtl, which enables the use of sparse matrices in mixed models, making this analysis computationally feasible.

TWAS-GSEA was written to analyse the output of the FUSION's [**FUSION.assoc_test.R** ](https://github.com/gusevlab/fusion_twas/blob/master/FUSION.assoc_test.R) script, though it could be used analyse any gene level association results. 

## Getting started

### Prerequisites

* Install R and the following packages:
  
```R
# Install packages from the CRAN
install.packages(c('data.table','optparse','WGCNA','Matrix','VGAM','gdata','lme4qtl','lme4','matrixcalc','pbkrtest','foreach','doMC'))

# Install packages from Bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite(c('GWASTools','biomaRt','qusage'))

# Install packages from GitHub
library(devtools)
install_github("variani/lme4qtl")

```

* Perform TWAS using FUSION:
  * Instructions on how to perform a TWAS are available [here](http://gusevlab.org/projects/fusion/).

* Impute gene expression levels in a reference sample:
  * Instructions on how to impute gene expression levels are [here](https://github.com/opain/Predicting-TWAS-features).



### Input files

##### --twas_results

The output of [**FUSION.assoc_test.R** ](https://github.com/gusevlab/fusion_twas/blob/master/FUSION.assoc_test.R) or a file containing the following columns FILE, ID, P0, P1, TWAS.Z, TWAS.P.  Per chromosome files should be combined into a single file. An example is available [here](http://gitlab.psycm.cf.ac.uk/mpmop/gene-expression-risk-scoring/blob/master/ukbiobank-2017-1160-prePRS-fusion.tsv.GW). Gene IDs are expected to be gene symbols (this can be changed using --use_alt_id parameter).

##### --pos

A file containing the start and stop coordinates of each feature. This should be the .pos file used to perform the TWAS. Gene IDs are expected to be gene symbols (this can be changed using --use_alt_id parameter).

##### --expression_ref

A file containing feature predictions in the target sample. This is output of the FeaturePred script. The first two columns are FID and IID, then each column contains feature predictions for each individual. An example is available here. The gene expression column names must match the values in the FILE column in the --twas_results file. IFRisk ignores the substring before the last '/' and the '.wgt.RDat' string when matching. For example, the column name for the gene expression corresponding to the first value of the example TWAS results should be 'CMC.LOC643837'. The file can whitespace or comma delimited. If the file name ends .gz, the file will be assumed to gzipped.

##### --gmt_file (for gene set analysis)

A standard .gmt file which contains gene set names in the first column, a second column which can be ignored by the analysis, and then a series of entrez ids. This file must be tab delimited. An example can be found here.

##### --prop_file (for gene property analysis)

The first column should have the header 'ID' and contain gene ids. These are assumed to be entrez IDs. Each column after wards should contain values for each gene, with a header stating the gene property. 

### Optional parameters

##### --n_cores

Number of cores for parallel computing.

Default = 1

##### --covar

Covariates you would like to include. The covariate data must be in the TWAS file, except gene length. Covariate names must be comma separated (e.g. GeneLength,NSNP,MODELCV.R2). Specify 'none' if you don't want include any covariates.

Default = 'none'

##### --weights

Variable used to weight observations. The variable must be in the TWAS file.

Default = NA

##### --use_alt_id

Specify column name in TWAS results file to be used for matching with the gmt or property file.

Default = NA

##### --cor_window

Size of window for correlations between genes.

Default = 5e6

##### --min_Ngenes

Minimum number of available genes required in gene set for analysis.

Default = 5

##### --qqplot

Specify as F if you do not want a qqplot.

Default = T

##### --probit_P_as_Z

Specify as F if you want to used abs(TWAS.Z) as the outcome.

Default = T

##### --p_cor_method

Select method for correction of multiple testing. Options are the same as the method option for the p.adjust function.

##### --outlier_threshold

Threshold for truncating Z scores.

Default = 3

##### --save_CorMat

Specify T if you would like to save the correlation matrix.

Default = F

##### --input_CorMat

RDS file containing previously made correlation matrix.

Default = F

##### --allow_duplicate_ID

Specify T if you would like to include multiple copies of the same gene ID. Otherwise only the version of the ID with the best R-squared will be retained.

Default = F

##### --self_contained

Specify T if you would like to perform self contained analysis.

Default = F

##### --competitive

Specify F if you do not want to perform competitive analysis.

Default = T

##### --max_r2

Specify the R-squared threshold between genes for pruning.

Default = 1

##### --min_r2

Specify the R-squared threshold between genes assuming independence.

Default = 0.0001

##### --linear_p_thresh

Linear model p-value threshold for mixed model analysis. Default behavior is to use a multiple testing corrected p-value threshold of 0.1.

Default = NA

##### --output

Output file prefix for results files. This must be specified.

Default = NULL



### Output files

##### '.txt'

These space delimited files will contain all results for either competitive linear models, competitive mixed models, or self-contained mixed models.

##### '.sig.txt'

These files will contain a breakdown of the genes within gene sets achieving significance for either competitive linear models, competitive mixed models, or self-contained mixed models.

##### '.png'

A QQ-plot, comparing the observed association results to a null distribution.

##### '.log'

This is a log file containing general information on the time taken, any errors, the number of genes at different stages and more.



## Example

##### When using default settings:

```sh
Rscript TWAS-GSEA.V1.2.R \
	--twas_results ukbiobank-2017-1160-prePRS-fusion.tsv.GW \
	--pos CMC.BRAIN.RNASEQ.pos \
	--gmt_file GO_STRICT_PC_10-2000_MAGMA.gmt \
	--expression_ref CMC.BRAIN.RNASEQ_GeneX_all_MINI.csv \
	--output demo
```

## Help
If you have any questions or comments use the [google group](https://groups.google.com/forum/#!forum/twas-related-r-scripts), or email oliver.pain@kcl.ac.uk.


