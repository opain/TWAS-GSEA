# TWAS-based Gene Set Enrichment Analysis (TWAS-GSEA)

TWAS-GSEA is a tool for performing gene set or gene property analysis based TWAS results. It uses a similar method to MAGMA, in that it uses a mixed model to test for enrichment, specifying a gene-gene correlation matrix as a random effect to avoid bias due to non-independent observations. It uses the fantastic package lme4qtl, which the use of sparse matrices in mixed models, making this analysis computationally feasible.

TWAS-GSEA was written to analyse the output of the FUSION's [**FUSION.assoc_test.R** ](https://github.com/gusevlab/fusion_twas/blob/master/FUSION.assoc_test.R) script, though it could be used analyse any gene level association results. 



## Getting started

### Prerequisites

* Install the following R packages:
  * data.table
  * optparse
  * GWASTools
  * WGCNA
  * Matrix
  * VGAM
  * biomaRt
  * qusage
  * lme4qtl
  * lme4
  * matrixcalc
  * pbkrtest
  * foreach
  * doMC

* Perform TWAS using FUSION:
  * Instructions on how to perform a TWAS are available [here](http://gusevlab.org/projects/fusion/).
* Impute gene expression levels in the target sample:
  * Instructions on how to impute gene expression levels are [here](http://gusevlab.org/projects/fusion/).



### Input files

##### --twas_results

The output of [**FUSION.assoc_test.R** ](https://github.com/gusevlab/fusion_twas/blob/master/FUSION.assoc_test.R) or a file containing the following columns FILE, P0, P1, TWAS.Z, TWAS.P.  Per chromosome files should be combined into a single file. An example is available [here](http://gitlab.psycm.cf.ac.uk/mpmop/gene-expression-risk-scoring/blob/master/ukbiobank-2017-1160-prePRS-fusion.tsv.GW). 

##### --expression_ref

A file containing gene expression values in a reference sample (e.g. the FUSION 1000 genomes reference).  The first two columns are should FID and IID, then each column should contain gene expression data. An example is available [here](http://gitlab.psycm.cf.ac.uk/mpmop/gene-expression-risk-scoring/blob/master/CMC.BRAIN.RNASEQ_GeneX_all_MINI.csv). The gene expression column names must match the values in the FILE column in the --twas_results file. IFRisk ignores the substring before the last '/' and the '.wgt.RDat' string when matching. For example, the column name for the gene expression corresponding to the first value of the [example TWAS results](http://gitlab.psycm.cf.ac.uk/mpmop/gene-expression-risk-scoring/blob/master/ukbiobank-2017-1160-prePRS-fusion.tsv.GW) should be 'CMC.LOC643837'. The file must be readable by the fread function in R.

##### --gmt_file (for gene set analysis)

A standard .gmt file which contains gene set names in the first column, a second column which can be ignored by the analysis, and then a series of entrez ids. This file must be tab delimited. An example can be found here.

##### --prop_file (for gene property analysis)

The first column should have the header 'ID' and contain gene ids. These are assumed to be entrez IDs. Each column after wards should contain values for each gene, with a header stating the gene property. 



### Output files

##### '.competitive.txt' 

This is a space delimited file containing the results of the competitive gene set or property analysis.

##### '.self_contained.txt'

This is a space delimited file containing the results of the self-contained gene set analysis.

##### '.png'

A QQ-plot, comparing the observed association results to a null distribution.

##### '.log'

This is a log file containing general information on the time taken, any errors, the number of genes at different stages and more.



### Optional parameters

##### --use_twas_id

Specify as T if gene property file contains IDs in the TWAS file ID column instead of entrez IDs.

Default = F

##### --n_cores

Number of cores for permutation testing.

Default = 1

##### --cor_window

Size of window for correlations between genes.

Default = 5e6

##### --covar

Covariates you would like to include. Specify 'none' if you don't want any. The covariate data must be in the TWAS file, except gene length.

Default = GeneLength,NSNP,MODELCV.R2

##### --min_Ngenes

Minimum number of available genes required in gene set for analysis.

Default = 5

##### --qqplot

Specify as F if you do not want a qqplot.

Default = T

##### --probit_P_as_Z

Specify as F if you want to used abs(TWAS.Z) as the outcome.

Default = T

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

Default = 0.9

##### --min_r2

Specify the R-squared threshold between genes assuming independence.

Default = 0.0001

##### --output

Output file prefix for results files. This must be specified.

Default = NULL



## Example

##### When using default settings:

```R
Rscript TWAS-GSEA.V1.0.R \
	--twas_results ukbiobank-2017-1160-prePRS-fusion.tsv.GW \
	--gmt_file GO_STRICT_PC_10-2000_MAGMA.gmt \
	--expression_ref CMC.BRAIN.RNASEQ_GeneX_all_MINI.csv \
	--output demo
```



## Help

This script was written by Dr Oliver Pain under the supervision of Dr Richard Anney whilst at the MRC Centre for Neuropsychiatric Genetics and Genomics, Cardiff University.

If you have any questions or comments please email Ollie (paino@cardiff.ac.uk).







