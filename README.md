# SparsePro for fine-mapping with summary statistics and functional annotations

SparsePro is a command line tool for efficient and accurate fine-mapping. For genome-wide fine-mapping, please refer to [SparsePro_gw](https://github.com/zhwm/SparsePro_gw).

## Overview 

<img align="right" src="doc/Fig1.png" width=55% height=55%>

Identifying causal variants from genome-wide association studies (GWASs) is challenging due to widespread linkage disequilibrium (LD) and possible existence of multiple causal variants in the same genomic locus. Functional annotations of the genome may help to prioritize variants that are biologically relevant and thus improve fine-mapping of GWAS results.

To fine-map causal variants, SparsePro takes two lines of evidence. First, with GWAS summary statistics and matched LD information, we jointly infer both the causal status of each effect group for traits as well as variant representation within each effect group. Second, we estimate the functional enrichment of causal variants and prioritize variants according to relevant functional information. As outputs, our method yields variant-level PIP estimates, set-level posterior summary as well as enrichment estimates of functional annotations.

## Installation

SparsePro was developed under Python 3.9.7 environment but should be compatible with older versions of Python 3. The following Python modules are required:

* [numpy](http://www.numpy.org/) (version==1.21.3)
* [scipy](http://www.scipy.org/) (version==1.7.1)
* [pandas](https://pandas.pydata.org/getpandas.html) (version==1.3.4)

To install SparsePro:

```
git clone https://github.com/zhwm/SparsePro.git
cd SparsePro
pip install -r requirements.txt 
gunzip dat/ld.txt.gz
``` 

To test the installation and display basic usage:
```
$> python sparsepro_zld.py -h
usage: sparsepro_zld.py [-h] [--zld ZLD] --zdir ZDIR --N N --save SAVE --prefix PREFIX [--verbose] [--anno ANNO] [--K K] [--pthres PTHRES] [--cthres CTHRES]
                        [--ethres ETHRES] [--aW AW] [--h2]

SparsePro Commands:

optional arguments:
  -h, --help       show this help message and exit
  --zld ZLD        locus finemapping mode: path to matched zscore and ld lists (or annotation file)
  --zdir ZDIR      path to zscores files
  --N N            GWAS sample size
  --save SAVE      path to save result
  --prefix PREFIX  prefix for result files
  --verbose        options for displaying more information
  --anno ANNO      name of the annotation column in zld file
  --K K            largest number of effect groups
  --pthres PTHRES  p value threshold for enrichment
  --cthres CTHRES  coverage level for credible sets
  --ethres ETHRES  entropy level for credible sets
  --aW AW          significant enriched file
  --h2             use previous h2 file as zld file

```

## Input files

Example input files are included in the [dat](dat/) directory.

SparsePro takes in a summary file of loci to be finemapped, z-scores files, LD files and annotations as inputs.

1. **a summary file** contains two mandatory columns: names of z-score file and ld files. Optionally, names for annotation files can be included in the third column. An example can be found at [dat/zldanno.txt](dat/zldanno.txt).

2. **zscore files** that contains two mandatory columns: variant IDs and z-scores. An example can be found at [dat/C1.txt](dat/C1.txt).

3. **LD files** that contains Pearson correlation coefficient matrix. **Please make sure the REF/ALT alleles used in calculating LD are the same as the GWAS study!!** An example can be found at [dat/ld.txt](dat/ld.txt.gz). **We also provide a script for matching alleles at [match_bim_ss.py](match_bim_ss.py)**

4. (optional) **annotation file** with entries indicating annotation status for variants. An example can be found at [dat/anno.txt](dat/anno.txt).


## Usage

### Fine-mapping with GWAS summary statistics only

```
python sparsepro_zld.py --zld dat/zldanno.txt --zdir dat --N 353570 --save result_no --prefix result_no --verbose --K 5 
```

### Fine-mapping with GWAS summary statistics and estimating annotation enrichment

```
python sparsepro_zld.py --zld dat/zldanno.txt --zdir dat --N 353570 --save result_no --prefix result_no --verbose --K 5 --anno anno --pthres 1e-5
```

### Fine-mapping with both GWAS summary statistics and functional annotations selected by G-test

```
python sparsepro_zld.py --zld result_no/result_no.h2.txt --zdir dat --N 353570 --save result_anno --prefix result_anno --verbose --K 5 --anno anno --aW result_no/result_no.W1e-05
```

### Fine-mapping with both GWAS summary statistics and all functional annotations

```
python sparsepro_zld.py --zld result_no/result_no.h2.txt --zdir dat --N 353570 --save result_anno --prefix result_anno --verbose --K 5 --anno anno --aW result_no/result_no.W1.0
```



## Output interpretation

1. **variant-level PIP** (pip) file contain three columns: variant ID, z-scores and PIP. 

```
$> head result_no/C1.txt.pip
rs138862362	0.330197	0.0
rs190750807	0.0329643	0.0
rs9609016	-0.4942	0.0
rs117492340	-2.12801	0.0002
rs572332074	4.38025	0.3072
rs202128203	-1.84605	0.0001
rs55679829	-0.286135	0.0
rs62227035	-1.00044	0.0001
rs5753228	1.15621	0.0001
rs5753229	1.50095	0.0001
```

2. **set-level summary** (cs) file contains three columns: the cs column contains variants included in each effect group; the pip column contains variant representation probabilities in each effect group; the last column contains corresponding effect sizes.

```
$> head result_no/C1.txt.cs
cs	pip	beta
rs557364786	0.9887	0.0118
rs117728004	1.0	0.0107
rs80055673/rs117981957	0.904/0.0546	0.0092/0.0084
rs182440662/rs183772757	0.9135/0.0701	0.0085/0.0077
```

3. **G-test statistics** (wsep) file contains eight columns of G-test statistics.

```
$> head result_no/result_no.wsep
index	k0	k1	r0	r1	W	W_se	p
Conserved_LindbladToh	39069.0	2277.0	30.96	13.47	2.01	0.33	5.58e-07
DHS_Trynka	33750.0	7596.0	14.6	29.84	2.21	0.32	3.27e-12
H3K27ac_Hnisz	14085.0	27261.0	3.15	41.28	1.91	0.58	4.66e-05
H3K4me3_Trynka	32076.0	9270.0	11.42	33.01	2.3	0.34	7.09e-13
Transcr_Hoffman	15795.0	25551.0	14.38	30.05	0.26	0.32	5.15e-01
TSS_Hoffman	39492.0	1854.0	35.59	8.85	1.67	0.37	4.81e-04
UTR_3_UCSC	39789.0	1557.0	41.91	2.52	0.43	0.65	7.91e-01
UTR_5_UCSC	40671.0	675.0	43.88	0.55	-0.28	1.36	1.00e+00
non_synonymous	40878.0	468.0	37.29	7.14	2.82	0.41	1.80e-06
```


4. **joint estimates of enrichment weights** (W) file contains nine columns of joint estimates of enrichment weight for relevant annotations in the [annotation file](dat/anno.txt).

```
$> head result_no/result_no.W1e-05 
index	k0	k1	r0	r1	W	W_se	p	sigidx
Conserved_LindbladToh	32956.15	8389.85	30.96	13.47	0.54	0.33	1.59e-01	0
DHS_Trynka	28604.32	12741.68	14.6	29.84	1.52	0.32	1.46e-06	1
H3K4me3_Trynka	25264.55	16081.45	11.42	33.01	1.51	0.34	3.65e-06	3
non_synonymous	40256.56	1089.44	37.29	7.14	1.96	0.41	3.47e-04	8
```

5. **heritability estimates** (h2.txt) file contains heritability estimates used for specifying prior and can be provided to `--zld` with `--h2`.

```
$> head result_no/result_no.h2
z	ld	anno	h2	pval	varb	K
C1.txt	ld.txt	anno.txt	5.31e-04	1.13e-12	1.43e-04	5.00e+00
C2.txt	ld.txt	anno.txt	5.25e-04	3.00e-14	1.63e-04	5.00e+00
C3.txt	ld.txt	anno.txt	4.95e-04	6.92e-10	1.08e-04	5.00e+00
C4.txt	ld.txt	anno.txt	5.84e-04	1.42e-14	1.67e-04	5.00e+00
C5.txt	ld.txt	anno.txt	4.94e-04	1.14e-12	1.43e-04	5.00e+00
C6.txt	ld.txt	anno.txt	5.03e-04	1.20e-11	1.30e-04	5.00e+00
C7.txt	ld.txt	anno.txt	5.44e-04	5.43e-13	1.47e-04	5.00e+00
C8.txt	ld.txt	anno.txt	4.65e-04	2.86e-13	1.51e-04	5.00e+00
C9.txt	ld.txt	anno.txt	7.51e-04	2.77e-20	2.41e-04	5.00e+00
```

## Citations

If you use this software, please cite:

[Wenmin Zhang, Hamed Najafabadi, Yue Li. SparsePro: an efficient fine-mapping method integrating summary statistics and functional annotations. bioRxiv 2021.10.04.463133](https://doi.org/10.1101/2021.10.04.463133)

