# alpha-significance-threshold
Compute significance threshold for a genome-wide association study using estimated effective number of independent tests.

## Usage

The script needs a PLINK formatted genotype dataset (BED/BIM/FAM). 

Rscript num-effective-tests.R -b <plink_file_root> -w <window_size_bp>

The computation proceeds by splitting genome into windows, computing SNP correlation matrix in each window, and performing an eigendocomposition of this matrix. Afterwards, the estimated effective number of SNPs are computed according to 3 methods. Finally, the numbers are added to get a genome-wide estimate.

Larger window sizes lead to slower computations and higher memory use, but may slightly reduce the final number (giving a more correct estimate, since splitting by windows ignores correlations between SNPs in newarby windows).

One heuristic for window size is that the window size should be at least twice estimated LD distance for the genome.

For rice, it should be 0.5Mb or more.

## Output

A CSV file with effective number of independent SNPs estimated according to several methods


## Citations


Li, J., & Ji, L. (2005). Adjusting multiple testing in multilocus analyses using the eigenvalues of a correlation matrix. Heredity, 95(3), 221–227.
https://www.nature.com/articles/6800717


Gao, X., Starmer, J., & Martin, E. R. (2008). A multiple testing correction method for genetic association studies using correlated single nucleotide polymorphisms. Genetic Epidemiology, 32(4), 361–369.

Cheverud JM (2001). A simple correction for multiple comparisons in interval mapping genome scans. Heredity 87: 52–58.

## See also

Galwey, N. W. (2009). A new measure of the effective number of tests, a practical tool for comparing families of non-independent significance tests. Genetic Epidemiology, 33(7), 559–568.

R package poolr:
Cinar, O. & Viechtbauer, W. (2022). The poolr package for combining independent and dependent p values. Journal of Statistical Software, 101(1), 1–42. ⁠https://doi.org/10.18637/jss.v101.i01⁠





