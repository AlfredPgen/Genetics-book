--- 
title: "A handbook for Computational Genetics"
author: "Alfred Pozarickij"
date: "2019-04-21"
site: bookdown::bookdown_site
output: bookdown::gitbook
documentclass: book
bibliography: [book.bib, packages.bib]
biblio-style: apalike
link-citations: yes
github-repo: rstudio/bookdown-dem
description: "This is the book describing theoretical and practical approaches in analysis of genome data."
---

# Preface {-}

Placeholder



<!--chapter:end:index.Rmd-->

# (PART) Quantitative Genetics {-}

# Population parameters

## Mean

## Variance

## Covariance

## Genetic correlation

## Additivity

## Dominance/Recesivness

## Codominance

## Infinitesimal model

### Omnigenic model

## Henetic relationships by Malecot


## Genotype simulations

## Phenotype simulations

<!--chapter:end:01-Descriptives.Rmd-->

# Sequencing technologies

The first step in any genetic analysis is to map sequence reads, callibrate base qualities, and call variants. 

Prior to mapping, evaluate base composition along reads. Calculate the proportion of A, C, G, T bases along each read. Flag runs with evidence of unusual patterns of base composition compared to the target genome.
Evaluate machine quality scores along reads. Calculate average quality scores per position. Flag runs with evidence of unusual quality score distributions.
Calculate the input number of reads and number of bases for each sequenced sample

## Genotype calling algorithms

## Sequence alignment

Sequence alignment is a method of arranging sequences of DNA, RNA, or protein to identify regions of similarity. The similarity being identified, may be a result of functional, structural, or evolutionary relationships between the sequences.

If we compare two sequences, it is known as pairwise sequence alignment. If we compare more than two sequences, it is known as multiple sequence alignment.

## Sequence assembly

## SNP annotation

### CNV annotation

## Gene prediction



<!--chapter:end:02-Sequencing-technologies.Rmd-->


# Genome-wide association analysis

Placeholder


## DNA processing quality control
## Batch effects
### Calculation of missing call rate for samples and SNPs
### Calculation of missing call rates by batch
### Testing for allele frequency differences in batches
## Sample quality control
### Cryptic relatedness
### Population stratification
### Heterozygosity and missingness outliers
### Differential missingness
### Sex chromosome anomalies
## Marker quality control
### Genotyping concordance
### Mendelian errors
### Genotype call rate
### Minor allele frequency
### Hardy-Weinberg equilibrium outliers
### Additional QC for regions like MHC
### Ambigious nucleotides
### Non-matching nucleotides
### Quality control prior to meta-analysis
## X-chromosome quality control
## Single marker regression
### Allelic test
### Genotypic test
### Additive model
### Dominant model
### Recessive model
### Categorical phenotype
### Multi-allelic GWAS
## Two-stage approach
## Haplotype GWAs design
## Joint analysis (all independent markers)
### Genomic control
#### λ<sub>1000<sub>
##### Stratified λ<sub>GC<sub>
## Multimarker single gene-based approaches
## VEGAS
## Multimarker gene-set approaches (a.k.a. pathway analysis)
## fastBAT
## MAGMA
## VEGA
## Extensions to binary and categorical phenotypes
### Threshold model
## Analysis of rare variants
### Collapsing methods based on pooling multiple rare variants (burden or adaptive burden tests)
#### Sum test
#### Cohort Allelic Sums test (CAST)
#### Combined Multivariate Collapsing (CMC)
#### Weighted Sum test (WSS)
#### Kernel Based Adaptive Cluster (KBAC)
#### Replication Based Test (RBT)
#### ARIEL test
#### The EREC test
### Methods treating rare variant effects as random (Variance-component tests)
#### The SSU approach
#### C-alpha test
#### SKAT
### Methods based on model selection
#### Seq-aSum
#### Seq-aSum-VS
#### Variable Threshold Test (VT)
#### RARECOVER
#### Selective grouping method
#### Step-Up
### Combination of collapsing and random effects methods
#### SKAT-O
#### SKAT-C
#### Fisher method
#### MiST
### EC test
### Family-based tests
## Analysis of X, Y and mitochondrial chromosomes
### Dosage compensation
## Analysis of copy number variants
### Common variation
### Analysis of rare variants
## Analysis of multi-ethnic samples
## Analysis of indirect genetic effects
## Exome analysis
## Whole-genome analysis
### Deep whole genome sequencing
### Low coverage whole genome sequencing
## Analysis of multiple traits
## Mixed-model association analysis
### EMMAX
### Fast-LMM
### GEMMA
### BOLT-LMM
### Caveats
## Penalized regression GWAS
## Bayesian GWAS
## Machine learning for GWAS
## Expected increase in GWAS loci with sample size
## The joint effect of genotypes over all traits
## Adjustment for winner's curse
## Adjustment for assortative mating
## Adjustment for attenuation bias
## Enrichment in candidate genes

<!--chapter:end:03-Genome-wide-association-analyses.Rmd-->


# Heritability

Placeholder


## Realized heritability
## Liability vs observed scale
### Evolvability
### Reliability
## Twin studies
## Haseman-Elston regression
## GREML
## GREML in family data
## GREML in WGS or imputed data
### GREMLd
### Bivariate GREML
## LD-score regression (LDSC)
## LDAK
## Heritability by chromosome, MAF bin, or functional category

<!--chapter:end:04-Heritability.Rmd-->


# Genomic prediction

Placeholder


## Unweighted sum of risk alleles
## Polygenic risk scores
## Gene-based polygenic score (POLARIS)
## Pathway-based polygenic risk score
## LD adjusted PRS
### LDpred with functional annotation
## Annopred
## Pleiopred
## Prediction including GxE
## BLUP
### GBLUP
### sBLUP
## Bayesian Zoo
### B
### C
### S
### N
### NS
### R
## Reproducing kernel Hilbert space
## Machine learning methods

<!--chapter:end:05-Genomic-prediction.Rmd-->

# Pleiotropy

## Fisher's geometric model

## Direct

## Indirect


<!--chapter:end:06-Pleiotropy.Rmd-->

# Gene-set analysis

## Gene-set conditional analysis

## Gene-set interaction analysis

<!--chapter:end:07-Pathway-analysis.Rmd-->

# Functional annotation

<!--chapter:end:08-Functional-annotation.Rmd-->


# Causal inference

Placeholder


## Gene-knockout
## Conditioning
## COJO
## mtCOJO
## Finemapping
## Mendelian Randomization
### Summary data-based MR (SMR)
### Generalised summary-data-based MR (GSMR) and HEIDI
### Joint analysis of GWAS and eQTL data
### Tissue-specific MR

<!--chapter:end:09-Inferring-causality.Rmd-->

# Combining multiple datasets

## Meta-analysis

### Meta-analysis of gene-level associations (common)?

### Meta-analysis of rare variants

RAREMETAL and RAREMETALWORKER

## Mega-analysis


## Z-statistic to estimated SNP effect

After sample-size-weighted meta-analysis, Z-statistics can be transformed into unstandardized regression coefficinets using the following equation:

$$\LARGE \hat{\beta_j} = Z_j\frac{\hat{\sigma_Y}}{\sqrt{2N_j MAF_j(1-MAF_j)}}$$
for SNP *j* with minor allele frequency MAF*j*, sample size N*j*, Z-statistic Z*j*𝑍𝑗, and standard deviation of the phenotype $\hat{\sigma_Y}$𝑌. 

<!--chapter:end:10-Analysis-of-multiple-datasets.Rmd-->

# Gene-environment interaction

It is generally accepted that complex diseases are caused by an interplay of genetic and environmental factors, creating a challenge for understanding the disease mechanisms. Understanding the interplay between genes and environmental factors is important, as genes do not operate in isolation but rather in complex networks and pathways influenced by environmental factors.
Identification of gene-environment interactions has important implications for understanding underlying disease etiology and developing disease prevention and intervention strategies.
In addition to providing insights into disease etiology, exploiting gene-environment (G-E) interaction can help discover novel susceptibility loci for complex diseases, where genetic effects are modified and masked by the effects of environmental factors.
From a public health perspective, G-E interaction is useful because findings based on interactions can help develop strategies for targeted intervention; conducting an intervention focusing on a subset of the population identified by G-E interactions can provide efficiency in disease prevention.
There are several challenges of G-E interaction analysis that include replication issues. While more powerful statistical methods for detecting interactions are helpful, ultimately studies with larger sample sizes are needed to identify interactions through consortium-based studies to achieve adequate power for G-E analysis.

This chapter aims to describe methods that extend single marker regression by taking into consideration the effect of gene-environment interaction.

## Single step methods

There are several disease risk models for the joint effects of G and E, and interpretations of G-E interactions depend on the underlying disease risk models.

### Case-control

### Prospective likelihood-based approach

### Retrospective likelihood approach

To address the limitations of case-only approaches that can only test for multiplicative interactions (not for the main effects of G and E), Umbach and Weinberg (1997) generalized the case-only design idea to use a log-linear model based on case-control data. They showed the maximum-likelihood estimates for all parameters of a logistic regression model can be obtained using a log-linear model. Along the same line, Chatterjee and Carroll developed a general method using a retrospective likelihood that exploits the G-E independence assumption to test for multiplicative interaction, but can use both cases and controls to estimate all of the parameters in a general logistic regression model. Basically, this method employs a retrospective likelihood that explicitly models the conditional probability of G given E mediated by an association parameter θ that can be constrained to be zero when the G-E independence assumption holds. This likelihood can be used for testing both multiplicative and additive interactions; recently, Han et al. developed a likelihood ratio test that exploits the G-E independence assumption using a retrospective likelihood. Their numerical investigation of power suggests that the incorporation of the independence assumption can enhance the efficiency of the test for additive interaction by 2- to 2.5-fold.

#### Multiplicative scale

A multiplicative model incorporating GxE effects is one of the most commonly used models via logistic regression: logit(Pr(*D* = 1|*G*,*E*)) = *β<sub>o</sub>* + *β<sub>g</sub>G* + *β<sub>e</sub>E* + *β<sub>ge</sub>GE*, where *G* is a genotype of a SNP, *E* is an environmental risk factor, and *D* is the disease status. 

Assuming binary factors for both G and E, a 2 × 2 table for a disease risk for each combination of G and E values can be constructed based on this model (Table X).

#### Additive scale

An additive model is shown as logit(Pr(*D* = 1|*G*,*E*)) = *b<sub>o</sub>* + *b<sub>g</sub>G* + *b<sub>e</sub>E* + *b<sub>ge</sub>GE*, where the effects of G and E are additive on the disease risk scale, but not on the logit scale. 

### Case-only

This design depends on an assumption of G–E independence in the underlying population. Under the assumption of G-E independence in the underlying population (i.e., controls), a multiplicative interaction test statistic becomes equivalent to testing the association between G and E among cases. One major limitation of the case-only design is that while the case-only method has improved power over the traditional methods when G and E are independent in the underlying population, this method has an increased type I error if the independence assumption is violated. In addition, the regression parameters for the main effects of G and E cannot be estimated using this method because the case-only test is only for evaluating a multiplicative interaction.

### Empirical Bayes and Bayesian Model Averaging

### Other Bayesian approaches

#### Interactions using haplotypes

First, haplotypes are biologically relevant. There is strong evidence that several mutations within a gene may interact ( cis-interaction) to cause a disease. Second, the power of single-marker-based methods in association studies depends on LD between the tested marker and the disease susceptibility locus. LD information contained in flanking markers is generally not incorporated, which can result in a reduction in power. In addition, haplotype-based methods can be more powerful when multiple disease susceptibility alleles occur within a single gene.

Lake et al.  [21]  proposed a likelihood-based method in the generalized linear model framework, which has been commonly employed in haplotype-based association studies, because it is publicly available and easy to implement with its R package. This approach, however, is limited by ignoring the interacting effects between haplotype blocks. Subsequently, several methods have been developed to study haplotype-related interactions, but these methods do not consider all potential haplotypes and interactions simultaneously  [22–27] . Recently, Guo and Lin  [15]  proposed a generalized linear model with regularization to detect interacting haplotype effects. However, their method applies an omnibus test and consequently does not provide inference on the effects of specific haplotypes and their interactions. Another concern for haplotype-based methods is that large numbers of haplotypes inferred from genotype data  [28–30]  often lead to high degrees of freedom for corresponding statistical tests and thus reduce power  [14, 15, 31–34].  If interacting effects are considered, such problems become severer. Unfortunately, few statistical methods have been developed to tackle these farreaching problems in the study of haplotype interactions.  One challenge in haplotype-based association studies is that for haplotypes comprising multiple markers, there might be many rare haplotypes. Because of their low frequencies, the parameter estimates related to rare haplotypes will have large variances, leading to unstable models. Schaid  [6]  described several approaches to handle rare haplotypes. One approach is to combine all rare haplotypes into one group or rare haplotypes with common ancestral haplotypes, and another is to exclude them from the model, which is equivalent to including them in the baseline group. However, both approaches yield results that may be difficult to interpret. In addition, it has been argued that rare haplotypes may account for a substantial fraction of the multifactorial inheritance of common diseases  [35–39],  thus the aforementioned approaches may miss the rare haplotypes having true effects. Another approach is to include the effects of each rare haplotype in the model but shrink their effects towards the common mean or towards the effect of a similar haplotype  [6, 15, 40].  Recently, Guo and Lin  [15]  adopted a lasso penalty in their generalized linear model to allow assessment of the effects of rare haplotypes by decreasing the coefficients of unassociated haplotypes to zero so that the associated ones, especially those that are rare, can be estimated. 


The estimate of haplotype dosage is the estimate of the number of copies of a specific haplotype for a subject. For the haplotypes that can be unambiguously resolved based on the observed genotype data, the values of haplotype dosage of a haplotype for a subject can be 0, 1, or 2. But for the haplotypes that cannot be unambiguously resolved, the values of haplotype dosage of a haplotype for a subject would be non-integer, ranging from zero to two, reflecting the possibility that haplotypes are based on the subject’s genotypes. For each subject, the sum of haplotype dosage across all haplotypes is equal to two.
 
## GxE in the context of family studies

## Multi stage methods

Several approaches exist. In general, these methods suggest selecting a subset of SNPs based on the marginal effects of SNPs or G-E correlation tests in the first stage and conducting standard G-E interaction tests in the second stage, where the independence between the test statistics used in the two stages is required to provide a valid screening procedure.

For calculating power for G-E interactions, the powerGWASinteraction R package is available (https://cran.r-project.org/web/packages/powerGWASinteraction/index.html), which includes a power calculation tool for four two-stage screening and testing procedures.

## Joint tests

### Gene-based

## Set-based interaction tests

There are multiple reasons for using set-based gene-environment interaction tests. 

1. Multiple comparison adjustments for a large number of markers across the genome could result in power loss.
2. Closely located SNPs are correlated because of linkage disequilibrium. Multiple tests for GxE in these single-marker-based GxE models are even more dependent, as interaction terms in these models share the same environmental variable. Dependence among multiple tests can result in incorrect Type 1 error rates and causes bias in standard multiple comparison adjustments and this bias is often difficult to correct.
3. The single-marker GxE test does not interrogate the joint effects of multiple SNPs that have similar biological functions. When the main effects of multiple SNPs in a set are  associated with a disease/trait, the classical single marker regression interaction test can be biased.

Lin et al. (2013) developed a method to analyse GxE for a set of markers using generalized linear mixed models. The method tests for SNP-set by environment interactions using a variance component test, and because a set of variants will likely be correlated, the main SNP effect estimates under the null hypothesis are obtained using ridge regression. 
Their software is called GESAT. Here, they model GxE effects as random, as opposed to the classical approach of treating *β<sub>j</sub>*’s as fixed effects followed by a test with *p* degrees of freedom. The latter approach can suffer from power loss when *p* is moderate/large, and numerical difficulties when some genetic markers in the set are in high LD. 
The model allows to adjust for the main effects of all SNPs while simultaneously testing for the interactions between the SNPs in the region and environmental variable. For unbalanced designs when a binary environmental exposure has a low frequency in one category, GESAT is most advantageous over single marker regression GxE test. Such unbalanced designs can occur due to case–control sampling and the strong association of an environmental factor with disease. When the effect size is modest, GESAT performs better that single marker regression GxE test, but when the effect size is strong, the opposite is true. 
Their simulations suggest that the power of GESAT seems fairly robust to the dependence between G and E.

The same approach can be applied to investigating various other biological problems. For example, we can test for the interactions between gene expressions in a pathway or network and an environmental variable by simply replacing G by gene expressions in a gene-set.

Existing methods for assessing common variants by environment interactions such as Gene-Environment Set Association Test (GESAT) (Lin et al., 2013) have several limitations when applied for rare variants. GESAT estimates the main effects of the common variants by applying a L<sub>2</sub> penalty on the genotypes scaled to unit variance; this assumes that the main effects of the scaled genotypes are comparable in magnitudes, which may not hold in the case of rare variants. GESAT also assumes that the regression coefficients of the rare variants by environment interactions are independent of each other, and suffers from power loss when most rare variants in a gene interact with the environmental factor and the interaction effects have the same direction.

Similar idea from GESAT was later extended and applied to analyse rare variants. GESAT have several limitations when applied for rare variants. GESAT estimates the main effects of the common variants by applying a L<sub>2</sub> penalty on the genotypes scaled to unit variance; this assumes that the main effects of the scaled genotypes are comparable in magnitudes, which may not hold in the case of rare variants. GESAT also assumes that the regression coefficients of the rare variants by environment interactions are independent of each other, and suffers from power loss when most rare variants in a gene interact with the environmental factor and the interaction effects have the same direction. The new proposed test iSKAT is optimal in a class of variance component tests and is powerful and robust to the proportion of variants in a gene that interact with environment and the signs of the effects. This test properly controls for the main effects of the rare variants using weighted ridge regression while adjusting for covariates.

A naive approach to assess rare variants by environment interactions is to extend the burden test by fitting a model with both the summary genetic burden variable, environment, and their interaction, and performing a one degree of freedom test for the interaction. However, when there are multiple causal variants with their main effects having different magnitudes and/or signs, such a burden rare variant by environment test fails, and may lead to inflated Type 1 error rates. This is because adjusting for the main effects of the multiple causal variants using a single summary genetic burden variable is inappropriate. Likewise, a naive approach to assess rare variants by environment interactions using SKAT by including the main effects of rare variants as part of covariates and applying SKAT to the interaction terms is problematic. This is because SKAT only allows adjustment of a small number of covariates and cannot handle the presence of a large number of rare variants in a region. Furthermore since the rare variants are observed in low frequency, a model with all the rare variants as main effects will be highly unstable and may not even converge.

Both GESAT and iSKAT are able to incorporate multiple environmental variables.

No gene-based GxE test exists for the analysis of the sex chromosomes as of writing this text.

The rareGE R package (https://www.hsph.harvard.edu/han-chen/software/) provides various functions for detecting G-E interaction as well as for testing the joint effect of a gene and G-E interaction under a set-based framework. 
The SIMreg R package (http://www4.stat.ncsu.edu/~jytzeng/software_simreg.php) offers functions for testing a set-based G-E interaction by using genetic similarity to aggregate information across SNPs, and incorporating adaptive weights depending on allele frequencies to accommodate rare and common variants.

https://onlinelibrary.wiley.com/doi/abs/10.1111/biom.12428

https://academic.oup.com/biostatistics/article-abstract/18/1/119/2555344?redirectedFrom=fulltext

### Combining multiple environmental factors

### Multi-trait multi-GxE tests

### Variance heterogeneity

All methods described in this chapter so far conduct a location (mean) based test. In other words, the main concern is whether there is a difference in the mean effect between genotypes. Testing for scale (variance) heterogeneity could be used to identify novel variants involved in the control of the phenotype possibly through interactions with other genetic variants or environmental factors. This is advantageous when the interacting variable was not collected and the interaction term may not be directly modeled or when a phenotype was transformed. Therefore, the analysis of variances of the trait is based on SNP information only. The presence of interaction between a genotype and certain factor is expected to alter the trait's variance in the group of subjects carrying such genotype.

Yang et al. (2012) noticed that the mapping of variance-controlling loci is prone to inflated test statistics when the minor allele frequency (MAF) is small.

One way of interpreting an increase in variability is as a decrease in stability. Waddington (1942) described the concept of canalization, whereby natural selection favors the relative constancy of some attributes, for example, well-formed organs and limbs, and thereby leads to the evolution of heritable architectures that buffer the impact of environmental or background genetic variation that would otherwise cause development to go astray. These architectures create virtual “canals” down which developmental programs flow. For a canalized phenotype, which modern usage expands to include nondevelopmental traits, the “zone of canalization” is the range of underlying liability over which potentially disruptive variation may be absorbed without serious consequence to the expressed trait value (Lynch and Walsh 1998). A well-studied example of a stabilizing architecture is that provided by heat-shock protein 90 (Hsp90), which buffers genetic and stochastic variation in the development of plants and flies (Rutherford and Lindquist 1998; Queitsch et al. 2002; Sangster et al. 2008).

But in absorbing variation, such stabilizing architectures also hide it from view, and a sensitizing change in the stabilizer that shifts liability outside the zone of canalization can have a dramatic effect on the phenotype. Such shifts release the combined effects of previously “cryptic” genetic variation: now decanalized, the phenotype is more sensitive to internal (including genetic) and external environment, and as a result varies more greatly between individuals (Dworkin 2005; Hornstein and Shomron 2006). In this vein, decanalization has been proposed to explain why the genetic architectures of some diseases in human populations seem more amenable than others to genetic dissection through genome-wide association (Gibson and Goldstein 2007). Specifically, whereas some disease phenotypes in homogeneous populations may be heavily canalized and thereby harder to dissect, others may have been decanalized by modern living conditions (e.g., inflammatory diseases) or modern admixture, while yet others are simply too recent in evolutionary history for buffering networks to have evolved (e.g., response to HIV). Increased variability can also be adaptive. In natural populations disruptive selection favors diversity, with increased “capacitance” (Rice 2008) or “bet-hedging” (Beaumont et al. 2009) spreading risk over a variable fitness landscape. Feinberg and Irizarry (2010) recently proposed a heritable and selectable mechanism for this based on stochastic epigenetic variation. In controlled populations, variability can be increased through directional selection. For example, in a Drosophila selection experiment Clayton and Robertson (1957) reported increased bristle number variance, which is consistent with the idea that genotypes associated with higher environmental variance have a greater chance of being selected under directional selection (Hill and Zhang 2004). Moreover, genetic differences have been observed for phenotypic variability in body weight for chickens (Rowe et al. 2006) and snails (Ros et al. 2004) and litter size in rabbits (Ibanez-Escriche et al. 2008), sheep (Sancristobal-Gaudy et al. 1998), and pigs (Sorensen and Waagepetersen 2003). In natural populations with stabilizing selection we should expect to find alleles minimizing variance for fitness traits (Lande 1980; Houle 1992), whereas directional selection during domestication will favor alleles that increase variance. One may therefore expect to find vQTL in experimental crosses between wild and domestic animals (see Andersson 2001). Nonetheless, genetic buffering that leads to phenotypic robustness need not require an evolutionary explanation to be observed, nor to be useful in medicine and agriculture. Plainly, detecting vQTL and inferring how they arose are separate questions; here we concentrate on the first.

Given that a genotypic effect at one locus may be affected by the genotype at an adjunct locus or some environmental factor, the genotypic effect at that locus is a composite of multiple distributions with different means, resulting in an inflated variance for that genotype (Struchalin et al., 2010; Deng and Par´e, 2011). In addition, it has been shown that vQTL can be induced by linkage disequilibrium with nearby functional locus with mean effects only (Cao et al., 2014).

### Levene's test

Levene’s test (Levene, 1960) is known for its simplicity and robustness to modeling assumptions, and it is perhaps the most popular method for evaluating variance heterogeneity between *k* groups. It's performance is best when the genetic loci have only variance effects on quantitative traits (Hong et al., 2016). More recently it has been employed as an indirect test for interaction effects.

### Brown-Forsythe test (also known as Levene's median test)

$$\LARGE T^2 = \frac{(N-k) \sum_{j=0}^{k-1} n_j(Z_j.-Z_..)^2}{(k-1)\sum_{i=1}^{N}(Z_i-Z_gi.)}$$

The Brown-Forsythe is a statistical test for the equality of group variances and is based on
an ANOVA of the absolute deviation from the median. It has earlier been shown to be robust to deviations from normality of the phenotypic distribution in GWAS applications.

### Bartlett's test

$$\LARGE T^2 = \frac{(N-k)ln(\sigma_P^2)- \sum_{j=0}^{k-1} (n_j-1)ln(\sigma_j^2)}{1+\frac{1}{3(k-1)}(\sum_{j=0}^{k-1}(\frac{1}{n_j-1}-\frac{1}{N-k}))}$$
where *k* is the number of genotypes tested,

### Bartlett's test with prior rank transformation to normality

### Generalized Levene's scale tests

The Levene's test for variance heterogeneity was expanded to include sample correlation and group membership uncertainty (genotype imputation probabilities). Following a two-stage regression framework, it was shown that the least absolute deviation regression must be used in the stage 1 analysis to ensure a correct asymptotic χ<sup>2</sup><sub>k−1</sub>/(k−1) distribution of the generalized scale (gS) test statistic. 
gS has good type 1 error control in the presence of sample correlation, small samples, unbalanced group sizes, and non-symmetric outcome data.

Same study showed that the least absolute deviation (LAD) regression approach to obtain group-median-adjusted residuals is needed to ensure robust performance of gS (and possibly LAD should be used in all scale dependent tests).

Several improvements including adjustment for covariates could be further explored (check original paper)

### The Lepage test

The Lepage test is a rank-based non-parametric test for either mean or variance difference, which combines the Kruskal–Wallis test statistic for mean difference and the Fligner–Killeen test statistics for variance difference, and has been shown to be powerful and robust when the locus has both mean and variance effects (Lepage, 1971; Hollander and Wolfe, 1999).

### The *D*-test

Aschard et al. (2013) introduced a general distribution-based test or *D*-test, which can be thought as an extension of the classical Kolmogorov–Smirnov test (Brittain, 1987). There are two versions of *D*-test, namely, the constrained *D*-test (*D*<sub>c</sub>) and the unconstrained *D*-test (*D*<sub>u</sub>). The *D*<sub>c</sub> test is designed for the situation when the effect of the variant is monotonic, while the *D*<sub>u</sub> test has the ability to detect genetic effects in a broader range of situations (Aschard et al., 2013).

### Regression using the squared Z-score

### Gamma regression models

### Two-step screening on residual variance heterogeneity

#### YGVH

One caviet involving VH tests is that the absence of VH for a SNP can not be interpreted as the absence of the SNP involvment of the SNP in the interaction network.

When there are only two genotype classes, the type I error rate can be very large if the number of minor genotypes contains fewer than 100 observations when using regression on the squared Z-score, and this cannot be overcome by increasing the total sample-size. The Levene and Brown–Forsythe tests also show such an inflation of false positives, but use of a Gamma regression model, which accounts for the fact that the squared Z-score follows a chi-square distribution, overcomes this problem. Populations with three genotypes will, in practice, be more robust when the allele substitution model implemented in most GWAS-software is used (i.e., when regression on all three genotypes is used to estimate the additive effect). Inflated type I error rates are then observed only when the intermediate-size genotype class (i.e., in practice most often the heterozygotes) contains fewer than 100 individuals.

## Mean-variance QTL (joint test)

It has been suggested that many loci have both mean and variance effects, while some of the mean or variance effects alone would not be strong enough to be detected (Shen et al., 2012; Cao et al., 2014).

One approach to detect both mean and variance effects includes double generalized linear model (DGLM). DGLM can be used to model effects of the mean and variance simultaneously, thus enhacing the power to detect variants associated with the trait. This approach is able to not only incorporate genetic and covariate effects for the mean but also set of such effects for the residual variance.

Extensions to family-based data have been developed (Cao et al., 2014).

### Conditional quantile regression 

In contrast to OLS regression, quantile regression tests for the effect differences across the range of quantiles. More specifically, conditional quantile regression (CQR), can be used to model effect size changes across the sample distribution. This method is similar to other varaince heterogeneity methods in a way that does not require interacting effects to be known and hence included in the model. The downside of this is that further investigation is required to identify environmental factors causing statistical interactions.

If the marginal and interaction effects are in opposite directions, then the power to detect potential interactions is reduced. Hence, failure to identify differences in variance by genotype effects does not rule out interactions effects.

Tests for variance heterogeneity detect variance inconsistency rather than variance structure (i.e. direction of change). Therefore, variance heterogeneity is not specific to interaction signals, but includes conditions where no consistent direction for increasing or decreasing variance per genotype is observed. Modeling relationship between variance and genotypes assuming a structure (i.e. linear trend) could help improve power to detect variants with potential interactions if the assumtions are met.

### Sliced inverse regression

Jiang and Liu, 2014

### Semiparametric model for vQTL mapping

The semiparametric model is designed to identify the combined differences in higher moments among genotypic groups. The pairwise likelihood is constructed based on conditioning procedure, where the unknown reference distribution is eliminated.

### Variance heterogeneity for related individuals

## RELIEF and other machine learning tools

### Multidimensionality reduction

## Meta-analytic GxE approaches

## Gene-environment correlation

## Other challenges

There are several challenges of G-E interaction analysis. One main challenge is replication issues. While various GWAS findings of the main effects of SNPs have been replicated by independent studies for many complex diseases (http://www.ebi.ac.uk/gwas/), relatively few interactions have been reproduced. It is likely that the sample sizes of GWAS that have required measurements on environmental exposures are not yet adequate to reliably identify G-E interactions of modest magnitude. In addition, differences in the underlying distribution of environmental exposures across various studies as well as difficulties in accurately measuring environmental exposures can also lead to reduced power of detecting G-E interactions. While more powerful statistical methods for detecting interactions are helpful, ultimately studies with larger sample sizes are needed to identify interactions (e.g., through consortium-based studies) to achieve adequate power for G-E analysis. A reasonable goal for the future will be to at least identify parsimonious models that adequately describe the risks of diseases associated with a combination of genetic and environmental risk factors. The lack of reporting of interaction in current studies so far indicates that linear logistic models, i.e., multiplicative models, in general may be a good starting point for building models for evaluating the joint effects of genetic and environmental factors.

When using a test designed to detect variance heterogeneity, it is of vital importance to ensure that phenotype transformations do not lead to false positive associations. For an autosomal diallelic SNP in particular, it is difficult to make valid biological inferences about a quantitative trait on the basis of variance heterogeneity among the genotype-specific distributions unless the locations (in particular, the means) of the three distributions are equal. When the means are not all equal, the variance heterogeneity might in fact be explained by a distribution in which the variance is a function of the mean. There exists a transformation of the data that will tend to make the variances equal and there is a method that will at a minimum decrease the variance heterogeneity.

GWAS analyses to detect variance heterogeneity is inherently sensitive to unbalanced data. The major reason for this is that the distribution of the variance often deviates from normality as it: 
1) is bounded at zero 
2) has a distribution skewed to the right 
3) has a variance depending on its mean 
Such deviations lead to violations of, e.g., the Gauss–Markov assumptions in a regression model.

One disadvantage of using variance heterogeneity methods is that the underlying interacting factors are unknown. Given that shifts in variance could be caused by the incomplete LD between the causal polymorphism and the tested marker, multiple functional alleles, gene-gene or gene-by-environment interactions, it could be difficult to identify the root cause of phenotypic variability. Transformations on a phenotype can also result in variance heterogeneity (Sun et al., 2013). This transformation can occur knowingly for statistical purposes, for example, log(Y), or unknowingly, for example, choosing a phenotype measurement that does not directly represent the true underlying biological outcome of a gene.

### Controlling for covariate interactions

In order to eliminate potential confounders, it is common to add covariates such as sex and ethnicity in the model. While this provides control for their potentially confounding influences on the main effects of the genotype and the environment, in the case of GxE studies, this practice does not control for the effects these variables might have on the G×E interaction. Rather, to properly control for confounders, researchers need to enter the covariate-by-environment and the covariate-by-gene interaction terms in the same model that tests the G×E term.

The GxE interaction will be biased in the model without the covariate-by-environment and the covariate-by-gene interaction terms if:

1) The covariate is related to the genetic variable and the covariate-by-environment interaction coefficient is nonzero.
2) The covariate is related to emvironmental variable and the covariate-by-gene interaction coefficient is nonzero.

Finally, it should be noted that even if a G×E result ‘disappears’ after properly controlling for covariates, this does not necessarily mean that the original G×E hypothesis was wrong. For example, the genetic polymorphism might cause changes in the covariate which in turn moderates the environmental variable, in which case the covariate is a mediating mechanism by which the gene moderates the environmental variable. That said, this possibility applies to all models that statistically control for covariates in regression, and the traditional interpretation of ‘disappearing’ effects after controlling for a covariate is that the true causal pathway is ambiguous and alternative (confounding) explanations cannot be ruled out. That said, in some cases, a particular causal pathway can be discarded as impossible or unlikely. In such cases, investigators can be more definitive about ruling out certain hypotheses. For example, changes at a genetic polymorphism will not lead to changes in ethnicity, and so a G×E hypothesis can be safely discarded if it is mediated by an ethnicity-by-environment interaction.

## Candidate gene-by-environment interaction studies

This research area has been a hot topic in genetics, with hundreds of publications reporting positive G×E discoveries over the last 15 years, but there has been increasing skepticism about the validity of many of these findings. This skepticism is based on a number of substantive and statistical concerns:

1) A low replication rate among attempted direct replications of GxE findings.
2) The possibility that GxE findings capitalized on chance from among many unreported analyses.
3) A publication bias toward positive findings.
4) Small sample sizes that exacerbate the already low statistical power for detecting interactions.
5) The low prior probability that a specified environmental variable interacts with a specified candidate gene polymorphism.

## Higher order interactions

It is difficult and potentially misleading to interpret two-way interactions in the presence of three-way interactions. In such a model, the lower-order two-way interactions become conditional interactions, and the regression betas and p-values are interpreted as the predicted two-way interactions when the other (omitted) variable is coded as 0.




<!--chapter:end:11-GxE.Rmd-->

# Gene-gene interaction

When the combined phenotypic effect of alleles at two or more loci deviates from the sum of their individual effects, this is referred to as a genetic interaction, or epistasis. 
There are some situations where data and theory have suggested that it might be particularly important to account for genetic interactions. One is when the aim is to predict the phenotypes of individuals on the basis of their genotype. If interactions lead to extreme phenotypes for some genotypes, these phenotypes are unlikely to be captured by additive models, particularly if they are rare. Another case is the prediction of long-term selection response. Under an additive model, both the additive variance and the response are expected to be nearly constant over the first few generations. As generations proceed, allele frequencies change to alter the additive variance and, consequently, the response to selection. This change is more rapid for traits regulated by fewer loci with larger effects than for traits regulated by many loci with smaller effects. It is known that genetic interactions can contribute to the additive genetic variance in a population. The contribution, however, varies depending on the joint allele frequencies across all the interacting loci as well as on the types and strengths of the genetic interactions. The changes in the additive variance, and hence the response, during ongoing selection are therefore more complex in the presence of genetic interactions.

Most genetic variance in a population is expected to be additive even in the presence of extensive epistasis. The lack of empirical knowledge about the pervasiveness and strength of epistasis in the genetic architectures of complex traits makes it largely unknown how much of the observed additive genetic variance in quantitative genetics studies is due to genetic interactions.

Epistatic gene action, namely when the effect of an allele at one locus varies depending on the genotype at another locus, is therefore not directly proportional to the level of epistatic variance in a population. This is because it will usually contribute to both the additive and epistatic genetic variances (Goodnight, 1988; Cheverud and Routman, 1995; Mackay, 2014). To what extent epistatic gene action will contribute additive genetic variance is determined by allele frequencies, the type of genetic interactions, and how the genetic models used are parameterized (Cheverud and Routman, 1995; Hill et al., 2008; Huang and Mackay, 2016).

## Single step methods

## Multi stage methods

## Machine learning methods

## Variance heterogeneity

The details for this method have been described in the previous chapter.

To identify an individual locus that makes direct contributions to the trait variance, a statistical test is used to identify significant differences in the phenotypic variance between the groups of individuals that carry alternative alleles at the locus. When such a variance difference exists between the genotypes at a locus, the locus displays a genetic variance-heterogeneity. 

The concept that the trait variability could also be under direct genetic control was introduced already in the 1940s when Waddington presented the idea of canalization, where he suggested that natural selection could act to produce traits that are robust to environmental and genetic perturbations (Waddington, 1942). He partly based his ideas on the observation that natural populations often are less variable than artificial populations of the same species. More recently, Hill and Mulder (2010) proposed that ‘the environmental variation’ can be regarded as a phenotype in its own right. One can then invoke much of the quantitative genetics methodology to search for genetic determinants of this phenotype. They consequently called this phenomenon ‘genetic control of the environmental variation’, a terminology which implies that it is the randomness, or instability, of the trait that is genetically controlled. Several studies have recently mapped individual loci where the different alleles affect not only the mean, but also the variance of traits (Hill and Mulder, 2010; Rönnegård and Valdar, 2012; Shen et al., 2012). These loci can be detected since the variability of the measured trait differs between groups of individuals that carry alternative alleles at the locus. A simple example would be two groups of humans, where the group of individuals homozygote for a certain allele include both very short and very tall individuals, while the second group that is homozygote for the alternative allele include individuals of similar height. This would lead to genetic variance heterogeneity between the two groups of individuals. Note that the mean height does not have to be different between the groups in order for this to occur.

<!--chapter:end:12-GxG.Rmd-->

# Other omics

## Transcriptome-wide association studies

### cis eQTLs

### trans eQTLs

### 3-D structure of the genome

## Phenome-wide association studies

## Metabolomics

## Epigenomics

<!--chapter:end:13-Multi-omics.Rmd-->

# Quantitative trait loci mapping



<!--chapter:end:14-QTL-mapping.Rmd-->


# Additional points of interest

Placeholder


## Kinship matrix
### Path coefficients
## Genetic relationship matrix
## Animal models
## Phasing
### Switch rate
## Haplotyping
## Statistical power
### Quanto
### GCTA-GREML
### Mendelian Randomisation
### Twin design
### When marker is a disease susceptibility locus
### When marker is not disease susceptability locus
## Multiple comparisons
### Effective number of independant variants
## Biases
### Ascertainment
### Attenuation
### Selection
## Family studies
### Transmission disequilibrium tests
## Twin studies
## Adoption studies
## Equifinality (many genes give same trait)
## Gene dosage
### Allelic dosage
## Allelic heterogeneity
## Genetic heterogeneity
## Genomic imprinting
## Penetrance/phenocopy
## Endophenotypes
## Ploidy
## Extended phenotype
## Genome sizes
## cis-eQTL vs trans-eQTL

<!--chapter:end:15-Extra.Rmd-->

# (PART) Population Genetics {-}

# Genetic drift

<!--chapter:end:16-Genetic-drift.Rmd-->

# Mutation


## Mutation age

<!--chapter:end:17-Mutation.Rmd-->

# Selection

## Directional

## Balancing

### Frequency-dependent selection 1

### Frequency dependent selection 2


## Background selection

<!--chapter:end:18-Selection.Rmd-->

# Migration

<!--chapter:end:19-Migration.Rmd-->

# Diversity

<!--chapter:end:20-Diversity.Rmd-->

# Admixture

<!--chapter:end:21-Admixture.Rmd-->

# Linkage disequilibrium

<!--chapter:end:22-LD.Rmd-->

# In breeding and heterosis

<!--chapter:end:23-Inbreeding-and-heterosis.Rmd-->

# Assortative mating

<!--chapter:end:24-Assortative-mating.Rmd-->

# Identity

## IBS

### Long runs of IBT

## IBD

## IBT


<!--chapter:end:25-Identity.Rmd-->

# Neutral theory of molecular evolution

## Nearly neutral theory of molecular evolution

<!--chapter:end:26-Neutral-theory.Rmd-->

