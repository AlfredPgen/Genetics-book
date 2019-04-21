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

## DNA processing quality control

## Batch effects

https://www.bioconductor.org/packages/devel/bioc/vignettes/GWASTools/inst/doc/DataCleaning.pdf

The overall goal of this step is to check the quality of the sample batches. Substantial quality control is done by the genotyping centers prior to releasing the genotype data. However, it is possible that quality control for batches is still lower than desired. If a lower quality batch is detected then it may be necessary to re-run the genotyping for that batch. We can check the batch quality by comparing the missing call rates between batches and looking for significant allele frequency differences between batches.

### Calculation of missing call rate for samples and SNPs

The ﬁrst step is to calculate the missing call rates for each SNP and for each sample. A high missing call rate for a sample is often indicative of a poorly performing sample. It has been seen that samples from DNA that has undergone whole-genome ampliﬁcation (WGA) have a relatively higher missing call rate. Similarly a high missing call rate for a SNP is indicative of a problem SNP. Experience from the GENEVA studies has shown that there seem to be a subset of SNPs from which genotype calls are more diﬃcult to make than others. We calculate the missing call rates in a two step process: ﬁrst the missing call rates over all samples and SNPs are calculated, then the missing call rates are calculated again, ﬁltering out SNPs and samples that have an initial missing call rate greater than 0.05. The initial SNP missing call rate over all samples is saved in the SNP annotation data ﬁle as missing.n1. The analogous idea is applied to the samples: missing.e1 is saved in the sample annotation ﬁle and corresponds to the missing call rate per sample over all SNPs, excluding those SNPs with all calls missing. The missing.n2 is calculated as the call rate per SNP over all samples whose missing.e1 is less than 0.05. Again, similarly for the samples, missing.e2 is calculated for each sample over all SNPs with missing.n2 values less than 0.05. It is important to remember that the Y chromosome values should be calculated for males only, since we expect females to have no genotype values for the Y chromosome, although an occasional probe on the Y chromosome is called in a female.
If any samples have a high missing rate, we recommend further investigation of what may be causing the missing calls; the samples with a missing call rate greater than 0.05 should be ﬁltered out due to low sample quality.

### Calculation of missing call rates by batch

The missing call rate by batch is calculated to check that there are no batches with comparatively lower call rates. Usually a“batch”is a plate containing samples that were processed together through the genotyping chemistry. In this case all samples were run on diﬀerent plates (as controls for another dataset).

### Testing for allele frequency differences in batches

In this step, the chi-square test for diﬀerences in allelic frequency is performed between each batch individually and a pool of all the other batches in the study. We then look at the mean χ<sup>2</sup> statistic over all SNPs for each batch as a function of the ethnic composition of samples in a batch.
Next we test for association between batches and population groups, using a χ<sup>2</sup> contingency test. Then we look at the relationship between the ethnic composition of each batch and the previously calculated χ<sup>2</sup> test of allelic frequency between each batch and a pool of the other batches. The point is to look for batches that diﬀer from others of similar ethnic composition, which might indicate a batch eﬀect due to genotyping artifact. In this experiment, there are only a few batches and wide variations in race among batches, so it is diﬃcult to interpret the results. In larger GWAS experiments, we generally observe a U-shaped curve of allelic frequency test statistic as a function of ethnic composition.
The χ<sup>2</sup> test is not suitable when the 2×2 tables for each SNP have very small values. For arrays in which many SNPs have very low minor allele frequency, Fisher’s exact test is more appropriate. 

## Sample quality control

### Cryptic relatedness

### Population stratification

Sometimes finding an association can be confounded by population stratification. This is because a condition may be more prevalent in one group of people than in a different group, resulting in a spurious association between the condition or trait being tested for and any genetic characteristics which vary between the two different groups of people.

While it is good practice for studies to be based on as homogeneous a group of test subjects as possible, it has been noted in [Price, 2006] that even the mild variation in genetic characteristics among those who classify themselves as belonging to one ethnic group or another can be problematic enough to confound a study done over thousands of genetic markers.

Hidden population stratification may be thought of as a non-zero F<sub>st</sub> between unknown groupings of samples.

### Heterozygosity and missingness outliers

### Differential missingness

### Sex chromosome anomalies

## Marker quality control

### Genotyping concordance

In genotyping studies where DNA is directly assayed for positions of variance, concordance is a measure of the percentage of SNPs that are measured as identical. Samples from the same individual or identical twins theoretically have a concordance of 100%, but due to assaying errors and somatic mutations, they are usually found in the range of 99% to 99.95%. Concordance can therefore be used as a method of assessing the accuracy of a genotyping assay platform.

### Mendelian errors

### Genotype call rate

### Minor allele frequency

### Hardy-Weinberg equilibrium outliers

### Additional QC for regions like MHC

### Ambigious nucleotides

If the base and target data were generated using different genotyping chips and the chromosome strand (+/-) for either is unknown, then it is not possible to match ambiguous SNPs (i.e. those with complementary alleles, either C/G or A/T) across the data sets, because it will be unknown whether the base and target data are referring to the same allele or not. While allele frequencies can be used to infer which alleles match, it is recommended to remove all ambiguous SNPs since the allele frequencies provided in base GWAS are often those from resources such as the 1000G project, and so aligning alleles according to their frequency could lead to systematic biases. 

### Non-matching nucleotides

When there is a non-ambiguous mismatch in allele coding between the data sets, such as A/C in the base and G/T in the target data, then this can be resolved by ‘flipping’ the alleles in the target data to their complementary alleles. 

### Quality control prior to meta-analysis

Allele Frequency Plots (AF Plots): looking for errors in allele frequencies and strand orientations by visually inspecting a plot of the sample allele frequency of filtered SNPs against the frequency in the 1000 Genomes phase 1 version 3 European panel3 for example.
P value vs Z-statistic Plots (PZ Plots): looking for the consistency between the reported P values and the P values implied by the coefficient estimates and standard errors in
individual cohort.
Quantile-Quantile Plots (QQ Plots): looking for the cohort-level QQ plots to look for evidence of unaccounted-for stratification.
Predicted vs Reported Standard-Error Plots (PRS Plots): maaking sure that the standard errors reported in individuals cohorts are approximately consistent with the reported sample size, allele frequency, and phenotype distribution. 
Use of bivariate LD score regression to verify that the estimated genetic correlations between all large cohorts (defined as N > 10,000) are large and positive.

## X-chromosome quality control

The X chromosome plays an important role in complex human traits and diseases, especially those with sexually dimorphic characteristics. Special attention needs to be given to the analysis of X due to its unique inheritance pattern and X-inactivation. These analytical complications have resulted in exclusion or mishandling of the X chromosome in the majority of genome-wide association studies (GWAS) to date.

## Single marker regression

Summary statistics can be obtained using one of the following tests: Correlation/Trend Test, Armitage Trend Test, Exact Form of Armitage Test, (Pearson) Chi-Squared Test, (Pearson) Chi-Squared Test with Yates’ Correction, Fisher’s Exact Test, Odds Ratio with Confidence Limits, Analysis of Deviance (e.g. different variance heterogeneity tests), F-Test, Logistic Regression, Linear Regression.

### Allelic test

### Genotypic test

### Additive model

Here, the genotype is coded in terms of the number of specific allele at a given locus.

### Dominant model

Here, the genotype with at least 1 copy of a specific allele at a given locus is coded as 1 and other genotypes as 0.

### Recessive model

Here, the genotype with at least 2 copies of a specific allele at a given locus is coded as 1 and other genotypes as 0.

### Categorical phenotype

### Multi-allelic GWAS

## Two-stage approach

## Haplotype GWAs design

## Joint analysis (all independent markers)

### Genomic control

In an ordinary GWAS, genomic control (GC) is used to shrink any existing inflation of the test scores (-log10 *p*-values). When testing for the single genetic effect in the GWAS, the null distribution of the test statistic for the nominal p-values is χ<sup>2</sup> with 1 degree of freedom. Since most of the SNPs are not expected to be associated with the trait, the sample distribution of the chi-squares across the genome should resemble the null distribution. If there is inflation, the chi-squares are adjusted using λ, i.e. the inflation factor estimated by comparing the distribution of the sample χ<sup>2</sup>’s and χ<sup>2</sup> distribution with 1 degree of freedom.

#### λ<sub>1000<sub>

Since λ scales with sample size, some have found it informative to report λ<sub>1000</sub>. This is equivalent to a study of 1000 cases and 1000 controls and can be calculated by rescaling λ with 1 + (λ - 1) x (1/case + 1/control) x 500, where case and control refers to the number of cases and controls respectively.

##### Stratified λ<sub>GC<sub>

Because the strength of λ<sub>GC</sub> deviation depends on allele frequency and very large sample sizes become available, it could be useful to report λ<sub>GC</sub> values based on certain MAF bins.

## Multimarker single gene-based approaches

https://cran.r-project.org/web/packages/aSPU/aSPU.pdf

## VEGAS

## Multimarker gene-set approaches (a.k.a. pathway analysis)

## fastBAT

## MAGMA

## VEGA

## Extensions to binary and categorical phenotypes

### Threshold model

## Analysis of rare variants

Check Lee et al (2014) for a review

Due to the low frequencies of rare variants, classical single marker tests commonly used in genome-wide association studies (GWAS) for studying common variants effects are not applicable.
In view of the lack of power of single marker analysis of rare variants, methods investigating rare variation are typically region-based tests where one tests for the cumulative effects of the rare variants in a region. These region-based methods can be broadly classified into three classes: burden tests, non-burden tests and hybrid of the two. The key difference between burden and non-burden tests is how the cumulative effects of the rare variants are combined for association testing. For the commonly used simple burden tests, one summarizes the rare variants within a region as a single summary genetic burden variable, e.g. the total number of rare variants in a region, and tests its association with a trait. Burden tests implicitly assume all the rare variants in the region under consideration are causal and are associated with the phenotype in the same direction and magnitude. Hence, they all share the limitation of substantial power loss when there are many non-causal genetic variants in a region and/or when there are both protective and harmful variants.
Several region-based non-burden tests have been proposed by aggregating marginal test statistics (Neale et al., 2011; Basu and Pan, 2011; Lin and Tang, 2011). One such test is the sequence kernel association test (SKAT) (Wu et al., 2011), where one summarizes the rare variants in the region using a kernel function, and then test for association with the trait of interest using a variance component score test. SKAT is robust to the signs and magnitudes of the associations of rare variants with a trait. It is more powerful than the burden tests when the effects are in different directions or the majority of variants in a region are null, but is less powerful than burden tests when most variants in a region are causal and the effects are in the same direction. Several hybrids of the two methods have been proposed to improve test power and robustness (Lee et al., 2012; Derkach et al., 2013; Sun et al., 2013).

### Collapsing methods based on pooling multiple rare variants (burden or adaptive burden tests)

#### Sum test

The most powerful multi-marker test when there are no causal variants with effects in opposite directions and when there are few or no non-causal RVs. Otherwise, it suffers from substantial loss of power.

#### Cohort Allelic Sums test (CAST)

#### Combined Multivariate Collapsing (CMC)

#### Weighted Sum test (WSS)

#### Kernel Based Adaptive Cluster (KBAC)

#### Replication Based Test (RBT)

#### ARIEL test

#### The EREC test

### Methods treating rare variant effects as random (Variance-component tests)

#### The SSU approach

Has good power in the presence of opposite association directions and small fraction of causal RVs.

#### C-alpha test

#### SKAT

Has good power in the presence of opposite association directions and non-causal RVs.
It was recently suggested that using SKAT in the presence of RVs and common variants (CVs) may be less optimal because RVs are weighted to have much more importance than CVs (Ionita-Laza et al., 2013). 

### Methods based on model selection

The model-selection approaches perform in the middle of random eﬀect and collapsing methods. One issue common to model-selection methods is that  model selection approaches use dimension-reduction strategies to substantially reduce the number of parameters one would require to ﬁt these large number of RVs. Hence, any model we can construct will never be the true model that generated the data we observe. In other words, the set of models is clearly misspeciﬁed, and model selection is best seen as a way of approximating, rather than identifying, full reality (Burnham and Anderson (2002), pp. 20-23).

#### Seq-aSum

#### Seq-aSum-VS

The Seq-aSum-VS approach classiﬁes RVs based on the direction of association (‘+1’ for positive association, ‘-1’ for negative association and ‘0’ for no association) and implements a sequential variable selection scheme to select the best model for association between the SNP-set and the disease. The only diﬀerence between the Seq-aSum approach and the Seq-aSum-VS approach is that the variable selection (‘0’ allocation for a variant) is not implemented in the former. The Seq-aSum-VS approach starts with putting all the RVs in the ‘+1’ group and proceeds by moving each RV sequentially to the other two groups and ﬁnally chooses the allocation (‘+1’,‘-1’, or ‘0’ ) with highest likelihood to the RV. The process of choosing the best model in Basu and Pan (2011)’s method can be compared to a stepwise regression, where one may not always ﬁnd the best model due to this selection scheme. This is especially true if a particular allocation results in a slightly higher likelihood than the other two allocations. In this case, choosing the allocation with highest likelihood for a SNP might not be optimal, rather it might be more eﬃcient to allow multiple allocations for a RV and construct a test that takes into account multiple plausible models for the disease-RV association. Moreover, the performance of the sequential search often depends on the ordering of the variants in this search mechanism. A model-averaging approach could potentially reduce the dependency on the ordering of the variants in this sequential search.

#### Variable Threshold Test (VT)

#### RARECOVER

#### Selective grouping method

#### Step-Up

### Combination of collapsing and random effects methods

According to Basu and Pan (2011), the model selection methods, especially Seq-aSum-VS approach, performed very well when there were both protective and deleterious causal RVs and very few non-causal RVs, but the performance of the Seq-aSum-VS approach was not very impressive in the presence of a moderate or large number of non-causal RVs. These and other ﬁndings (Basu and Pan, 2011) have led to combining the strengths of collapsing and random eﬀect methods.

#### SKAT-O

#### SKAT-C

#### Fisher method

#### MiST

### EC test

Exponentially combines score statistics. Powerful when a very small proportion of variants are causal.

### Family-based tests

https://www.omicsonline.org/open-access/literature-reviews-on-methods-for-rare-variant-association-studies-2161-0436-1000133.pdf

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

Can only be applied to limited numbers of samples
Most complete ascertainment of variation

### Low coverage whole genome sequencing

Can be applied to moderate numbers of samples
Very complete ascertainment of shared variation
Less complete ascertainment of rare variants

## Analysis of multiple traits

## Mixed-model association analysis

GWAS using mixed models is appealing for several reasons:
1) More powerful in very large GWAS
2) Reduces the need for sample exclusion
3) Amplifies effective sample sizes via conditioning on polygenic predictions from genome-wide SNPs.

When using mixed models for association analysis, care must be taken to consider non-additive effects when retaining related individuals.

### EMMAX

### Fast-LMM

### GEMMA

### BOLT-LMM

When BOLT-LMM was used on very large sample sizes, analyses revealed subtleties in the interpretation of LD score regression intercepts as a means of differentiating polygenicity from confounding; the attenuation ratio was proposed to be possibly a more suitable metric as sample sizes increase.

### Caveats

First, chi-squared-based tests (such as BOLT-LMM) can incur inflated type I error rates when used to analyze highly unbalanced case–control traits (case fractions <10%). There are two solutions for this.
1) Increase MAF.
2) Use saddlepoint approximation (SAIGE software).

Second, conditioning on genome-wide signal can produce loss of power under case–control ascertainment.

## Penalized regression GWAS

## Bayesian GWAS

## Machine learning for GWAS

## Expected increase in GWAS loci with sample size

## The joint effect of genotypes over all traits

## Adjustment for winner's curse

One type of adjustment for winer's curse (using empirical Bayes), depends on the assumption that SNP effects are drawn randomly from the following mixture distribution 

## Adjustment for assortative mating

## Adjustment for attenuation bias

## Enrichment in candidate genes

Atwell et al. (2010) introduced a method for evaluating the enrichment of strong, but not necessarily genome-wide significant, signals for SNPs in candidate genes. An enrichment of such signals indicates that the analysis identifies true signals rather than random noise.

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

Placeholder


## Single step methods
### Case-control
### Prospective likelihood-based approach
### Retrospective likelihood approach
#### Multiplicative scale
#### Additive scale
### Case-only
### Empirical Bayes and Bayesian Model Averaging
### Other Bayesian approaches
#### Interactions using haplotypes
## GxE in the context of family studies
## Multi stage methods
## Joint tests
### Gene-based
## Set-based interaction tests
### Combining multiple environmental factors
### Multi-trait multi-GxE tests
### Variance heterogeneity
### Levene's test
### Brown-Forsythe test
### Bartlett's test
### Bartlett's test with prior rank transformation to normality
### Generalized Levene's scale tests
### The Lepage test
### The *D*-test
### Regression using the squared Z-score
### Gamma regression models
### Two-step screening on residual variance heterogeneity
#### VH
#### YGVH
## Mean-variance QTL (joint test)
### Conditional quantile regression 
### Sliced inverse regression
### Semiparametric model for vQTL mapping
### Variance heterogeneity for related individuals
## RELIEF and other machine learning tools
### Multidimensionality reduction
## Meta-analytic GxE approaches
## Gene-environment correlation
## Other challenges
### Controlling for covariate interactions
## Candidate gene-by-environment interaction studies
## Higher order interactions

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

