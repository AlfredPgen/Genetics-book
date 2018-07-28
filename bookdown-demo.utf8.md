--- 
title: "A handbook for Computational Genetics"
author: "Alfred Pozarickij"
date: "2018-07-28"
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

The scope of this book is to provide an outline of computational methods available for the analysis of genetic data.

First chapter introduces methods to infer population parameters. 

Next X chapters focus on population genetics.

In this book, the amount of mathematics and statistics is kept to a minimum. Only methods designed specifically to address issues in genetics are shown. I have produced another book [insert link here], which is intended to familiarise the reader with commonly used approaches and develop some intuition behind them. By no means the list is comprehensive and only serves as a quick guide. The internet provides much more information regarding this topic.

Rather than providing references at the end of each chapter, I decided to combine them into supplementary text [insert link here]. References are arranged according to different topics discussed in this book. I tried to do my best to cause as little confusion as possible.

Finally, don't hesitate to contact me if I have not included your favorite method (apozarickij@gmail.com). I would be more than happy to hear about it.








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

<!--chapter:end:01-Descriptives.Rmd-->

# Sequencing technologies

<!--chapter:end:02-Sequencing-technologies.Rmd-->

# Genome-wide association analysis

## Genotype calling algorithms

## DNA processing quality control

## Sample quality control

### Cryptic relatedness

### Population stratification

### Heterozygosity and missingness outliers

### Differential missingness

### Sex chromosome anomalies

## Marker quality control

### Genotyping concordance

### Switch rate

### Genotype call rate

### Minor allele frequency

### Hardy-Weinberg equilibrium outliers

### Additional QC for regions like MHC

### Ambigious nucleotides

### Non-matching nucleotides

## X-chromosome quality control

## Single marker regression

### Trend test

### Alleles test

## Two-stage approach

## Haplotype GWAs design

### Genomic control

## Multimarker gene-set based approaches

https://cran.r-project.org/web/packages/aSPU/aSPU.pdf

## Extensions to binary and categorical phenotypes

### Threshold model

## Analysis of rare variants

Check Lee et al (2014) for a review

Due to the low frequencies of rare variants, classical single marker tests commonly used in genome-wide association studies (GWAS) for studying common variants effects are not applicable.
In view of the lack of power of single marker analysis of rare variants, methods investigating rare variation are typically region-based tests where one tests for the cumulative effects of the rare variants in a region. These region-based methods can be broadly classified into three classes: burden tests, non-burden tests and hybrid of the two. The key difference between burden and non-burden tests is how the cumulative effects of the rare variants are combined for association testing. For the commonly used simple burden tests, one summarizes the rare variants within a region as a single summary genetic burden variable, e.g. the total number of rare variants in a region, and tests its association with a trait. Burden tests implicitly assume all the rare variants in the region under consideration are causal and are associated with the phenotype in the same direction and magnitude. Hence, they all share the limitation of substantial power loss when there are many non-causal genetic variants in a region and/or when there are both protective and harmful variants.
Several region-based non-burden tests have been proposed by aggregating marginal test statistics (Neale et al., 2011; Basu and Pan, 2011; Lin and Tang, 2011). One such test is the sequence kernel association test (SKAT) (Wu et al., 2011), where one summarizes the rare variants in the region using a kernel function, and then test for association with the trait of interest using a variance component score test. SKAT is robust to the signs and magnitudes of the associations of rare variants with a trait. It is more powerful than the burden tests when the effects are in different directions or the majority of variants in a region are null, but is less powerful than burden tests when most variants in a region are causal and the effects are in the same direction. Several hybrids of the two methods have been proposed to improve test power and robustness (Lee et al., 2012; Derkach et al., 2013; Sun et al., 2013).

### Collapsing methods based on pooling multiple rare variants

#### Sum test

The most powerful multi-marker test when there are no causal variants with effects in opposite directions and when there are few or no non-causal RVs. Otherwise, it suffers from substantial loss of power.

#### Cohort Allelic Sums test (CAST)

#### Combined Multivariate Collapsing (CMC)

#### Weighted Sum test (W-Sum)

#### Kernel Based Adaptive Cluster (KBAC)

#### Replication Based Test (RBT)

#### ARIEL test

#### The EREC test

### Methods treating rare variant effects as random

#### The SSU approach

Has good power in the presence of opposite association directions and non-causal RVs.

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

## Analysis of X, Y and mitochondrial chromosomes

## Analysis of copy number variants

### Common variation

### Analysis of rare variants

## Analysis of multi-ethnic samples

## Analysis of indirect genetic effects

## GWAS vs whole-genome association

## Analysis of multiple traits

## Mixed-model association analysis

## Penalized regression GWAS

## Bayesian GWAS

## Machine learning for GWAS


<!--chapter:end:03-Genome-wide-association-analyses.Rmd-->

# Heritability

## Realized heritability

### Evolvability

### Reliability

## Twin studies

## GCTA

## LD-score regression

## LDAK

<!--chapter:end:04-Heritability.Rmd-->

# Genomic prediction

## Polygenic risk scores

## Gene-based polygenic score (POLARIS)

## Pathway-based polygenic risk score

## LD adjusted PRS

### LDpred with functional annotation

## Annopred

## Pleiopred

## BLUP

## Bayesian Zoo

## Reproducing kernel Hilbert space

## Machine learning methods


<!--chapter:end:05-Genomic-prediction.Rmd-->

# Pleiotropy

## Direct

## Indirect


<!--chapter:end:06-Pleiotropy.Rmd-->

# Pathway-analysis

<!--chapter:end:07-Pathway-analysis.Rmd-->

# Functional annotation

<!--chapter:end:08-Functional-annotation.Rmd-->

# Causal inference

## Gene-knockout

## Conditioning

## Finemapping

## Mendelian Randomization

<!--chapter:end:09-Inferring-causality.Rmd-->

# Combining multiple datasets

## Meta-analysis

## Mega-analysis

<!--chapter:end:10-Analysis-of-multiple-datasets.Rmd-->

# Gene-environment interaction

Identification of gene-environment interactions has important implications for understanding underlying disease etiology and developing disease prevention and intervention strategies.

## Single step methods

### Case-control

### Case-only

### Empirical Bayes and Bayesian Model Averaging

## Multi stage methods

## Joint tests

## Set-based interaction tests

There are multiple reasons for using set-based gene-environment interaction tests. 

1. Multiple comparison adjustments for a large number of markers across the genome could result in power loss.
2. Closely located SNPs are correlated because of linkage disequilibrium. Multiple tests for GxE in these single-marker-based GxE models are even more dependent, as interaction terms in these models share the same environmental variable. Dependence among multiple tests can result in incorrect Type 1 error rates and causes bias in standard multiple comparison adjustments and this bias is often difficult to correct.
3. The single-marker GxE test does not interrogate the joint effects of multiple SNPs that have similar biological functions. When the main effects of multiple SNPs in a set are  associated with a disease/trait, the classical single marker regression interaction test can be biased.

Lin et al. (2013) developed a method to analyse GxE for a set of markers using generalized linear mixed models. The method tests for SNP-set by environment interactions using a variance component test, and because a set of variants will likely be correlated, the main SNP effect estimates under the null hypothesis are obtained using ridge regression. 
Their software is called GESAT. Here, they model GxE effects as random, as opposed to the classical approach of treating *BETAj*’s as fixed effects followed by a test with *p* degrees of freedom. The latter approach can suffer from power loss when *p* is moderate/large, and numerical difficulties when some genetic markers in the set are in high LD. 
The model allows to adjust for the main effects of all SNPs while simultaneously testing for the interactions between the SNPs in the region and environmental variable. For unbalanced designs when a binary environmental exposure has a low frequency in one category, GESAT is most advantageous over single marker regression GxE test. Such unbalanced designs can occur due to case–control sampling and the strong association of an environmental factor with disease. When the effect size is modest, GESAT performs better that single marker regression GxE test, but when the effect size is strong, the opposite is true. 
Their simulations suggest that the power of GESAT seems fairly robust to the dependence between G and E.

The same approach can be applied to investigating various other biological problems. For example, we can test for the interactions between gene expressions in a pathway or network and an environmental variable by simply replacing G by gene expressions in a gene-set.

Existing methods for assessing common variants by environment interactions such as Gene-Environment Set Association Test (GESAT) (Lin et al., 2013) have several limitations when applied for rare variants. GESAT estimates the main effects of the common variants by applying a L2 penalty on the genotypes scaled to unit variance; this assumes that the main effects of the scaled genotypes are comparable in magnitudes, which may not hold in the case of rare variants. GESAT also assumes that the regression coefficients of the rare variants by environment interactions are independent of each other, and suffers from power loss when most rare variants in a gene interact with the environmental factor and the interaction effects have the same direction.

Similar idea from GESAT was later extended and applied to analyse rare variants. GESAT have several limitations when applied for rare variants. GESAT estimates the main effects of the common variants by applying a L2 penalty on the genotypes scaled to unit variance; this assumes that the main effects of the scaled genotypes are comparable in magnitudes, which may not hold in the case of rare variants. GESAT also assumes that the regression coefficients of the rare variants by environment interactions are independent of each other, and suffers from power loss when most rare variants in a gene interact with the environmental factor and the interaction effects have the same direction. The new proposed test iSKAT is optimal in a class of variance component tests and is powerful and robust to the proportion of variants in a gene that interact with environment and the signs of the effects. This test properly controls for the main effects of the rare variants using weighted ridge regression while adjusting for covariates.

A naive approach to assess rare variants by environment interactions is to extend the burden test by fitting a model with both the summary genetic burden variable, environment, and their interaction, and performing a one degree of freedom test for the interaction. However, when there are multiple causal variants with their main effects having different magnitudes and/or signs, such a burden rare variant by environment test fails, and may lead to inflated Type 1 error rates. This is because adjusting for the main effects of the multiple causal variants using a single summary genetic burden variable is inappropriate. Likewise, a naive approach to assess rare variants by environment interactions using SKAT by including the main effects of rare variants as part of covariates and applying SKAT to the interaction terms is problematic. This is because SKAT only allows adjustment of a small number of covariates and cannot handle the presence of a large number of rare variants in a region. Furthermore since the rare variants are observed in low frequency, a model with all the rare variants as main effects will be highly unstable and may not even converge.

Both GESAT and iSKAT are able to incorporate multiple environmental variables.

No gene-based GxE test exists for the analysis of the sex chromosomes as of writing this text.

## Combining multiple environments

## Variance heterogeneity

### Levene's test

### Two-step screening on residual variance heterogeneity

## Conditional quantile regression 

## RELIEF and other machine learning tools

### Multidimensionality reduction

## Meta-analytic GxE approaches

<!--chapter:end:11-GxE.Rmd-->

# Gene-gene interaction

## Single step methods

## Multi stage methods

## Machine learning methods

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

# Additional points to consider

## Kinship matrix

### Path coefficients

## Genetic relationship matrix

## Animal models

## Phasing

## Haplotyping

## Statistical power

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

## Exended phenotype

## Genome sizes








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

