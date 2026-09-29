# Plotting_GWAS_results
A collection of R scripts to aid in the visualisation of GWAS results. Software versions are cited within individual scripts.

## EffectSize
R script `GWASeffectSizeAndImplicatedGenes.R` will plot three panels (i) SNP-ID-Allele, (ii) EffectSize(ConfidenceInterval), (iii) Genes panel (ii) is coloured by GWAS phenotypic group/"Factors". Package versions are cited in the script. 
![EffectSizeExample](Images/EffectSize.pdf)

## Fine-mapping
R script `Multipanel_LocusZoom.R` will plot three panels (i) regional association plot overlayed with scaled recombination rate, where -log10(p-value) has been coloured for pair-wise LD, (ii) genes in region, (iii) cCREs. For this script, pair-wise LD (r²) relative to the lead SNP is supplied using the output of PLINK, e.g., (`--r2`, `--ld-snp chrN:BP`, `--ld-window-kb 100000`, `--ld-window 999999`, `--ld-window-r2 0`),for each GWAS locus (`--snp chrN:BP`, `--window number`). Package versions and any required files, e.g., for recombination rate etc are cited in the script.
![LocusZoomExample](Images/LocusZoomExample.pdf)

R script `PlotCategoricalDataDistributionAndKeySNPalelleFreq.R` will plot a bar chart style plot for the distribution of individuals across categorical groups and overlay this with allele frequency per category for a selected SNP. Inputs are taken from PLINK as described in the script. Package versions are cited in the script. 

R script `PlotContinuousDataDistributionAndKeySNPalelleFreq.R` will plot a distribution plot of individuals per decile for a continuous phenotype. This is overlayed with allele frequency per decile for a selected SNP. Inputs are taken from PLINK as described in the script. Package versions are cited in the script. 
![AlleleFrequencyEaxamples](Images/AlleleFreqExamples.pdf)

