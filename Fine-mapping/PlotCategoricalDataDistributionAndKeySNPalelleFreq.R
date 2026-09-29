#####################################################################################
# Summarising Allele Frequency Data 
# CATEGORICAL PHENOTYPE; e.g., 1, 2, 3, 4, 5, 6

# by Jennifer Meadows - July 2026

#####################################################################################
# Plot were generated with R version 4.5.2 (2025-10-31) -- "[Not] Part in a Rumble"

# Load required packages
library(dplyr) #v1.2.0
library(janitor) #v2.2.1 
library(ggplot2) #v4.0.2

#####################################################################################

## 1. Use PLINK ped file to extract allele counts for SNP(s) of interest. PLINK map file is for sanity checking
# A per SNP map and ped can be generated from PLINK using the "--snp" command, e.g., 
# plink --bfile FILENAME --dog --snp SNP-ID --recode --out plink  


df <- read.table("plink.ped", header=FALSE, sep = " ")
map.df <- read.table("plink.map", header=FALSE, sep = "\t")
# Note, different seps in above

# Add a header so can call column later
names(df) <- c("FID", "IID", "P", "M", "Sex", "Pheno", "a","b")


# For a given SNP count alleles from columns a and b. 
# Note! Need to update to true alleles. In the below case, this SNP has alleles "A" and "C"

df %>%
  group_by(Pheno) %>%
  summarise(
    a_1 = sum(a == "A", na.rm = TRUE),
    b_1 = sum(b == "A", na.rm = TRUE),
    total_1 = a_1 + b_1,
    .groups = "drop"
  )

df %>%
  group_by(Pheno) %>%
  summarise(
    a_2 = sum(a == "C", na.rm = TRUE),
    b_2 = sum(b == "C", na.rm = TRUE),
    total_2 = a_2 + b_2,
    .groups = "drop"
  )

#####################################################################################
## 2. Use above data to make table containing counts of dogs per category and frequency of effect allele

## Expecting format similar "AlleleFreq_Category.txt" below
# Category	CountDog	Freq
# 1	602	0.102990033
# 2	757	0.064729194
# 3	504	0.05952381
# 4	309	0.040453074
# 5	139	0.057553957
# 6	111	0.018018018

## Read in the allele frequencies

Summary_SNP <- read.table("AlleleFreq_Category.txt", header = TRUE, sep = "\t")

scale_factor <- max(Summary_SNP$CountDog) / max(Summary_SNP$Freq)

#####################################################################################
## 3. Plot
# Update plot axis and title names to match your SNP of interest.
# The plot below is for GWAS phenotype "Item 146", SNP name "chr1:75148477", and the effect alelle "A"

p1 <- ggplot(Summary_SNP, aes(x = Category)) +
  geom_col(aes(y = CountDog),
           fill = "cornflowerblue",
           alpha = 0.5) +
  geom_line(aes(y = Freq * scale_factor, group = 1),
            color = "black",
            linewidth = 1) +
  scale_y_continuous(
    name = "Dogs (n)",
    sec.axis = sec_axis(~ . / scale_factor, name = "Frequency (A)")
  ) +
  ggtitle("Item 146 - chr1:75148477") + 
  theme_minimal()

p1a <- ggplot(Summary_SNP, aes(x = Category)) +
  geom_col(aes(y = CountDog),
           fill = "cornflowerblue",
           alpha = 0.5) +
  geom_line(aes(y = Freq * scale_factor, group = 1),
            color = "black",
            linewidth = 1) +
  scale_y_continuous(
    name = "Dogs (n)",
    sec.axis = sec_axis(~ . / scale_factor, name = "Item 146 - chr1:75148477\nFrequency (A)")
  ) +
  theme_minimal()  +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10))

