#####################################################################################
# Summarising Allele Frequency Data 
# CONTINUOUS PHENOTYPE

# by Jennifer Meadows - April 2026

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

# Here split the phenotype into the desired number of bins. Here it is 10.
df <- df %>%
  arrange(Pheno) %>%
  mutate(decile_group = ntile(Pheno, 10))


#Count the genotypes per decile 
# For a given SNP count alleles from columns a and b. 
# Note! Need to update to true alleles. In the below case, this SNP has alleles "G" and "T"

df %>%
  group_by(decile_group) %>%
  summarise(
    a_1 = sum(a == "G", na.rm = TRUE),
    b_1 = sum(b == "G", na.rm = TRUE),
    total_a = a_1 + b_1,
    .groups = "drop"
  )

df %>%
  group_by(decile_group) %>%
  summarise(
    a_2 = sum(a == "T", na.rm = TRUE),
    b_2 = sum(b == "T", na.rm = TRUE),
    total_2 = a_2 + b_2,
    .groups = "drop"
  )

# Use this information to make and allele frequency table
# In the example below, "DecileAlleleFreq.txt" has 2 columns, with headers "decile_group"	and "Freq"

#decile_group	Freq
# 1	0.995884774
# 2	0.991769547
#	3 0.981481481
# 4	0.981481481
# 5	0.985596708
# 6	0.971193416
# 7	0.967078189
# 8	0.965020576
# 9	0.962962963
# 10 0.94214876

#####################################################################################
## 2 PLOTTING
# Read in the allele frequencies
df_freq <- read.table("DecileAlleleFreq.txt", header = TRUE, sep = "\t")


# Define your ranges
p_min <- min(df$Pheno, na.rm = TRUE)
p_max <- max(df$Pheno, na.rm = TRUE)
f_min <- min(df_freq$Freq, na.rm = TRUE)
f_max <- max(df_freq$Freq, na.rm = TRUE)

# And then make the plot - here is the example is for GWAS phenotype "CCDF1" and SNP "chr15:26540641" Alele "G"
p5 <- ggplot() +
  geom_point(
    data = df, 
    aes(x = decile_group, y = Pheno), 
    color = "cornflowerblue", 
    position = position_jitter(width = 0.2, height = 0)
  ) +
  geom_line(
    data = df_freq, 
    aes(
      x = decile_group, 
      y = p_min + (Freq - f_min) * (p_max - p_min) / (f_max - f_min)
    ), 
    color = "black", 
    size = 1
  ) +
  scale_x_continuous(name = "Decile",breaks = 1:10) +
  scale_y_continuous(
    name = "Phenotype",
    sec.axis = sec_axis(
      ~ ( . - p_min ) * (f_max - f_min) / (p_max - p_min) + f_min,
      name = "Frequency (G)"
    )
  ) +
  ggtitle("CCDF1 - chr15:26540641") +
  theme_minimal()


#####################################################################################
## 3. If have multiplots, including some generated with "PlotCategoricalDataDistributionAndKeySNPalelleFreq.R"
# combine with patchwork and then save as desired

## e.g., 
# pAll <- (p + p2) +
#  plot_layout(ncol = 2)


# ggsave("/path/Figure_AlleleFreq.pdf", plot = pAll,
#       width = 400, height = 100, units = "mm")

# ggsave("/path/Figure_AlleleFreq.svg", plot = pAll,
#       width = 400, height = 100, units = "mm")