#####################################################################################
## GWAS effect size plot 
# by Jennifer Meadows - September 2026

# Plot has three panels, 
# (i) SNP-ID-Allele, (ii) EffectSize(ConfidenceInterval), (iii) Genes
# panel (ii) is coloured by GWAS phenotypic group/"Factors", 

#####################################################################################
# Plot were generated with R version 4.5.2 (2025-10-31) -- "[Not] Part in a Rumble"

# Load required packages
library(ggplot) #v4.0.2
library(tidyverse) #v2.2.0
library(patchwork) #v1.3.2
library(grid) #v4.5.2

dat <- read.table("EffectSizes_PlusGenes.txt", header = TRUE, sep = "\t")

# "EffectSizes_PlusGenes.txt" contains 6 columns, with headers as indicated below. 
# "Factor" is the grouping you wish to colour by, and genes are those brought
# forward via fine mapping or other means.

# Factor	SNP	Beta	Beta_Upper	Beta_Lower	Genes
# 1	rs123-A	0.4	0.48	0.32	gene1
# 2	rs456-G	1.4	1.67	1.13	gene1,gene6
# 3	rs980-G	1.04	1.2	0.88	gene17
# 3	rs348-C	0.49	0.58	0.4	gene34


# Preserve input order
dat <- dat %>%
  mutate(row_id = row_number(),
         Factor = as.factor(Factor))

# Common y scale: reverse so the first input row is at the top
y_limits <- c(nrow(dat) + 0.5, 0.5)

#####################################################################################
# (i) Left panel - SNP-ID-Allele

p_left <- ggplot(dat, aes(y = row_id, x = 1, label = SNP)) +
  geom_text(
    hjust = 1,
    size = 3,
    family = "sans"
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = y_limits,
    expand = c(0, 0)
  ) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(0, "pt")
  )

#####################################################################################
# (ii) Centre panel - EffectSize(ConfidenceInterval)

p_center <- ggplot(dat, aes(y = row_id)) +
  geom_errorbar(
    aes(
      xmin = Beta_Lower,
      xmax = Beta_Upper,
      colour = Factor
    ),
    orientation = "y",
    width = 0.3,
    linewidth = 0.7
  ) +
  geom_point(
    aes(x = Beta, colour = Factor),
    size = 2.5
  ) +
  geom_vline(
    xintercept = c(-0.5, 0.5, 1.0, 1.5),
    colour = "grey85",
    linewidth = 0.5
  ) +
  geom_vline(
    xintercept = 0,
    colour = "grey55",
    linewidth = 0.5,
    linetype = "dotted",
  ) +
  scale_colour_manual(
    values = c("1" = "#993404",
        "2" = "#FE9929",
        "3" = "#FEE391"),
    guide = "none"
  ) +
  scale_y_continuous(
    limits = y_limits,
    expand = c(0, 0)
  ) +
  labs(
    x = "Beta",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(0, "pt"),
    legend.position = "right"
  )

#####################################################################################
# (iii) Right panel - genes

# x is deliberately placed at zero, but the panel is expanded using the 
# actual maximum label width below.

p_right <- ggplot(dat, aes(y = row_id, x = 0, label = Genes)) +
  geom_text(
    hjust = 0,
    size = 3,
    family = "sans"
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = y_limits,
    expand = c(0, 0)
  ) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(0, "pt")
  )

#####################################################################################
# Combine panels

p_combined <- p_left + p_center + p_right +
  plot_layout(
    widths = c(2.0, 3.5, 7.5),
    guides = "keep"
  )

# plot and save in desired format
p_combined


