# LocusZoom-style regional association plots for candidate loci identified in the GWAS analyses.
# Association statistics are visualised together with pairwise linkage disequilibrium (LD; r²) relative to the lead SNP and local recombination rates. 
# The regional association plots are supplemented with Gene transcript annotations and Candidate cis-Regulatory Elements (cCREs). 

setwd("$PATH_TO_FINEMAPPING/LocusZoom")

library(dplyr) #v.1.2.0
library(janitor) #v.2.2.1
library(ggplot2) #v.4.0.2, for plotting 
library(patchwork) #v.1.3.2, arranging plots 
library(tidyr) #v.1.3.2
library(ggtranscript) #v1.0.0
library(grid) #v.4.5.2
library(scales)  #v.1.4.0, for plotting (fixing x axis) 


################################################################
####################### Prepp Locus Zoom #######################
################################################################
###### reading in the summary stats from GWAS ######
## The summary statistics have been subsetted into only containing markers from chromosomes of interest to save space and computational power. 
Pvals <- read.table("GWASoutput/3a1.ChrOfInterest_F.txt", header = T)

# If the file contains markers on more than one chromosome, the summary stats needs to be subsetted again, 
# Only keeping the chromosome of interest for investigated locus. 
#Ex
#chr <- subset(Pvals, Pvals$Chr=="10")
#Pvals <- chr


######reading in LD file ######
LD <- read.table("LD_Zoom/LD_3a_chr10_F.ld", header = T)


# Recombination maps can be found: 
# Kidd, J.M. Fine-scale recombination rates inferred using the canFam4 assembly are strongly correlated with previous maps of dog recombination. Mamm Genome 37, 12 (2026). https://doi.org/10.1007/s00335-025-10178-0
# https://doi.org/10.5281/zenodo.17095604

# Reading in the recombination map for specified chromosome. 
Recomb  <- read.table("RecombinationMAP/dog_average_canFam4_chr10_map.txt", header = TRUE)
#Changing headers: 
names(Recomb) <- c("POS","RateMb","RatecM")


## Merge dataframes: 
P_LD <- merge(Pvals,LD, by.x="Marker", by.y="SNP_B", all.x=TRUE)
P_LD_Recomb <- merge(P_LD,Recomb,by.x="BP",by.y="POS", all.x=TRUE)

# Keep only necessary columns, add -log10P: 
Keep <- select(P_LD_Recomb, "BP", "Pvalue", "R2", "RateMb")
Keep$NLogP <-(-log10(Keep$Pvalue))


# Choose Region Of Interest: 

# Actual region of interest, what do we want to plot: 
chr_use <- "10"
win_start <- 29330964
win_end <- 30330964

# Regions extending what we want to plot, to make sure the recombination extends to the edges of the plot: 
# Its advisable to choose this empirically, start with 0, do the plot below,
# add window until an additional recombination point can be seen on either side. 
win_start_extra <- (win_start - 500000)
win_end_extra <- (win_end+ 500000)


# Subset the data for just the extra chunk we want to extend outside of plot 
ToPlot <- Keep %>%   filter(BP >win_start_extra &  BP <win_end_extra )

#plot used for empirically deciding win_*_extra: 
plot(ToPlot$BP, ToPlot$RateMb)

## Make the bins for plotting R2
# Bin R2; make sure 1.0 is included
ToPlot$R2_cat <- cut(
  ToPlot$R2,
  breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, 1.01),
  labels = c("0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1"),
  right = FALSE,
  include.lowest = TRUE
)


ToPlot$R2_cat <- as.character(ToPlot$R2_cat)
ToPlot$R2_cat[is.na(ToPlot$R2)] <- "NA"


ToPlot$R2_cat <- factor(
  ToPlot$R2_cat,
  levels = c("0.8-1", "0.6-0.8", "0.4-0.6", "0.2-0.4", "0-0.2", "NA")
)


dummy_row <- ToPlot[1, ] %>%  # take structure from first row
  mutate(
    BP = NA_real_,
    P = NA_real_, 
    R2 = NA_real_,
    R2_cat = factor("0.6-0.8", levels = levels(ToPlot$R2_cat))
  )


ToPlot_WD <- bind_rows(ToPlot, dummy_row)



## Scale the Recombination Rate for plottin in same range
scale_factor <- max(ToPlot$NLogP, na.rm = TRUE) / max(ToPlot$RateMb, na.rm = TRUE)


# For plotting Recombination Rate
ToPlot.na <-na.omit(ToPlot)


###############################################################
####################### Plot Locus Zoom #######################
###############################################################


zoom <- ggplot(ToPlot, aes(x = BP, y = NLogP, colour = R2_cat)) +
  geom_hline(yintercept = 5, colour = "grey", linetype = 2) +  # line at -log10(P) = 5
  geom_hline(yintercept = 6.39794, colour = "red", linetype = 2) +  # line at -log10(P) = SIG
  geom_point(alpha = 0.8, size = 2) +
  geom_point(
    data = ToPlot[which.max(-log10(ToPlot$P)), ],
    aes(x = BP, y = NLogP),
    shape = 18, colour = "purple", size = 5, inherit.aes = FALSE) + # colouring the most significant SNP purple and changing shape. 
  geom_point(data = ToPlot.na, aes(x = BP, y = RateMb * scale_factor),color="grey", alpha=0.7) +
  geom_line(data = ToPlot.na, aes(x = BP, y = RateMb * scale_factor),color="grey", alpha=0.7) +
  scale_colour_manual(
    values = c(
      "0-0.2" = "darkblue",
      "0.2-0.4" = "lightblue",
      "0.4-0.6" = "green",
      "0.6-0.8" = "orange",
      "0.8-1" = "red",
      "NA" = "grey"
    ),
    name = expression("LD (" ~ r^2 ~ ")"),
    drop = FALSE    # keep all levels in the legend, even if some are empty
  ) +
  coord_cartesian(xlim = c(win_start, win_end),ylim = c(0,10), expand = FALSE) +
  labs(x = "", y = expression(-log[10](P))) +
  scale_y_continuous(
    name = expression(-log[10](p)),
    sec.axis = sec_axis(~ . / scale_factor, name = "Recomb. Rate\n(cM/Mb)"),
    limits = c(0, 10)) +
  scale_x_continuous(name="Mb", labels = label_number(scale = 1e-6))+
  theme_minimal() +
  theme(
    legend.title = element_text(),
    legend.position = "right"
  )+
  labs(title="3a1 chromosome 10 Female")

######################################################################
####################### Prepp Gene Transcripts #######################
######################################################################


# Gene models for CanFam4 could be accessed from: 
# https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9615/106

# The generation of Gene and Exon files were produced by ChengCheng Song. 

## Reding in files: 
#Autosomal genes 
genes <- read.table("canFam4.NCBI.genes.data.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
exons <- read.table("canFam4.NCBI.exon.data.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#Chr X genes
#genes <- read.table("canFam4.NCBI.genes.data_X.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#exons <- read.table("canFam4.NCBI.exon.data_X.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)


# Changing names: 
colnames(genes) <- c("Gene", "Chrom", "Start", "End", "Coding", "Strand")
colnames(exons) <- c("Gene", "ExonStart", "ExonEnd")

#Subset genes for the region if interest: 
genes_sub <- genes %>%
  filter(Chrom == chr_use, End >= win_start, Start <= win_end) %>%
  mutate(
    y = factor(Gene, levels = rev(unique(Gene))),
    x1 = ifelse(Strand == "-", End, Start),
    x2 = ifelse(Strand == "-", Start, End)
  )

#Subset exons for the region if interest: 
exons_sub <- exons %>%
  inner_join(genes_sub %>% select(Gene, Chrom, Strand, Coding, y), by = "Gene") %>%
  transmute(
    Gene, Chrom, Strand, Coding, y,
    start = ExonStart,
    end = ExonEnd
  ) %>%
  filter(end >= win_start, start <= win_end)



#####################################################################
####################### Plot Gene Transcripts #######################
#####################################################################
gene_plot <- ggplot() +
  geom_range(
    data = exons_sub,
    aes(xstart = start, xend = end, y = y, fill = "Exons"),
    height = 0.25,
    colour = "black",
    linewidth = 0.2
  ) +
  geom_segment(
    data = genes_sub,
    aes(x = x1, xend = x2, y = y, yend = y),
    colour = "black",
    linewidth = 0.5,
    arrow = arrow(length = unit(0.16, "cm"), ends = "last", type = "closed")
  ) +
  scale_fill_manual(
    values = c(Exons = "black"),
    name = "Coding"
  ) +
  guides(
    fill = guide_legend(override.aes = list(colour = "black", fill = "black"))
  ) +
  coord_cartesian(xlim = c(win_start, win_end), expand = FALSE) +
  scale_y_discrete(drop = FALSE) +
  scale_x_continuous(name="Mb", labels = label_number(scale = 1e-6))+
  theme_classic() +
  theme(
    axis.line.x = element_line(colour = "grey80", linewidth = 0.5),
    axis.line.y = element_line(colour = "grey80", linewidth = 0.5),
    axis.ticks = element_line(colour = "grey80"),
    axis.text.x = element_text(colour = "grey40"),
    axis.text.y = element_text(colour = "grey40"),
    axis.title.y = element_text(colour = "black"),
    axis.title.x = element_text(colour = "black"),
    legend.position = "none") +
  labs(x = "Base pair", y="Genes")



##########################################################
####################### Prepp cCRE #######################
##########################################################

# UU Brain cCRE unpublished,
#EpicDog elements located on CanFam3.1 (Keun Hong Son et al. ,Integrative mapping of the dog epigenome: Reference annotation for comparative intertissue and cross-species studies.Sci. Adv.9,eade3399(2023).DOI:10.1126/sciadv.ade3399)
#EpicDog elements converted from CanFam3.1 to CanFam4 co-ordinates using the UCSC liftover tool by Matthew Christmas, (https://genome.ucsc.edu/cgi-bin/hgLiftOver). 

# Reading in the data 
Epic_brain_cCRE <- read.table("EPIC_active_merged_brain.bed", header=FALSE, sep = "\t")
Epic_other_cCRE <- read.table("EPIC_active_merged_other.bed", header=FALSE, sep = "\t")
UU_brain_cCRE <- read.table("UUbrain_merged.bed", header=FALSE, sep = "\t")


# Set numbers per dataset so that can plot in that order
# if want addtional datasets, add them here
Epic_other_cCRE$dataset <- 1
Epic_brain_cCRE$dataset <- 2
UU_brain_cCRE$dataset <- 3

# bind the cCRE datasets
dat <- bind_rows(Epic_brain_cCRE, Epic_other_cCRE, UU_brain_cCRE)

# choose chr of interest! Note! need to add a chr here
dat_chr_use <- dat[dat$V1 == "chr10",]

##########################################################
####################### Plot cCRE ########################
##########################################################
# plot for dataset colour, region and order 
# note, there are no headers so will plot start-stop as V2-V3, V1 is chr

ccre_plot <- ggplot(dat_chr_use) +
  geom_segment(
    aes(x = V2, xend = V3, y = dataset, yend = dataset, colour = factor(dataset)),
    linewidth = 4
  ) +
  scale_colour_manual(values = c("1" = "olivedrab", "2" = "gold", "3" = "goldenrod")) +
  scale_y_continuous(
    breaks = c(1, 2, 3),
    labels = c("EPIC other", "EPIC brain", "UU brain"),
    limits = c(0.8, 3.2),
    expand = c(0, 0)
  ) +
  coord_cartesian(xlim = c(win_start, win_end), expand = FALSE) +
  scale_x_continuous(labels = label_number(scale = 1e-6))+
  theme_minimal() +
  theme(
    panel.grid.major.y = element_line(colour = "grey85"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none") +
  labs(x = paste0("Chr", chr_use, " Mb"), y = "cCREs")




########################################################################
####################### Addig all plots together #######################
########################################################################

#removing legend from ccre and gene plot 
gene_plot <- gene_plot + theme(legend.position = "none")
ccre_plot <- ccre_plot + theme(legend.position = "none")

zoom / gene_plot / ccre_plot +  plot_layout(widths = c(1, 1, 1), heights = unit(c(5, 2, 1), c('cm', 'null')))


