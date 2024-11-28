## R code for TE analyses with Command Line and PDF Output
## Shujun Ou
## Date: 11/11/2023

# Define a function to print the help message.
print_help <- function() {
  cat("
Usage: Rscript div_plot.R [data_file] [species_accession]

[data_file]          : Name of the file produced by div_table2.pl (e.g., 'species.div_long')
[species_accession]  : Species & accession for graph labeling (e.g., 'Maize_B73')

This script generates a PDF depicting the sequence divergence of TEs relative to the repeat lib. 

R and its packages can be installed with conda.
e.g., conda config --add channels conda-forge && conda config --add channels CRAN && conda create -n test r-essentials r-tidyr r-dplyr r-ggplot2 r-scales r-here

")
}

# Accept command line arguments.
args <- commandArgs(trailingOnly = TRUE)

# Check if the user has requested help or has not provided the correct number of arguments.
if (length(args) != 2 || any(args %in% c('-h', '--help'))) {
  print_help()
  quit(save = "no", status = 0)  # Exit the script after displaying help message
}

# Load necessary libraries.
library(tidyr)
library(dplyr)
library(ggplot2)
library(scales)
#library(here)

# Set the working directory to the script's location
script_dir <- dirname(Sys.getenv("R_SCRIPT"))
#setwd(here())

data_file <- args[1]  # Data file name (e.g., 'testmaize.txt')
species_accession <- args[2]  # Species accession label (e.g., 'Maize_B73')

#ADDED
# Define the color mapping for TE categories
TE_colors2 <- c("DNA/Helitron"="#550000",
                "DNA/transposon"="#AA0000",
                "LTR/BARE-2"="#D40000",
                "LTR/Copia"="#FF0000",
                "LTR/Gypsy"="#FF2F52",
                "LTR/LARD"="#FF5500",
                "LTR/TR_GAG"="#FF8000",
                "LTR/TRIM"="#FFAA00",
                "LTR/unknown"="#FFD400",
                "LINE/L1"="#FFFF00",
                "LINE/RTE"="#D4FF2B",
                "LINE/unknown"="#AAFF55",
                "nonTIR/helitron"="#80FF80",
                "SINE/unknown"="#55FFAA",
                "TIR/CACTA"="#2BFFD4",
                "TIR/MITE"="#00FFFF",
                "TIR/Mutator"="#00D4FF",
                "TIR/PIF_Harbinger"="#00AAFF",
                "TIR/Tc1_Mariner"="#0080FF",
                "TIR/unknown"='#0055FF',
                "TIR/hAT"="#0000D4",
                "repeat_fragment"="#CFCFCF",
                "centromeric_repeat"="#ADAAAB",
                "knob"="#A9A6A7",
                "low_complexity"="#A09D9E",
                "plastid"="#989596",
                "rDNA"="#8C898A",
                "subtelomere"="#807C7D",
                "unknown"='#5B5859')


# DNA_transposon, LTR/LARD, LTR/TR_GAG, LTR/BARE-2, LTR/TRIM, SINE/unknown, TIR/MITE, nonTIR/helitron, repeat_fragment

# Read data and format genome names
long_div <- read.csv(data_file, sep="\t", header=TRUE)

# Filter out rows with 'supfam' values not in TE_colors2
long_div <- long_div[long_div$supfam %in% names(TE_colors2), ]

# Assuming 'genome_size' is a column in 'long_div' and contains the size in base pairs
genome_size_bp <- mean(long_div$genome_size)  # Replace with appropriate filtering if needed

# Convert to Mb and format
genome_size_Mb <- genome_size_bp / 1e6
formatted_genome_size <- format(round(genome_size_Mb), big.mark = ",", scientific = FALSE)

# Create the label using the species accession argument
genome_label <- paste(species_accession, " (", formatted_genome_size, " Mb)", sep="")

# Specify the output file name for the PDF
output_file_name <- paste(species_accession, "_divergence_plot.pdf", sep="")

# Start PDF output
pdf(file = output_file_name, width = 6, height = 4)

# Plot TE divergence
div_p <- long_div %>%
  filter(!grepl('SINE|LINE', supfam)) %>%
  filter(div < 40) %>%
  ggplot(aes(fill = supfam, y = pcnt, x = div)) + 
    geom_bar(position = "stack", stat = "identity", color = "black", linewidth = 0.1) +
    scale_fill_manual(values = TE_colors2) +
    theme_classic() + 
    xlab("Divergence (%)") +
    ylab("Percent of the genome (%)") + 
    labs(fill = "", title = genome_label) +  
    theme(axis.text = element_text(size = 11), axis.title = element_text(size = 12))

# Display the plot
print(div_p)

# Close the PDF device
dev.off()
