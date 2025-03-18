library(tidyverse)
library('tximport')
library('DESeq2')
library(janitor)

# Sample info for DeSeq2
samples <- read_csv('ref_data/sample_info.csv')
# Transform condition to factor as this is what DESeq2 expects
samples$condition <- factor(samples$condition, levels = c("control", "gordoni", "catalase"))

rownames(samples) <- samples$sample

tx2gene <- read_csv('ref_data/simplemine_results.csv', na = 'N.A.') |>
  clean_names() |> 
  select(c('public_name', 'transcript')) |> 
  drop_na(transcript) |> 
  relocate(transcript, public_name)
dir <- 'raw_counts/quants/'

# Associates each sample to the appropriate file
files <- data.frame(
  sampleName = samples$sample,
  fileName = paste0(dir, samples$sample, "_quant/quant.sf") 
)

txi <- tximport(files$fileName, type='salmon', tx2gene = tx2gene)

ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ condition)
dds <- DESeq(ddsTxi)

# zip-10 comparisions of all conditions
gordoni_vs_control <- results(dds, contrast = c("condition", "gordoni", "control"))
catalase_vs_control <- results(dds, contrast = c("condition", "catalase", "control"))
gordoni_vs_catalase <- results(dds, contrast = c("condition", "gordoni", "catalase"))

# Apply filtering criteria: fold change > 1.5 or < -1.5 and adjusted p-value < 0.05
filter_results <- function(results) {
  results <- as.data.frame(results)
  filtered_results <- results[
    (results$log2FoldChange > log2(1.5) | results$log2FoldChange < -log2(1.5)) & results$padj < 0.05,
  ]
  return(filtered_results)
}

write.csv(as.data.frame(results), file = "zip_10.csv")
