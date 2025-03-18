library(tidyverse)
library(RColorBrewer)
library(readxl)
library(paletteer)
library(ggplot2)
library(ggsci)
library(ggrepel)
library(ggallin)

#List of WormCat categories for annotating heatmap
annotations <- read_xlsx("wormcat annotations.xlsx")
colnames(annotations) <- gsub(" ", "_", colnames(annotations))

#fix regulated by multple stresses to be part of stress response
annotations$Category_1[annotations$Category_3 == "Unassigned: regulated by multiple stresses"] <- "Stress response"

#stress response genes only
stress_resp_genes <- annotations %>%
  filter(Category_1 == "Stress response") %>%
  select(!Category_1)

#Cleaning up categories
stress_resp_genes$Category_2 <- gsub("Stress response: ", "", stress_resp_genes$Category_2)
stress_resp_genes$Category_2[stress_resp_genes$Category_2 == "Unassigned"] <- "Multiple stressors"
stress_resp_genes$Category_3 <- str_remove_all(stress_resp_genes$Category_3, ".*: ")

#Getting sequencing data
n2_op50 <- read_csv(file = "FC_n2_op_gord_sig.csv") %>%
  select(!c(2,5))
gene_annotations <- read_csv(file = "all_genes.csv")%>%
  select(1,2)

plotdata <- inner_join(n2_op50, stress_resp_genes, by = join_by(ID == Wormbase_ID)) %>%
  left_join(., gene_annotations, by = join_by(Sequence_ID == Sequence.Name)) %>%
  rename(gene = Public.Name)

#Shapes for volcano plot
myshapes <- c(15, 16, 17, 3, 4, 18, 25, 11, 8)
categoryshape <- as.data.frame(table(plotdata$Category_2)) %>%
  arrange(desc(Freq)) %>%
  select(Var1) %>%
  unlist()
names(myshapes) <- categoryshape

#Volcano Plot
png(file = "N2 on strep.png", width = 6, height = 5.5, units = "in", res = 1200)
ggplot(data = plotdata, aes(x = log2FoldChange, y = -log10(padj), 
                            shape = Category_2, col = Category_3)) +
  geom_point(alpha = 0.5, size = 3.25) +
  geom_text_repel(data = plotdata, aes(x = log2FoldChange, y = -log10(padj),
                                       label = gene),
                  size = 2, segment.size = 0.1, 
                  segment.color = "black", max.overlaps = Inf, 
                  inherit.aes = FALSE) +
  theme(legend.key.size = unit(0.5, "cm"), legend.title = element_text(size = 6), legend.text = element_text(size = 5),
        axis.text = element_text(size = 6), axis.title = element_text(size = 8), axis.line = element_line(color = "black", size = 1),
        panel.background = element_blank(), panel.grid = element_line(color = "grey90", size = 0.25, linetype = "dashed"),
        plot.title = element_text(size = 10, hjust = 0.5), plot.subtitle = element_text(size = 10, hjust = 0.5) ) +
  scale_shape_manual(values = myshapes) +
  scale_color_igv() +
  scale_x_continuous(breaks = c(seq(-6, 6, 2)), limits = c(-6,6), trans = pseudolog10_trans) +
  scale_y_log10() +
  guides(shape = guide_legend(order = 1)) +
  labs(
    title = "Expression of significant stress response genes",
    subtitle = expression(paste("in N2 worms after 1 hour of ",italic("S. gordoni"), " infection")),
    x = "Fold Change",
    y = "P-value",
    shape = "Categories",
    col = "Sub-categories"
  )
dev.off()

#heatmap work
plotdatahm <- plotdata %>%
  left_join(., rawcounts, by = "ID")
  

#matrix for heatmap
hm_matrix <- plotdatahm %>%
  select(!ID:Automated_Description) %>%
  data.frame(., row.names = 1) %>%
  as.matrix()

#color code for row annotations based on category 2
row_info <- plotdata %>%
  select(c(gene, Category_2)) %>%
  data.frame(row.names = "gene")
list_colors <- pal_igv()(length(unique(plotdatahm$Category_2)))
names(list_colors) <- unique(plotdatahm$Category_2)
row_colors <- list(Category_2 = list_colors)
#Column names
sample_info <- c("OP50 A", "OP50 B", "OP50 C", "S. gordoni A", "S. gordoni B", "S. gordoni C")

pheatmap(hm_matrix, 
         scale = "row",
         cutree_rows = 3,
         show_rownames = TRUE,
         annotation_row = row_info,
         annotation_colors = row_colors,
         annotation_names_row = FALSE,
         labels_col = sample_info, angle_col = 45,
         filename = "stress genes N2 s gord.png",
         silent = TRUE
)

