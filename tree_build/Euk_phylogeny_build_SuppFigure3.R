#################################################################################################
##                                                                                             ##
##                  Building phylogeny from Single Cell Isolation analysis                     ##
##                                                                                             ##
#################################################################################################

library(phangorn)
library(ape)
library(ggtree)
library(treeio)
library(ggplot2)

tip_info <- read.csv("/Users/syrenawhitner/Desktop/fungi_500_tree/Supplementary tables - Sheet11.csv")

tip_info$tip_labels <- gsub("\\.aa\\.fasta\\.gz$", "", tip_info$tip_labels)

tr <- read.tree("/Users/syrenawhitner/Downloads/ufboot (2).contree")

tr <- root(tr, outgroup = c("Physo3_GeneCatalog_proteins_20110401"),
           resolve.root = TRUE)

tr <- ladderize(tr)

p<- ggtree(tr, color = "black") + theme_tree2() +
  geom_tree(linewidth = 0.5, color = "black") +
  geom_tiplab(size = 3, color = "black") +
  geom_text2(
    aes(subset = !isTip & !is.na(as.numeric(label)) & as.numeric(label) >= 70,
        label = round(as.numeric(label))),
    color = "black", hjust = -0.2, size = 3
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.05)))

print(p)





































