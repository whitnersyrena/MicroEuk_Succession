
###############################################################################################  
# Step 1 - run SQM in coassembly mode on single cell samples 
# check ID of bins (automatic output of SQM)
############################################################################################### 

#!/bin/bash
#SBATCH --job-name=SingleCell_sequential
#SBATCH --partition=icemhh
#SBATCH --account=icemhh

## 3 day max run time for public partitions, except 4 hour max runtime for the sandbox partition
#SBATCH --time=02-00:00:00 ## time format is DD-HH:MM:SS
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=400G ## max amount of memory per node you require

#SBATCH --error=job%A.err ## %A - filled with jobid
#SBATCH --output=job%A.out ## %A - filled with jobid
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE,TIME_LIMIT_80
#SBATCH --mail-user=swhitner@hawaii.edu

ml lang/Anaconda3/2024.02-1
source activate SqueezeMeta
cd /home/swhitner/cmaiki_koastore/swhitner/SqueezeMeta

SqueezeMeta.pl -m sequential -p SingleCell_sequential \
-s /home/swhitner/cmaiki_koastore/swhitner/SqueezeMeta/SingleCell.samples \
-f /home/swhitner/cmaiki_koastore/swhitner/kewalos/single_cell --fastnr -t 20


###############################################################################################  
# Step 2 - check contamination/completeness using EukCC of bins identified as "chytrids"
###############################################################################################  

eukcc folder --out outfolder --threads 12 /home/swhitner/cmaiki_koastore/swhitner/SqueezeMeta/SingleCell/results/bins/chytrid_bins


###############################################################################################  
# Step 3 -  run dRep to combine, and dereplicate all chytrid bins to get non-redundant bins 
# use EukCC info from this (eukcc.csv)
############################################################################################### 

dRep dereplicate drep_output -g /Users/syrenawhitner/Desktop/full_run/single_cell/SQM/*.fa --genomeInfo eukcc.csv -comp 1 -con 100 


###############################################################################################  
# Step 4 -  Kaiju on dereplicated/non-redundant draft MAG to assign taxonomy to contigs 
############################################################################################### 

cd /home/swhitner/cmaiki_koastore/swhitner/kaiju

kaiju -z 20 \
  -t nodes.dmp \
  -f /home/swhitner/cmaiki_koastore/swhitner/kaiju/nr_euk/kaiju_db_nr_euk.fmi \
  -i /home/swhitner/cmaiki_koastore/swhitner/SqueezeMeta/SingleCell/results/bins/chytrid_bins/eukcc_outfolder/drep_output/merged.0.fa \
  -o MAG1.kaiju \
  -v

kaiju-addTaxonNames \
  -t nodes.dmp \
  -n names.dmp \
  -i MAG1.kaiju \
  -o MAG1.kaiju.names \
  -r superkingdom,phylum,class,order,family,genus,species

kaiju2table -t nodes.dmp -n names.dmp -r genus -o kaiju_summary.tsv MAG1.kaiju


###############################################################################################  
# Step 5 -  remove contamination contigs (anything that isn't Chytridiomycota/Unclassified)
############################################################################################### 

awk -F'\t' '
  ($1 == "C" && $8 ~ /Chytridiomycota/) || $1 == "U" { print $2 }
' MAG1.kaiju.names > keep_contigs.txt

seqkit grep -f keep_contigs.txt /home/swhitner/cmaiki_koastore/swhitner/SqueezeMeta/SingleCell/results/bins/chytrid_bins/eukcc_outfolder/drep_output/merged.0.fa > /home/swhitner/cmaiki_koastore/swhitner/SqueezeMeta/SingleCell/results/bins/chytrid_bins/eukcc_outfolder/drep_output/CHYTRID_MAG.0.fa


###############################################################################################  
# Step 6 - run EucCC again 
###############################################################################################  

eukcc single --out EukCC_results --threads 12 CHYTRID_MAG_refined.fa --db /home/swhitner/cmaiki_koastore/swhitner/databases/eukccdb/eukcc2_db_ver_1.1


###############################################################################################  
# Step 7 - calculate read coverage of MAG
###############################################################################################  

# calculating read coverage of my MAG 
ml bio/Bowtie2/2.4.5-GCC-11.3.0
ml bio/SAMtools/1.17-GCC-12.2.0

# index assembly

cd /home/swhitner/cmaiki_koastore/swhitner/phyling/MY_MAG
bowtie2-build CHYTRID_MAG_refined.fa CHYTRID_MAG_index

# align reads to MAG - only used original reads that had contigs labeled as "chytridiomycota" 

bowtie2 -x CHYTRID_MAG_index -1 /home/swhitner/cmaiki_koastore/swhitner/kewalos/single_cell/7096_001cat_S171_R1_001.fastq.gz \
-2 /home/swhitner/cmaiki_koastore/swhitner/kewalos/single_cell/7096_001cat_S171_R2_001.fastq.gz -S CHYTRID_MAG.sam

samtools view -bS CHYTRID_MAG.sam | samtools sort -o CHYTRID_MAG.sorted.bam

samtools index CHYTRID_MAG.sorted.bam

samtools depth CHYTRID_MAG.sorted.bam | \
  awk '{sum+=$3} END { print sum/NR }'



###############################################################################################  
# Step 8 -  predict proteins from Draft Chytrid MAG
############################################################################################### 

source activate funannotate 

funannotate clean -i /home/swhitner/cmaiki_koastore/swhitner/SqueezeMeta/SingleCell/results/bins/chytrid_bins/eukcc_outfolder/drep_output/CHYTRID_MAG.0.fa \
-o /home/swhitner/cmaiki_koastore/swhitner/SqueezeMeta/SingleCell/results/bins/chytrid_bins/eukcc_outfolder/drep_output/CHYTRID_MAG_cleaned.0.fa

funannotate mask -i /home/swhitner/cmaiki_koastore/swhitner/SqueezeMeta/SingleCell/results/bins/chytrid_bins/eukcc_outfolder/drep_output/CHYTRID_MAG_cleaned.0.fa \
-o /home/swhitner/cmaiki_koastore/swhitner/SqueezeMeta/SingleCell/results/bins/chytrid_bins/eukcc_outfolder/drep_output/CHYTRID_MAG_masked.0.fa

funannotate predict \
  -i /home/swhitner/cmaiki_koastore/swhitner/SqueezeMeta/SingleCell/results/bins/chytrid_bins/eukcc_outfolder/drep_output/CHYTRID_MAG_masked.0.fa \
  -o chytrid_out_predict \
  --species rhizopus_oryzae \
  --cpus 24


###############################################################################################  
# Step 9 - run trophic analysis 
############################################################################################### 

catastrophy-pipeline --model v10 --outdir ./MY_MAG/results --ncpu 16 ./MY_MAG/*.fasta


###############################################################################################  
# Step 10 -  download JGI reference annotated genomes (predicted protein aa seqs) for Fungi
############################################################################################### 
 # This step there is no code for as taxa were manually selected and downlaoded 

 # rename all contig headers to they are unique for phyling - some of the downloaded genomes have redundant proteins with same names 
 for f in /home/swhitner/cmaiki_koastore/swhitner/phyling/ref_including_euks/*.fasta; do awk -v prefix="$(basename "$f")_" '/^>/{$0=$0"_"++i}1' "$f" > /home/swhitner/cmaiki_koastore/swhitner/phyling/ref_including_euks/renamed/$(basename "$f"); done


###############################################################################################  
# Step 11 -  Phyling pipeline for orthologous group alignment/filter/tree 
############################################################################################### 

source activate BUSCO 

# align
phyling align -I /home/swhitner/cmaiki_koastore/swhitner/phyling/renamed_files -o phyling_align_redo -m fungi_odb12 -t 16 -v

#!/bin/bash
#SBATCH --job-name=phyling_euk_align
#SBATCH --partition=icemhh
#SBATCH --account=icemhh

## 3 day max run time for public partitions, except 4 hour max runtime for the sandbox partition
#SBATCH --time=04-00:00:00 ## time format is DD-HH:MM:SS
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=600G ## max amount of memory per node you require

#SBATCH --error=job%A.err ## %A - filled with jobid
#SBATCH --output=job%A.out ## %A - filled with jobid
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE,TIME_LIMIT_80
#SBATCH --mail-user=swhitner@hawaii.edu

ml lang/Anaconda3/2024.02-1
source activate BUSCO
/home/swhitner/cmaiki_koastore/swhitner/phyling

phyling align -I /home/swhitner/cmaiki_koastore/swhitner/phyling/alL_euks_renamed -o all_euk_yalignk -m eukaryota_odb12 -t 16 -v

#tree 

#!/bin/bash
#SBATCH --job-name=full_euk_tree
#SBATCH --partition=icemhh
#SBATCH --account=icemhh

## 3 day max run time for public partitions, except 4 hour max runtime for the sandbox partition
#SBATCH --time=04-00:00:00 ## time format is DD-HH:MM:SS
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=600G ## max amount of memory per node you require

#SBATCH --error=job%A.err ## %A - filled with jobid
#SBATCH --output=job%A.out ## %A - filled with jobid
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE,TIME_LIMIT_80
#SBATCH --mail-user=swhitner@hawaii.edu

ml lang/Anaconda3/2024.02-1
source activate BUSCO
cd /home/swhitner/cmaiki_koastore/swhitner/phyling

phyling tree -I all_eu_alignk -o full_euk_tree -c -p -M iqtree -t 32


###############################################################################################  
# Step 12 -  Make tree in R --> also available as R script 
############################################################################################### 
``` R {
library(phangorn)
library(ape)
library(ggtree)
library(treeio)
library(ggplot2)

tip_info <- read.csv("/Users/syrenawhitner/Desktop/trees/tip_info.csv")
tip_info$tip_labels <- gsub("\\.aa\\.fasta\\.gz$", "", tip_info$tip_labels)

tr <- read.tree("/Users/syrenawhitner/Downloads/ufboot (2).contree")

tr <- root(tr, outgroup = c("Physo3_GeneCatalog_proteins_20110401"),
           resolve.root = TRUE)

tr <- ladderize(tr)

p_black <- ggtree(tr, color = "black") + theme_tree2() +
  geom_tree(linewidth = 0.5, color = "black") +
  geom_tiplab(size = 3, color = "black") +
  geom_text2(
    aes(subset = !isTip & !is.na(as.numeric(label)) & as.numeric(label) >= 70,
        label = round(as.numeric(label))),
    color = "black", hjust = -0.2, size = 3
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.05)))

print(p_black)

} ```



































