# marine_fungal_parasites
This pipeline was developed for the identification and annotation of fungal parasites that infect marine phytoplankton. 

#1. analyze amplicon data using QIIME2 (command line) 
#2. differential expression and basic statistics (R) 
#3. viaualize results (R) 

# Step #1 - Amplicon analysis using QIIME2
Import sequences (Casava One Eight format - sequences format as NTC-M1_S173_L001_R2_001.fastq.gz)
```
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path /Users/syrenawhitner/Desktop/kewalos/KML_2024/prelim/sequences/M2 \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path demux-paired-end.qza
```

Trim primers using cut adapt 
in this case we are using - ITSOF / fITS7R       
f - ACTTGGTCATTTAGAGGAAGT      
r - CAAAGATTCGATGAYTCAC   

```
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences demux-paired-end.qza \
  --p-front-f ACTTGGTCATTTAGAGGAAGT \
  --p-front-r CAAAGATTCGATGAYTCAC \
  --o-trimmed-sequences trimmed-seqs.qza \
  --verbose
```

View demux summary 
```
qiime demux summarize \
  --i-data trimmed-seqs.qza \
  --o-visualization trimmed-seqs.qzv
```


Denoise with dada2 - choose truncation this length based on demux summary / read quality 
```
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs trimmed-seqs.qza \
  --p-trunc-len-f 180 \
  --p-trunc-len-r 180 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza
```

Import and train the eukaryome database classifier - for this iteration we are using 97% similarity values of the EUKARYOME database
IF CLASSIFIER IS ALREADY TRAINED --> Skip to "assign taxonomy"

Import db reference sequences
```
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path /Users/syrenawhitner/Desktop/kewalos/KML_2024/prelim/databases/eukaryome/QIIME2_EUK_SSU_v1.9/QIIME2_EUK_SSU_v1.9.fasta \
  --output-path euk-ref-seqs.qza
```

Import db reference taxonomy
```
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path /Users/syrenawhitner/Desktop/kewalos/KML_2024/prelim/databases/eukaryome/QIIME2_EUK_SSU_v1.9/QIIME2_EUK_SSU_v1.9.tsv \
  --output-path euk-ref-taxonomy.qza
```

Train classifier 
```
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads euk-ref-seqs.qza \
  --i-reference-taxonomy euk-ref-taxonomy.qza \
  --o-classifier euk-classifier.qza
```

Assign taxonomy with your new classifier
```
qiime feature-classifier classify-sklearn \
  --i-classifier /Users/syrenawhitner/Desktop/kewalos/KML_2024/prelim/qiime2/M1/unite-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza
```

Tabulate and visualize the taxonomy information for the classified sequences. At this step, you can decide if there are various sequences you wish to remove from the rep-seqs/table files. 
```
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv
```

At this stage, we can export the tax.qza, feature table.qaz, and rep_seqs.qza and import into R for further analysis and visualization
Export feature table 
```
qiime tools export \
  --input-path table.qza \
  --output-path exported-feature-table
```

Convert feature table to TSV format
```
biom convert \
  --input-fp exported-feature-table/feature-table.biom \
  --output-fp exported-feature-table/feature-table.tsv \
  --to-tsv
```

Export representative sequences
```
qiime tools export \
  --input-path rep-seqs.qza \
  --output-path exported-rep-seqs
```

Export taxonomy
```
qiime tools export \
  --input-path taxonomy.qza \
  --output-path exported-taxonomy
```
THE FOLLOWING FOUR STEPS ARE ALL OPTIONAL VISUALISATIONS AND STATISTICS

Make a phylogenetic tree
```
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
```

Summarize feature table to choose an appropriate sampling depth for core diversity metrics
```
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table-summary.qzv
```

Core diversity metrics (sampling depth selected in previous step)
```
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth 15000 \
  --m-metadata-file sample-metadata.tsv \
  --output-dir core-metrics-results
```

Make a taxonomic barplot 
```
qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization taxa-bar-plots.qzv
```
