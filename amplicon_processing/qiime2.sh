##### importing samples already trimmed in dada2 into qiime2 for the taxonomic classification

##### import pre-trimmed sequences 
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path /Users/syrenawhitner/Desktop/filt/dada_trimmed_seqs \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path demux-paired-end.qza


##### generate table 
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux-paired-end.qza \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza


##### important docs to train classifier 
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path /Users/syrenawhitner/Desktop/kewalos/KML_2024/prelim/databases/eukaryome/QIIME2_EUK_SSU_v1.9/QIIME2_EUK_SSU_v1.9.fasta \
  --output-path euk-ref-seqs.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path /Users/syrenawhitner/Desktop/kewalos/KML_2024/prelim/databases/eukaryome/QIIME2_EUK_SSU_v1.9/QIIME2_EUK_SSU_v1.9.tsv \
  --output-path euk-ref-taxonomy.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads euk-ref-seqs.qza \
  --i-reference-taxonomy euk-ref-taxonomy.qza \
  --o-classifier euk-classifier.qza


##### classify sequences 
qiime feature-classifier classify-sklearn \
  --i-classifier /Users/syrenawhitner/Desktop/kewalos/KML_2024/Dec2024/qiime2/euk-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza


##### export feature table and rep seqs 
qiime tools export \
  --input-path table.qza \
  --output-path exported-feature-table


##### convert feature table to TSV format
biom convert \
  --input-fp exported-feature-table/feature-table.biom \
  --output-fp exported-feature-table/feature-table.tsv \
  --to-tsv


##### export representative sequences
qiime tools export \
  --input-path rep-seqs.qza \
  --output-path exported-rep-seqs


#####export taxonomy
qiime tools export \
  --input-path taxonomy.qza \
  --output-path exported-taxonomy
