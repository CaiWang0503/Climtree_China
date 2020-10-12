HOMEFOLDER="/Users/wang/Documents/climtree/climtree_2016data/"
SEQS="seqs/"
cd ${HOMEFOLDER}${SEQS}

####################################################################################################
# Set variable values
####################################################################################################
DAME="/usr/local/bin/DAMe/bin/"
SUMASIM=97
MINPCR_1=1
MINREADS_1=1
MINLEN=300
MAXLEN=320
POOLS=6
ARTHMINPROB=0.8
####################################################################################################
# 1.1 Use AdapterRemoval to trim Illumina sequencing adapters. (there are a few adapters still left in the raw data1801)
####################################################################################################
AdapterRemoval  --file1 C1_S8_L001_R1_002.fastq.gz --file2 C1_S8_L001_R2_002.fastq.gz --basename AdaptermvP1 --threads 3 --trimns
AdapterRemoval  --file1 C2_S9_L001_R1_002.fastq.gz --file2 C2_S9_L001_R2_002.fastq.gz --basename AdaptermvP2 --threads 3 --trimns
AdapterRemoval  --file1 C3_S10_L001_R1_002.fastq.gz --file2 C3_S10_L001_R2_002.fastq.gz --basename AdaptermvP3 --threads 3 --trimns
AdapterRemoval  --file1 P1_S13_L001_R1_001.fastq.gz --file2 P1_S13_L001_R2_001.fastq.gz --basename AdaptermvP4 --threads 3 --trimns
AdapterRemoval  --file1 P2_S14_L001_R1_001.fastq.gz --file2 P2_S14_L001_R2_001.fastq.gz --basename AdaptermvP5 --threads 3 --trimns
AdapterRemoval  --file1 P3_S15_L001_R1_001.fastq.gz --file2 P3_S15_L001_R2_001.fastq.gz --basename AdaptermvP6 --threads 3 --trimns


####################################################################################################
# 1.2 Use Sickle to trim away low-quality nucleotides from the 3' ends of the reads (does not try to fix errors)
####################################################################################################
sickle pe -f AdaptermvP1.pair1.truncated -r AdaptermvP1.pair2.truncated -o sickle_P1_R1.fq -p sickle_P1_R2.fq -s sickle_Single_P1.fq -t sanger > sickleP1.out
sickle pe -f AdaptermvP2.pair1.truncated -r AdaptermvP2.pair2.truncated -o sickle_P2_R1.fq -p sickle_P2_R2.fq -s sickle_Single_P2.fq -t sanger > sickleP2.out
sickle pe -f AdaptermvP3.pair1.truncated -r AdaptermvP3.pair2.truncated -o sickle_P3_R1.fq -p sickle_P3_R2.fq -s sickle_Single_P3.fq -t sanger > sickleP3.out
sickle pe -f AdaptermvP4.pair1.truncated -r AdaptermvP4.pair2.truncated -o sickle_P4_R1.fq -p sickle_P4_R2.fq -s sickle_Single_P4.fq -t sanger > sickleP4.out
sickle pe -f AdaptermvP5.pair1.truncated -r AdaptermvP5.pair2.truncated -o sickle_P5_R1.fq -p sickle_P5_R2.fq -s sickle_Single_P5.fq -t sanger > sickleP5.out
sickle pe -f AdaptermvP6.pair1.truncated -r AdaptermvP6.pair2.truncated -o sickle_P6_R1.fq -p sickle_P6_R2.fq -s sickle_Single_P6.fq -t sanger > sickleP6.out

# Now you can remove the AdapterRemoval output files
rm Adaptermv*.*

####################################################################################################
# 1.3 Use SPAdes to correct sequencing and PCR errors, via BayesHammer
####################################################################################################
spades.py --only-error-correction -1 sickle_P1_R1.fq -2 sickle_P1_R2.fq -o SPAdes_hammer_P1 -t 4
spades.py --only-error-correction -1 sickle_P2_R1.fq -2 sickle_P2_R2.fq -o SPAdes_hammer_P2 -t 4
spades.py --only-error-correction -1 sickle_P3_R1.fq -2 sickle_P3_R2.fq -o SPAdes_hammer_P3 -t 4
spades.py --only-error-correction -1 sickle_P4_R1.fq -2 sickle_P4_R2.fq -o SPAdes_hammer_P4 -t 4
spades.py --only-error-correction -1 sickle_P5_R1.fq -2 sickle_P5_R2.fq -o SPAdes_hammer_P5 -t 4
spades.py --only-error-correction -1 sickle_P6_R1.fq -2 sickle_P6_R2.fq -o SPAdes_hammer_P6 -t 4

rm sickle_*.fq

####################################################################################################
# 1.4 Use PandaSeq to merge the paired reads
####################################################################################################
pandaseq -f SPAdes_hammer_P1/corrected/sickle_P1_R1.00.0_0.cor.fastq.gz -r SPAdes_hammer_P1/corrected/sickle_P1_R2.00.0_0.cor.fastq.gz -A simple_bayesian -d bfsrk -F -N -T 7 -g pandaseq_log_P1.txt -w sickle_cor_panda_P1.fastq
pandaseq -f SPAdes_hammer_P2/corrected/sickle_P2_R1.00.0_0.cor.fastq.gz -r SPAdes_hammer_P2/corrected/sickle_P2_R2.00.0_0.cor.fastq.gz -A simple_bayesian -d bfsrk -F -N -T 7 -g pandaseq_log_P2.txt -w sickle_cor_panda_P2.fastq
pandaseq -f SPAdes_hammer_P3/corrected/sickle_P3_R1.00.0_0.cor.fastq.gz -r SPAdes_hammer_P3/corrected/sickle_P3_R2.00.0_0.cor.fastq.gz -A simple_bayesian -d bfsrk -F -N -T 7 -g pandaseq_log_P3.txt -w sickle_cor_panda_P3.fastq
pandaseq -f SPAdes_hammer_P4/corrected/sickle_P4_R1.00.0_0.cor.fastq.gz -r SPAdes_hammer_P4/corrected/sickle_P4_R2.00.0_0.cor.fastq.gz -A simple_bayesian -d bfsrk -F -N -T 7 -g pandaseq_log_P4.txt -w sickle_cor_panda_P4.fastq
pandaseq -f SPAdes_hammer_P5/corrected/sickle_P5_R1.00.0_0.cor.fastq.gz -r SPAdes_hammer_P5/corrected/sickle_P5_R2.00.0_0.cor.fastq.gz -A simple_bayesian -d bfsrk -F -N -T 7 -g pandaseq_log_P5.txt -w sickle_cor_panda_P5.fastq
pandaseq -f SPAdes_hammer_P6/corrected/sickle_P6_R1.00.0_0.cor.fastq.gz -r SPAdes_hammer_P6/corrected/sickle_P6_R2.00.0_0.cor.fastq.gz -A simple_bayesian -d bfsrk -F -N -T 7 -g pandaseq_log_P6.txt -w sickle_cor_panda_P6.fastq


####################################################################################################
# 1.5 gzip the final fastq files and remove working files
####################################################################################################
gzip sickle_cor_panda_P1.fastq
gzip sickle_cor_panda_P2.fastq
gzip sickle_cor_panda_P3.fastq
gzip sickle_cor_panda_P4.fastq
gzip sickle_cor_panda_P5.fastq
gzip sickle_cor_panda_P6.fastq
# Now you can remove the SPAdes output
rm -rf SPAdes_hammer*/
# remove pandaseq_log_txt files
rm pandaseq_log_*.txt

####################################################################################################
####################################################################################################
# STAGE TWO:  DAMe filtering the merged reads to keep only the sequences that appear in
#             multiple PCRs and in multiple copies
####################################################################################################
####################################################################################################


####################################################################################################
# 2.1 Place the libraries in a folder for that experiment (B) and then into separate folders
####################################################################################################

mkdir pool1/
mkdir pool2/
mkdir pool3/
mkdir pool4/
mkdir pool5/
mkdir pool6/

mv sickle_cor_panda_P1.fastq.gz pool1/
mv sickle_cor_panda_P2.fastq.gz pool2/
mv sickle_cor_panda_P3.fastq.gz pool3/
mv sickle_cor_panda_P4.fastq.gz pool4/
mv sickle_cor_panda_P5.fastq.gz pool5/
mv sickle_cor_panda_P6.fastq.gz pool6/


####################################################################################################
# 2.2 Sort through each fastq file and determine how many of each tag pair is in each fastq file. Each fastq file takes < 1 min
####################################################################################################

${DAME}sort.py -h #

cd pool1/
python ${DAME}sort.py -fq sickle_cor_panda_P1.fastq.gz -p ${HOMEFOLDER}Primers_COILeray.txt -t ${HOMEFOLDER}TagsA_climtree2016_COI.txt
cd ../pool2/
python ${DAME}sort.py -fq sickle_cor_panda_P2.fastq.gz -p ${HOMEFOLDER}Primers_COILeray.txt -t ${HOMEFOLDER}TagsA_climtree2016_COI.txt
cd ../pool3/
python ${DAME}sort.py -fq sickle_cor_panda_P3.fastq.gz -p ${HOMEFOLDER}Primers_COILeray.txt -t ${HOMEFOLDER}TagsA_climtree2016_COI.txt
cd ../pool4/
python ${DAME}sort.py -fq sickle_cor_panda_P4.fastq.gz -p ${HOMEFOLDER}/Primers_COILeray.txt -t ${HOMEFOLDER}/TagsB_climtree2016_COI.txt
cd ../pool5/
python ${DAME}sort.py -fq sickle_cor_panda_P5.fastq.gz -p ${HOMEFOLDER}/Primers_COILeray.txt -t ${HOMEFOLDER}/TagsB_climtree2016_COI.txt
cd ../pool6/
python ${DAME}sort.py -fq sickle_cor_panda_P6.fastq.gz -p ${HOMEFOLDER}/Primers_COILeray.txt -t ${HOMEFOLDER}/TagsB_climtree2016_COI.txt
####################################################################################################
# 2.3 Make an overview of the tag combinations.
####################################################################################################
cd ${HOMEFOLDER}seqs/
python ${DAME}splitSummaryByPSInfo.py -p ${HOMEFOLDER}/PSinfoAB_climtree_COI.txt -l 1 -s pool1/SummaryCounts.txt -o splitSummaryByPSInfo_P1.txt
python ${DAME}splitSummaryByPSInfo.py -p ${HOMEFOLDER}/PSinfoAB_climtree_COI.txt -l 2 -s pool2/SummaryCounts.txt -o splitSummaryByPSInfo_P2.txt
python ${DAME}splitSummaryByPSInfo.py -p ${HOMEFOLDER}/PSinfoAB_climtree_COI.txt -l 3 -s pool3/SummaryCounts.txt -o splitSummaryByPSInfo_P3.txt
python ${DAME}splitSummaryByPSInfo.py -p ${HOMEFOLDER}/PSinfoAB_climtree_COI.txt -l 4 -s pool4/SummaryCounts.txt -o splitSummaryByPSInfo_P4.txt
python ${DAME}splitSummaryByPSInfo.py -p ${HOMEFOLDER}/PSinfoAB_climtree_COI.txt -l 5 -s pool5/SummaryCounts.txt -o splitSummaryByPSInfo_P5.txt
python ${DAME}splitSummaryByPSInfo.py -p ${HOMEFOLDER}/PSinfoAB_climtree_COI.txt -l 6 -s pool6/SummaryCounts.txt -o splitSummaryByPSInfo_P6.txt

####################################################################################################
# 2.4 Generate a heatmap of the tag pair read information.
####################################################################################################
cd ${HOMEFOLDER}/seqs/pool1/
Rscript --vanilla --verbose ${HOMEFOLDER}/heatmap_for_teaching.R
cd ${HOMEFOLDER}/seqs/pool2/
Rscript --vanilla --verbose ${HOMEFOLDER}/heatmap_for_teaching.R
cd ${HOMEFOLDER}/seqs/pool3/
Rscript --vanilla --verbose ${HOMEFOLDER}/heatmap_for_teaching.R
cd ${HOMEFOLDER}/seqs/pool4/
Rscript --vanilla --verbose ${HOMEFOLDER}/heatmap_for_teaching.R
cd ${HOMEFOLDER}/seqs/pool5/
Rscript --vanilla --verbose ${HOMEFOLDER}/heatmap_for_teaching.R
cd ${HOMEFOLDER}/seqs/pool6/
Rscript --vanilla --verbose ${HOMEFOLDER}/heatmap_for_teaching.R


# then move heatmaps from inside the pool folders into folder_A/
cd ${HOMEFOLDER}/seqs/
mv pool1/heatmap.pdf heatmap_P1.pdf
mv pool2/heatmap.pdf heatmap_P2.pdf
mv pool3/heatmap.pdf heatmap_P3.pdf
mv pool4/heatmap.pdf heatmap_P4.pdf
mv pool5/heatmap.pdf heatmap_P5.pdf
mv pool6/heatmap.pdf heatmap_P6.pdf


####################################################################################################
# 2.5 Check if the PCRs were carried out correctly.
####################################################################################################
cd ${HOMEFOLDER}/seqs/
mkdir Filter_min${MINPCR_1}PCRs_min${MINREADS_1}copy # make a folder for the output

python ${DAME}filter.py -psInfo ${HOMEFOLDER}/PSinfoAB_climtree_COI.txt -x ${POOLS} -y ${MINPCR_1} -p ${POOLS} -t ${MINREADS_1} -l ${MINLEN} -o Filter_min${MINPCR_1}PCRs_min${MINREADS_1}copy
python ${DAME}RSI.py --explicit -o RSI_output_climtree2016_min${MINPCR_1}PCR_min${MINREADS_1}.txt Filter_min${MINPCR_1}PCRs_min${MINREADS_1}copy/Comparisons_${POOLS}PCRs.txt

####################################################################################################
# 2.6 Choose thresholds for filtering by minimum number of PCRS (pools) and copies.
####################################################################################################
cd ${HOMEFOLDER}seqs/Filter_min${MINPCR_1}PCRs_min${MINREADS_1}copy/
python ${DAME}plotLengthFreqMetrics_perSample.py -f FilteredReads.fna -p ${HOMEFOLDER}PSinfoAB_climtree_COI.txt -n 6

####################################################################################################
# 2.7 Filter the reads by presence in a user-defined minimum number of PCRS (pools) and copies.
####################################################################################################
python ${DAME}filter.py -h

MINPCR=2
MINREADS=30
cd ${HOMEFOLDER}/seqs/
mkdir Filter_min${MINPCR}PCRs_min${MINREADS}copies/
python ${DAME}filter.py -psInfo ${HOMEFOLDER}/PSinfoAB_climtree_COI.txt -x 6 -y ${MINPCR} -p 6 -t ${MINREADS} -l ${MINLEN} -o Filter_min${MINPCR}PCRs_min${MINREADS}copies

cd ${HOMEFOLDER}/seqs/Filter_min${MINPCR}PCRs_min${MINREADS}copies
python ${DAME}plotLengthFreqMetrics_perSample.py -f FilteredReads.fna -p ${HOMEFOLDER}/PSinfoAB_climtree_COI.txt -n 6

####################################################################################################
# 2.8 Prepare FilteredReads.fna file for OTU clustering and remove chimeras using vsearch
####################################################################################################

python ${DAME}convertToUSearch.py -h

MINPCR=2 # these commands are to make it possible to process multiple filter.py outputs, which are saved in different Filter_min folders
MINREADS=30

# vsearch requires a few seconds per library when applied to FilteredReads.forsumaclust.fna
cd ${HOMEFOLDER}/seqs/Filter_min${MINPCR}PCRs_min${MINREADS}copies
python ${DAME}convertToUSearch.py -i FilteredReads.fna -lmin ${MINLEN} -lmax ${MAXLEN}
gsed 's/ count/;size/' FilteredReads.forsumaclust.fna > FilteredReads.forvsearch.fna
vsearch --sortbysize FilteredReads.forvsearch.fna --output FilteredReads.forvsearch_sorted.fna
vsearch --uchime_denovo FilteredReads.forvsearch_sorted.fna --nonchimeras FilteredReads.forvsearch_sorted_nochimeras.fna
gsed 's/;size/ count/' FilteredReads.forvsearch_sorted_nochimeras.fna > FilteredReads.forsumaclust.nochimeras.fna

# remove vsearch uchime working files
rm FilteredReads.forsumaclust.fna
rm FilteredReads.forvsearch.fna
rm FilteredReads.forvsearch_sorted.fna
rm FilteredReads.forvsearch_sorted_nochimeras.fna
####################################################################################################
# 2.8 Analyse FilteredReads.fna for pairwise similarities to choose similarity threshold for Sumaclust
####################################################################################################

cd ${HOMEFOLDER}/seqs/Filter_min${MINPCR}PCRs_min${MINREADS}copies
python ${DAME}assessClusteringParameters.py -i FilteredReads.forsumaclust.nochimeras.fna -mint 0.8 -minR 0.6 -step 0.05 -t 4 -o COIclusterassess_mint08_minR06_step005_Filter_min${MINPCR}PCRs_min${MINREADS}copies.pdf
####################################################################################################
# 2.9 Sumaclust clustering and convert Sumaclust output to table format
####################################################################################################
# 97% sumaclust
SUMASIM=97 # if there is no SUMASIM value
cd ${HOMEFOLDER}/seqs/Filter_min${MINPCR}PCRs_min${MINREADS}copies
# cluster reads into OTUs. Output file assigns each read to a cluster.
sumaclust -t .${SUMASIM} -e FilteredReads.forsumaclust.nochimeras.fna > OTUs_${SUMASIM}_sumaclust.fna
# create OTU table and create a fasta file for OTU representative sequences
python ${DAME}tabulateSumaclust.py -i OTUs_${SUMASIM}_sumaclust.fna -o table_climtree2016_${SUMASIM}.txt -blast

mv table_climtree2016_${SUMASIM}.txt.blast.txt table_climtree2016_${SUMASIM}.fas # chg filename suffix to *.fas