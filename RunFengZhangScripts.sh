ReadDir="input/MySeqData/CrisprExperiment/220926_NB501012_0475_AH2GWNAFX5/Data/Intensities/BaseCalls/CTIBRU1"
LibrarySequences="input/library_sequences.csv"

cd "/home/ben/github/crispr-screen-analysis"

source venv/bin/activate

mkdir Results

############### OG

FileKey="OGA01"
ReadFile=$ReadDir"/OGA01_AH2GWNAFX5_S1_L001_R1_001.fastq.gz,"$ReadDir"/OGA01_AH2GWNAFX5_S1_L002_R1_001.fastq.gz,"$ReadDir"/OGA01_AH2GWNAFX5_S1_L003_R1_001.fastq.gz,"$ReadDir"/OGA01_AH2GWNAFX5_S1_L004_R1_001.fastq.gz"
echo $FileKey
python count_spacers.py --fastq=$ReadFile --input=$LibrarySequences
mv statistics.txt "Results/"$FileKey".statistics.txt"
mv library_count.csv "Results/"$FileKey".library_count.csv"
echo "Done with "$FileKey

############### A02-A06

FileKey="A02-A06_All"
ReadFile=$ReadDir"/prepA02_AH2GWNAFX5_S2_L001_R1_001.fastq.gz,"$ReadDir"/prepA02_AH2GWNAFX5_S2_L002_R1_001.fastq.gz,"$ReadDir"/prepA02_AH2GWNAFX5_S2_L003_R1_001.fastq.gz,"$ReadDir"/prepA02_AH2GWNAFX5_S2_L004_R1_001.fastq.gz,"$ReadDir"/prepA03_AH2GWNAFX5_S3_L001_R1_001.fastq.gz,"$ReadDir"/prepA03_AH2GWNAFX5_S3_L002_R1_001.fastq.gz,"$ReadDir"/prepA03_AH2GWNAFX5_S3_L003_R1_001.fastq.gz,"$ReadDir"/prepA03_AH2GWNAFX5_S3_L004_R1_001.fastq.gz,"$ReadDir"/prepA04_AH2GWNAFX5_S4_L001_R1_001.fastq.gz,"$ReadDir"/prepA04_AH2GWNAFX5_S4_L002_R1_001.fastq.gz,"$ReadDir"/prepA04_AH2GWNAFX5_S4_L003_R1_001.fastq.gz,"$ReadDir"/prepA04_AH2GWNAFX5_S4_L004_R1_001.fastq.gz,"$ReadDir"/prepA05_AH2GWNAFX5_S5_L001_R1_001.fastq.gz,"$ReadDir"/prepA05_AH2GWNAFX5_S5_L002_R1_001.fastq.gz,"$ReadDir"/prepA05_AH2GWNAFX5_S5_L003_R1_001.fastq.gz,"$ReadDir"/prepA05_AH2GWNAFX5_S5_L004_R1_001.fastq.gz,"$ReadDir"/prepA06_AH2GWNAFX5_S6_L001_R1_001.fastq.gz,"$ReadDir"/prepA06_AH2GWNAFX5_S6_L002_R1_001.fastq.gz,"$ReadDir"/prepA06_AH2GWNAFX5_S6_L003_R1_001.fastq.gz,"$ReadDir"/prepA06_AH2GWNAFX5_S6_L004_R1_001.fastq.gz"
echo $FileKey
python count_spacers.py --fastq=$ReadFile --input=$LibrarySequences
mv statistics.txt "Results/"$FileKey".statistics.txt"
mv library_count.csv "Results/"$FileKey".library_count.csv"
echo "Done with "$FileKey

############### A02

FileKey="A02"
ReadFile=$ReadDir"/prepA02_AH2GWNAFX5_S2_L001_R1_001.fastq.gz,"$ReadDir"/prepA02_AH2GWNAFX5_S2_L002_R1_001.fastq.gz,"$ReadDir"/prepA02_AH2GWNAFX5_S2_L003_R1_001.fastq.gz,"$ReadDir"/prepA02_AH2GWNAFX5_S2_L004_R1_001.fastq.gz"
echo $FileKey
python count_spacers.py --fastq=$ReadFile --input=$LibrarySequences
mv statistics.txt "Results/"$FileKey".statistics.txt"
mv library_count.csv "Results/"$FileKey".library_count.csv"
echo "Done with "$FileKey

############### A03

FileKey="A03"
ReadFile=$ReadDir"/prepA03_AH2GWNAFX5_S3_L001_R1_001.fastq.gz,"$ReadDir"/prepA03_AH2GWNAFX5_S3_L002_R1_001.fastq.gz,"$ReadDir"/prepA03_AH2GWNAFX5_S3_L003_R1_001.fastq.gz,"$ReadDir"/prepA03_AH2GWNAFX5_S3_L004_R1_001.fastq.gz"
echo $FileKey
python count_spacers.py --fastq=$ReadFile --input=$LibrarySequences
mv statistics.txt "Results/"$FileKey".statistics.txt"
mv library_count.csv "Results/"$FileKey".library_count.csv"
echo "Done with "$FileKey

############### A04

FileKey="A04"
ReadFile=$ReadDir"/prepA04_AH2GWNAFX5_S4_L001_R1_001.fastq.gz,"$ReadDir"/prepA04_AH2GWNAFX5_S4_L002_R1_001.fastq.gz,"$ReadDir"/prepA04_AH2GWNAFX5_S4_L003_R1_001.fastq.gz,"$ReadDir"/prepA04_AH2GWNAFX5_S4_L004_R1_001.fastq.gz"
echo $FileKey
python count_spacers.py --fastq=$ReadFile --input=$LibrarySequences
mv statistics.txt "Results/"$FileKey".statistics.txt"
mv library_count.csv "Results/"$FileKey".library_count.csv"
echo "Done with "$FileKey

############### A05

FileKey="A05"
ReadFile=$ReadDir"/prepA05_AH2GWNAFX5_S5_L001_R1_001.fastq.gz,"$ReadDir"/prepA05_AH2GWNAFX5_S5_L002_R1_001.fastq.gz,"$ReadDir"/prepA05_AH2GWNAFX5_S5_L003_R1_001.fastq.gz,"$ReadDir"/prepA05_AH2GWNAFX5_S5_L004_R1_001.fastq.gz"
echo $FileKey
python count_spacers.py --fastq=$ReadFile --input=$LibrarySequences
mv statistics.txt "Results/"$FileKey".statistics.txt"
mv library_count.csv "Results/"$FileKey".library_count.csv"
echo "Done with "$FileKey

############### A06

FileKey="A06"
ReadFile=$ReadDir"/prepA06_AH2GWNAFX5_S6_L001_R1_001.fastq.gz,"$ReadDir"/prepA06_AH2GWNAFX5_S6_L002_R1_001.fastq.gz,"$ReadDir"/prepA06_AH2GWNAFX5_S6_L003_R1_001.fastq.gz,"$ReadDir"/prepA06_AH2GWNAFX5_S6_L004_R1_001.fastq.gz"
echo $FileKey
python count_spacers.py --fastq=$ReadFile --input=$LibrarySequences
mv statistics.txt "Results/"$FileKey".statistics.txt"
mv library_count.csv "Results/"$FileKey".library_count.csv"
echo "Done with "$FileKey

# Summarize

python SummarizeReadCounts.py -g="input/Brunello_guides.txt" -r "Results"


deactivate
