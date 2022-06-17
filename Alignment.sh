# Set some variables.
ReferenceInputFile="input/Reference.AmbiguousBarcode.fasta"
ReadInputFile="input/sequenced_library/CRISPR20220513_S1_L001_R1_001.fastq.gz"
AlignmentPrefix="CRISPR20220513_S1_L001_R1_001.fastq.gz_alignment"
OutputSubdir=$AlignmentPrefix"/"

# Generate some variables
ReferenceIndex="${ReferenceInputFile/.fasta/.mmi}"
AlignmentBam=$OutputSubdir$AlignmentPrefix.bam

mkdir $OutputSubdir

# Make reference file
minimap2 -d $ReferenceIndex $ReferenceInputFile

# Do Alignment with Minimap.
# Pipe alignment into "samtools view" to convert to bam.
# Pipe into "samtools sort" to sort I guess.
# note: minimap2, and samtools operations can accept a # of threads to use.
minimap2 -a $ReferenceIndex $ReadInputFile | samtools view -b | samtools sort -o $AlignmentBam

# Index the alignment for viewing
samtools index $AlignmentBam

# You can view with IGV

