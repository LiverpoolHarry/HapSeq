Data analysis pipeline
Given the necessarily large number of libraries that must be processed for direct haplotyping many of the tasks that can be run manually for a single library project need to be automated. We have used Perl scripts to generate multiple commands in the following processes. All programmes were exclusively tested under Debian Linux. All were written for the specific set of samples that were being run and although highly parameterized they may need some editing for different samples in a different directory structure. In particular many scripts iterate over the 29 autosomes in a cow. This loop will need to be modified for other organisms.
Briefly the pipeline consisted of:
1)	make list of file ids and fastq files eg getFastqFiles.pl searches for the required fastq files and extracts sample name and file id from filename 
2)	 run runBWA.pl to map fastq with BWA and output bam files
3)	run runGATK.pl to split sample data by chromosome, merge multiple bam files for same sample, dedup and realign bam files
4)	run buildVcf.pl to use GATK to create a single vcf file for all bamfiles for a chromosome. buildVcf.pl calls removeUninformativeSNPfromVCF.pl to remove uninformative positions and dramatically reduce file size
5)	run runContigMetrics.pl to call the Picard Tools module ContigMetrics.jar to extract scaffolds from bamfiles to gff3 files
6)	run getTargetCut.pl to run samtools  targetcut for all bamfiles and parse output to gff3 files
7)	run runBuildRefhapInput.pl to run BuildRefhapInput.jar to generate input files for the SingleIndividualHaplotyper.jar and then run it to build SNP based haplotypes

The parameters for all Java programmes can be obtained by running them with no parameters.

Get list of Fastq files
The necessary fastq files are likely to be scattered across multiple directories. getFastqFiles.pl takes a list of directories and searches them for fastq.gz files that match a project identifier.  It also gets the md5 checksum for each file, which will be required later when files are submitted to a short read archive.  It prints a list of files to a file that will be used by the getBWA.pl script to pass fastq files to BWA for mapping. getFastqFiles.pl must be edited to insert the path names in the list of paths and the project identifier in the ‘find’ statement. 

Run BWA
A Perl script runBWA.pl takes a list of fastq files from a file and mapping parameters from a config file and maps them against a reference. The script generates a shellscript for each fastq file and uses ssh to send it to the next available machine for mapping, maintaining a load of one job per machine. BWA can be allocated a number of threads to use the available cores on each machine. We used 8 threads per machine and only allocated one job per machine at a time, but this can be changed by editing the $mappingJobsPerMachine parameter in the script.
The user name, directory containing the reference genome, BWA index files, reference genome fasta file, the file containing a list of fastq files, the BWA parameter value pairs, the BWA switches, the platform, a list of machines available by ssh to allocate the jobs to and the path to Picard tools must be set in a config file with default name BWAconfig.txt but a different name can be passed as a parameter of the Perl command line command.  Samtools is required and is assumed to be in the users path.
A log file is generated with all the BWA output. A PerlErrorLog.txt file records any errors thrown by the Perl script. These two files should be consulted in the event of no or unexpected output.
The script runs as much as possible in scratch to minimise network traffic. It assumes that a directory named /scratch exists and creates a subdirectory “/scratch/username/” if it does not already exist. The sequence of commands in the shell script is:
1.	mkdir –p /scratch/username/ 
2.	bwa aln
3.	bwa samse or bwa sampe depending whether it finds an R2 file for each R1 file. The script assumes that R1 will be present in a well formed fastq file name.
4.	Append a count of mapped lines to a file MappedLines.txt 
5.	Convert the sam output from bwa to bam with samtools
6.	Sort the bamfile using picard/SortSam.jar
7.	Index the bamfile with samtools
8.	Move output files from scratch back to the run directory 

BWA –q quality parameter
In the BWA manual the effect of the –q (quality) parameter is described as:
“Parameter for read trimming. BWA trims a read down to argmax_x{\sum_{i=x+1}^l(INT-q_i)} if q_l<INT where l is the original read length.”

SEQanswers pointed to a Perl script TrimBWAstyle.pl from UC Davis http://wiki.bioinformatics.ucdavis.edu/index.php/TrimBWAstyle.pl that implements the BWA algorithm as follows:

BWA quality trimming starts at 5' end and calculates its own quality statistic derived from the machine quality score using the equation BWAscore = q - (MachineScore - 33) so with q = 10 and machine score = 40 then the BWA score =3. Negative values of BWAscore are associated with high quality bases and positive scores with low quality bases.  A BWAscore of zero will be where q + 33 = machine score. It then sums the BWAscores from the 5' end until the sum goes negative. It then trims the read to the first base where it goes negative and uses that for mapping. The intention is that it does not just trim back to the first high quality base but to a region of high quality bases.

Run GATK to merge and dedup files by chromosome
Since the number of files and amount of data can be large we split each bam file output by BWA into multiple bamfiles, one for each chromosome, and then run the GATK dedupping and realignment on files for each chromosome. This generates a much larger number of output bam files.  The next script takes these output bam for simultaneous SNP calling over all bamfiles for a given chromosome at once.
The script mergeandDedupBam.pl takes a list of bam files and a config file default name “configMergeAndDedup.txt” but another filename can be passed as the only parameter to the script.  The list of bam files can be generated using getBamfiles.pl. 
The list of bam files is a tab separated list with a unique identifier for each file in column 0 and the bam file name in column 1. The unique identifier is assumed to start with the sample id separated from the rest of the name (the lane number in our case) with an underscore eg:  “47_L007”. Files with the same sample name and different suffixes will be merged prior to dedupping and realigning. 
The Perl script generates a shellscript that does the following for each chromosome:
1.	Makes a scratch directory if not present
2.	Takes all bam and bai files for a sample and copies them to scratch
3.	Extracts the reads associated with the current chromosomes with samtools
4.	Merges files for same sample using Picard merge
5.	Sorts resulting bam file using Picard SortSamFiles
6.	Removes duplicates using Picard  MarkDuplicates
7.	Realigns around Indels with GATK
8.	Sorts resulting bam file using Picard SortSamFiles
9.	Indexes bamfiles with Samtools
10.	Realigns around Indels with GATK
11.	Sorts resulting bam file using Picard SortSamFiles
12.	Indexes bamfiles with Samtools
13.	Copies output files back to the project directory
14.	Deletes intermediate files from scratch

The Perl script runs each shellfile on the next available machine. As the programmes used are all single threaded, multiple jobs can be sent to each available machine. The default is one, to stay friends with our sysadmin, but this can be changed in the Perl script by changing the parameter $gatkJobsPerMachine.
A logfile with the output from the Picard and GATK programs is generated for each sample x chromosome combination as well as a separate log for the dedupping by Picard tools and a single error log (PerlErrorLog.txt) for the Perl script. Each of these should be consulted in the event of no or unexpected output.

Run GATK to create VCF files
SNP are called on all bam files for all samples for a given chromosome jointly with buildVcf.pl.  The EMIT_ALL_CONFIDENT_SITES option is used so that genotypes for samples that are homozygous for reference alleles are written to the vcf file. This important since, given the low coverage, it cannot be assumed that a sample has the reference allele if there is no data reported. Using this option generates very large vcf files with most positions homozygous or missing in all samples. A Perl script removeUninformativeSNPfromVCF.pl is called by buildVcf.pl after VCF file construction, to remove these uninformative positions and dramatically reduce file size. This script uses bgzip and tabix to zip and index VCF files, these programmes are assumed to be in the path. removeUninformativeSNPfromVCF.pl  is assumed to be in the working directory.
A shellscript is generated to run the programme using all available machines but only one job at a time per machine.
The parameters to buildVcf.pl must be set by editing the perl script. The variables that must be set are all at the top of the script.
For each chromosome the script
1.	Makes a scratch directory
2.	Builds a list of bam files 
3.	Copies each bam and bai file to scratch
4.	Calls GATK unified genotyper
5.	Bgzips and tabix the output vcf file
6.	Runs removeUninformativeSNPfromVCF.pl
7.	Copies vcf file and index back to working directory
8.	Cleans up scratch 


Build scaffolds of reads.
 We have used two strategies for scaffold construction 1) A custom Picard Module “ContigMetrics.jar”; 2) targetcut in samtools.
Building scaffolds with ContigMetrics.jar
 ContigMetrics is not available with the Picard tools distribution. The jar file must be copied into the trunk/dist/ directory of Picard tools.
ContigMetrics is run on all merged bam files for each chromosome generated by the runGATK.pl script by using runContigMetrics.pl. This script must be edited to provide the path to the directory with the bam files to be processed and also the path to a bed file of repeat co-ordinates which will be used to filter out scaffolds wholly within repeats.  The script then builds a system command to extract contigs and scaffolds from each bam file. A large number of parameters can or must be set and these are set by the Perl script. 
The principal output for each sample and chromosome are a gff3 file of scaffolds sample_chr.gff3, a list of contigs sample_chr.contigs, a histogram of counts of each gap length between contigs sample_chr.contig_gaps, a histogram of contig lengths sample_chr.contig_lengths, percentiles of the gap length distribution sample_chr.percentiles.txt, a histogram of read coverage “sample_chr.prop.coverage”. A log file “ScaffoldLog.txt” has summary statistics for all bam files and these can be tabulated by chromosome and by sample by running getScaffoldSummary.pl and passing it the name of the ScaffoldLog.txt file as a parameter. Summary statistics include mean and N50 Scaffold and contig lengths, depth of reads over contigs, proportion of reference sequence covered, numbers of reads processed etc. See Supplementary data SummaryScaffoldStatistics.txt for the data from this project.

Building scaffolds with targetcut in Samtools
Samtools targetcut can be run for all bam files using getTargetCut.pl.  The script must be edited to insert paths to the directory containing bam files. The script parses the output from the targetcut command to produce a gff3 file of scaffold co-ordinates for each input bam file to use as input to BuildRefHap.jar to generate input file for the Single Individual Haplotyper program (SIH).

Running Single Individual Haplotyper program (SIH).
Input files for the SIH programme [1] are generated by BuildRefhapInput.jar from the gff3 scaffold files and the vcf SNP files. BuildRefhapInput.jar is called for all gff3 files by runBuildRefHapInput.pl. As well as generating the .frag and .allvars file for SIH it also compares all overlapping scaffolds to obtain consistency metrics; removes uninformative loci from the SNP data;  removes scaffolds with an excess of heterozygotes using user settable parameters (default 3 heterozygotes or > 20% heterozygotes); removes loci with more than two, two SNP haplotypes. Run  BuildRefhapInput.jar  without any parameters to get a list of required parameters: “java –jar BuildRefhapInput.jar “.
runBuildRefHapInput.pl  assumes that gff3 files have been moved to folders by chromosome with folder names in the format ChrN_Scaffolds. makeChrFolders.pl will search a directory and create the ChrN_Scaffolds folders and move gff3 files into them.
runBuildRefHapInput.pl does the following for each chromosome in turn:
1)	creates of list of all gff3 files
2)	calls BuildRefhapInput.jar which generates a gff3 file for each pair of input gff3 files with a line for each overlap and the consistency of the SNP. It also generates a .frags file and a .allvars file for SIH input.
3)	The overall consistency of all overlaps with and without a switch error correction is obtained by getConsistentProp2.pl which summarizes all the comparison gff3 file and is called by runBuildRefHapInput.pl.
4)	calls getPcentGenomeCoveredByScaffolds.pl to calculate the percentage of the genome that is covered by the complete set of scaffolds.
5)	Tars and zips all of the very large number of gff3 files that report consistency between scaffolds.
6)	Calls SingleIndividualHaplotyper.jar to obtain SNP haplotypes
7)	Once all chromosomes have been processed calls getHapBlockSummary.pl which removes uninformative sites from the .phase file which substantially reduces it s size. It also generates descriptive statistics of the haplotypes for each chromosome: mean, N50 and maximum lengths. 
BuildRefhapInput.jar writes descriptive statistics about the scaffolds to a log file including numbers of scaffolds removed because of heterozygotes, andnumber of loci removed because they had more than two haplotypes. Output getConsistentProp2.pl is appended to the same log file.


Additional Descriptive statistics and annotation
Additional descriptive statistics were obtained using:
1)	parseConsistencyStats.pl parses the consistency stats for each chromosome from the runBuildRefHapInput.pl logfile generated for each chromosome.
2)	parseHetPercentage.pl obtains the numbers and percentages of scaffolds from the runBuildRefHapInput.pl logfile that were deleted because they contained an excess of heterozygotes
3)	runValidateHaplotypes.pl runs ValidateHaplotypes.jar to compare scaffolds haplotypes against haplotypes generated by Beagle using Illumina high density SNP chip data.
4)	runValidateSIHhaplotypes.pl runs ValidateSIHhaplotypes.jar to compare SNP haplotypes generated from sequence data with SIH against haplotypes generated by Beagle using Illumina high density SNP chip data.
5)	annotateSNP.pl uses VCF tools and the Ensembl API. Takes a vcf file and a chromosome number as input and gets a list of genes from Ensembl, extracts the SNP within those genes from the vcf file and submits them to Ensembl for annotation. Writes a new vcf file with the consequences embedded in the annotation (INFO field).  runAnnotateSNP.pl runs annotateSNP.pl for each chromosome. countVariants.pl obtains the counts of each class of annotation in a vcf file.
6)	countPairsOfMisenseSNPOnScaffolds.pl searches genes for pairs of missense SNP within the same gene and counts the number that are captured by scaffolds.  countPairsOfMisenseSNPOnSnpHaplotypes.pl does the same for SNP haplotypes generated by SIH.  sumariseMissenseSNP.pl summarises the data from the individual chromosomes that are output by the first two programmes into a single whole genome table.
7)	getMedianGapBetweenHetsFromVCF.pl was used to obtain the mean and median gaps between informative SNP



References
1.	Suk E-K, McEwen GK, Duitama J, Nowick K, Schulz S, et al. (2011) A comprehensively molecular haplotype-resolved genome of a European individual. Genome Res 21: 1672–1685. doi:10.1101/gr.125047.111.

