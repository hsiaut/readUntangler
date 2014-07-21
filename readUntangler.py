#!//data/bifx/software/Python-2.7.1/python

#######################################################################
# readUntangler.py
#
#
# reads in fasta file
# folds fasta file to 80 chars per line and removes any escape characters
# uses blast to retrieve relevant reads from specified readset file
# uses bwa to align reads to input fasta file
# calls variants using freebayes or mpileup
# runs extractHAIRS to extrac the relevant reads and their variant calls
# filters out reads with a low number of calls (default: 6 for pacbio)
# runs hapCUT to generate the haplotype blocks
# reads hapCUT output file and generates fasta and gff files for the two haplotypes
#
#
#######################################################################
#
# Written by Tim Hsiau for Solazyme July 2013
#
########################################################################

##Libraries##
import subprocess, argparse, time, os

##script settings##
#location of binaries
#get location of script
scriptdir = os.path.dirname(__file__) +"/readUntanglerBinaries"
if scriptdir == "":
	scriptdir = "."
scriptdir += "/"


def run_bwa(ref, readset, outfile):
	"""
	runs bwa mem by: creating an index of the reference sequence
	running bwa mem
	piping output to samtools to get rid of unmapped reads
	writing to file
	"""

def clean_fasta(filename):
	"""
	cleans input fasta files by folding them
	and removing ^M escape characters
	^M is the same as \r
	"""
	foldedfile = filename+".folded.fa"
	foldcmd = "fold "+filename+" > "+foldedfile
	subprocess.call(foldcmd, shell=True)
	cleanedfile = filename+".cleaned.fa"
	sedcmd = "sed -e 's/\r//g' "+foldedfile+" > "+cleanedfile
	subprocess.call(sedcmd, shell=True)
	return cleanedfile

def readsetchooser(readsetname):
	#function used to switch between readsets
	#returns 0 if no readset found
	return ( "SRR317823_1.fastq" ,"SRR317823_2.fastq")

def run_bwa(reference, readset):
	"""
	creates the bwa index for the reference file
	runs bwa mem, current settings are 20 threads
	pipes bwa mem output to samtools, which gets rid of unmapped reads
	writes output to aligned.sam file
	"""
	bwaindexcmd = "bwa index "+reference
	subprocess.call(bwaindexcmd, shell=True)
	alignedsam = reference+".aligned.sam"
	bwacmd = "bwa mem -t 12 "+reference+" "+readset
	bwacmd += " | samtools view -SF 4 - > "+alignedsam
	subprocess.call(bwacmd, shell=True)
	return alignedsam

def sam_to_sorted_bam(reference, samfile):
	"""
	converts sam to sorted bam file
	uses the samtools in the local directory
	"""
	faidxcmd= "samtools faidx "+reference
	subprocess.call(faidxcmd, shell=True)
	bamfile = reference+".bam"
	subprocess.call("samtools view -bt "+reference+".fai "+samfile+" > "+bamfile, shell=True)
	sortedbam = reference+".sorted" #samtools automatically adds the ".bam" extension
	subprocess.call("samtools sort "+bamfile+" "+sortedbam, shell=True)

	#todo: if output file passes QC, delete intermediate bam file
	return sortedbam+".bam"

def run_var_caller(reference, sortedbam):
	"""
	calls freebayes version
	"""
	vcfoutfile=reference+".vcf"
	freebayescmd = scriptdir+"freebayes.binary --fasta-reference "+reference+" "+sortedbam+" > "+vcfoutfile
	subprocess.call(freebayescmd, shell=True)
	return vcfoutfile

def run_haplotyper(reference, vcffile, sortedbam, filterthreshold):
	"""
	runs extractHAIRS and hapcut
	"""
	fragoutfile = reference+".fragout"
	extractHAIRScmd=scriptdir+"extractHAIRS --VCF "+vcffile+" --bam "+sortedbam+" --indels 1 --ref "+reference+" > "+fragoutfile
	subprocess.call(extractHAIRScmd, shell=True)
	filteredfragfile = reference+".filtered.fragout"
	filterthreshold = str(filterthreshold)
	filtercmd = "awk '{ if (length($NF) > "+filterthreshold+" ) print;}' "+fragoutfile+" > "+filteredfragfile
	subprocess.call(filtercmd, shell=True)
	hapoutfile = reference+".hapout"
	hapcutcmd = scriptdir+"HAPCUT --fragments "+filteredfragfile+" --VCF "+vcffile+" --output "+hapoutfile+" --maxiter 11"
	subprocess.call(hapcutcmd, shell=True)
	
	return hapoutfile

def calls_to_gff(reference, hapout):
	haplotypercmd = "python "+scriptdir+"haplotyper.py "+reference+" "+hapout
	subprocess.call(haplotypercmd, shell=True)
	return 0

def output_main(args):
	"""
	main() function, files are all handled as strings of filenames
	"""
	#clean input file (fold and remove escape chars)
	reference = clean_fasta(args.infile)
	filterthreshold = args.threshold
	#look up proper readset using readset module
	readset = args.readset
	#if readset is in fasta format, inject fake quality scores
	
	#run bwa
	samfile = run_bwa(reference, readset)
	#convert sam to bam file, and sort
	sortedbam = sam_to_sorted_bam(reference, samfile)
	#run variant caller freebayes
	vcffile = run_var_caller(reference, sortedbam)
	#run hapcut suite
	hapoutfile = run_haplotyper(reference, vcffile, sortedbam, filterthreshold)
	#convert hapcut output to sequence and gff
	calls_to_gff(reference, hapoutfile)

##Read in arguments##

#input reference fasta required
#readset (this can be PE fastq, SE fastq, or fasta)
#output folder, optional (default is somesort of timestamp)
#TODO:variant caller, optional (default is freebayes)
#TODO:haplotyper, optional (default is hapCUT, not sure if there are any other realistic alternatives)
if __name__ == "__main__":
	msg = "maps reads using bwa, calls variants using freebayes, haplotypes using hapcut, and then produces haplotyped fasta and gff files"
	parser = argparse.ArgumentParser(description=msg)
	parser.add_argument("-i", "--infile",
						help="input consensus fasta file that you wish to haplotype",
						type=str, default=None, required=True)
	parser.add_argument("-r", "--readset",
						help="path to fastq file of reads",
						type=str, required=True)
	parser.add_argument("-d", "--directory",
						help="label the output directory, default is a timestamp",
						default=time.time(), required=False)
	parser.add_argument("-t", "--threshold",
						help="filter threshold, throw away reads that have fewer than this number of variants on them",
						default=1, required=False)
	args = parser.parse_args()
	output_main(args)

