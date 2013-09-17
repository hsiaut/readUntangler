#reads in output from HAPCUT, matches it to the reference file
# then outputs the sequences for the two haplotypes
# if more than one block is found, it arbitrarily chooses a block for each position
# and notes that a break was found

from pprint import pprint
import sys


VERBOSE = False

def parse_fasta(filename):
    f = open(filename)
    sequences = []
    for line in f:
        if line.startswith('>'):
            name = line[1:].rstrip('\n')
        else:
            sequences.append(line.strip())
    seqonly = "".join(sequences)
    return name, seqonly

def parse_hapcutout(filename):
	f = open(filename)
	haplotype1 = {}
	haplotype2 = {}
	haplotype_lists = [haplotype1, haplotype2]
	haplotype_ref =[{},{}]
	blocks_info = []
	blocks_flag = False
	for line in f:
		if line.startswith("BLOCK"):
			blocks_flag = True #this is used in the else statement to add the coordinates
			block={}
			splitline = line.split(" ")
			#example header line is
			# BLOCK: offset: 10187 len: 1644 phased: 1598 SPAN: 50008 MECscore 2974.57 fragments 602
			block["length"]=int( splitline[8] )
			blocks_info.append(block)
		elif line.startswith("*"):
			next
		else:
			#haplotype.out field definitions, the 1 and 0 are used to denote which column belongs to which haplotype
			#the allele in the left column is the one found in the reference sequence
			fields = line.split("\t")
			loc = int(fields[4])
			assignment1=int(fields[1])
			assignment2=int(fields[2])
			hap1 = fields[5]
			hap2 = fields[6]
			haplotype_lists[assignment1][loc]=hap1
			haplotype_ref[assignment1][loc]= True #reference sequence allele
			haplotype_lists[assignment2][loc]=hap2
			haplotype_ref[assignment2][loc]=False #not ref
			#if previous line was a block definition, we should use the coordinates of this variation as the start of the haplotype block
			if(blocks_flag == True):
				blocks_flag = False
				blocks_info[-1]["start"] = loc #exception, block not defined before variants
				blocks_info[-1]["stop"] = loc + blocks_info[-1]["length"]
				if VERBOSE:
					print blocks_info[-1]
	return haplotype_lists, haplotype_ref, blocks_info



referencefile=sys.argv[1]
haplotypefile=sys.argv[2]
haplotypes, haplotype_ref, blocks_info = parse_hapcutout(haplotypefile)
name, reference = parse_fasta(referencefile)

def make_haplotype_sequence(hapseq, haplotypes, other_haplotypes, haplotype_ref, blocks_info, one_or_two, name):
	haplist = list(hapseq)
	for k,v in haplotypes.items():
		if VERBOSE:
			print k, "replacing ", haplist[k-1], ' with ', v
		if(len(v) == 1 and len(other_haplotypes[k]) == 1):
			#pure SNP, we can simply substitute
			haplist[k-1] = v
		elif (len(v) == len(other_haplotypes[k])):
			#multi locus variant, but same length
			for i in range(0, len(v)):
				idx = i + k-1
				haplist[idx] = v[i]
		else:
			if VERBOSE:
				print "indel"
			#indel, we need to be more careful
			reference = haplotype_ref[k]
			start = k-1
			stop = k-1 + len(v)
			if ( reference ):
				#do nothing, allele already in reference
				a = 1+1
			else:
				#change to the variant seq
				ref_stop = start + len(other_haplotypes[k])
				if VERBOSE:
					print "replacing bases ", hapseq[start:ref_stop], " with ", v
				for iterator in range(start,ref_stop):
					haplist[iterator] = ''
				haplist[start] = v


	#make the gff string
	#next, we want to iterate over the dict again, this time we would like to compute the new coordinates
	## for each variant
	sortedkeylist = haplotypes.keys()
	sortedkeylist.sort()
	gff = ''
	one_or_two = str(one_or_two)
	for key in sortedkeylist:
		#key points to an entry in haplist that is a variant position, to get the new coordinates
		# we must find the length of the slice of haplist that comes before the variant position
		new_coord = len(''.join(haplist[:key]))
		length_of_var = len( haplotypes[key] )
		stop_coord = new_coord + length_of_var -1
		alleles = haplotypes[key]+"/"+other_haplotypes[key]
		#gff string
		gff += name+'-'+one_or_two+'\treadUntangler\t'+alleles+'\t'+str(new_coord)+'\t'+str(stop_coord)+'\t.\t.\t.\t.\n'
	block_num = 1
	for block in blocks_info:
		#add on haplotype block information
		newstart = len(''.join(haplist[:block["start"]]))
		newstop = len(''.join(haplist[:block["stop"]]))
		gff += name+'-'+one_or_two+'\treadUntangler\thaplotype block'+str(block_num)+'\t'+str(newstart)+'\t'+str(newstop)+'\t.\t.\t.\t.\n'   #score . . gffattribute
		block_num += 1
	hapseq = ''.join(haplist)
	return (hapseq, gff)

print "haplotype a"
(hapseq1, gff1) = make_haplotype_sequence(reference, haplotypes[1], haplotypes[0], haplotype_ref[1], blocks_info, 1, name)
print "==============haplotype b============"
(hapseq2, gff2) = make_haplotype_sequence(reference, haplotypes[0], haplotypes[1], haplotype_ref[0], blocks_info, 2, name)



f1=open(referencefile+".haplotype_1.fasta","w+")
f2=open(referencefile+".haplotype_2.fasta","w+")
gffile1=open(referencefile+".haplotype_1.gff", "w+")
gffile2=open(referencefile+".haplotype_2.gff", "w+")
print >>f1, ">"+name+"-1"
print >>f1, hapseq1
print >>f2, ">"+name+"-2"
print >>f2, hapseq2
print >>gffile1, gff1
print >>gffile2, gff2
f1.close()
f2.close()
gffile1.close()
gffile2.close()
