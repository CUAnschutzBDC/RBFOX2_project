import sys
from Bio import SeqIO
import gzip
import os
import argparse

### SETUP ###
# House keeping to read in arguments from the command line
parser = argparse.ArgumentParser()

parser.add_argument("-g", "--genome", dest = "genomefasta",
    help = "Path to the genome fasta file.",
    default = "/beevol/home/mossnico/GEO_Datasets/ref/star/human/GRCh38/GRCh38.primary_assembly.genome.fa",
    action = "store",
    metavar = "\b")

parser.add_argument("-i", "--input_dir", dest = "input_dir",
    help = "Path to the input directory. Default is /beevol/home/mossnico/.",
    default = "/beevol/home/mossnico/",
    action = "store",
    metavar = "\b")

options = parser.parse_args()



### INPUT FILES ###


# exon coordiante files
inclusion_coordinates = os.path.join(options.input_dir, 'Human_T2D_inclusionexon.coordinates.txt')
exclusion_coordinates = os.path.join(options.input_dir, 'Human_T2D_exclusionexon.coordinates.txt')
insensitive_coordinates = os.path.join(options.input_dir, 'Human_T2D_insensitive.coordinates.txt')

# sequence files
inc_upstream_fasta = os.path.join(options.input_dir, 'inclusionexons.upstream.fa')
inc_downstream_fasta = os.path.join(options.input_dir, 'inclusionexons.downstream.fa')
exc_upstream_fasta = os.path.join(options.input_dir, 'exclusionexons.upstream.fa')
exc_downstream_fasta = os.path.join(options.input_dir, 'exclusionexons.downstream.fa')
insensitive_upstream_fasta  = os.path.join(options.input_dir,'insensitive.upstream.fa')
insensitive_downstream_fasta  = os.path.join(options.input_dir,'insensitive.downstream.fa')


### FUNCTIONS ###

seq_dict = SeqIO.to_dict(SeqIO.parse(open(options.genomefasta, 'rt'), 'fasta'))

# translate coordinates to sequence using - cordstoseq(exon coordinate file, upstream sequence, downstream sequence)
def coordstoseq(coords, upstreamseqfile, downstreamseqfile):
  with open(coords, 'r') as coordfh, open(upstreamseqfile, 'w') as upstreamoutfh, open(downstreamseqfile, 'w') as downstreamoutfh:
      for line in coordfh:
        #Remove trailing newline characters and turn each line into a list
        line = line.strip().split('\t')
        chrm = line[1]
        strand = line[2]
        exonstart = int(line[3])
        exonstop = int(line[4])
        seqid = ('_').join(line)
          
        #If this is a positive strand gene, things are pretty straightfoward
        if strand == '+':
          upstreamintstart = exonstart - 200
          upstreamintend = exonstart
          downstreamintstart = exonstop
          downstreamintend = exonstop + 200
          upstreamseq = seq_dict[chrm].seq[upstreamintstart : upstreamintend].transcribe()
          downstreamseq = seq_dict[chrm].seq[downstreamintstart : downstreamintend].transcribe()
        #If it's a negative strand gene the upstream intron (in the RNA sense) is actually after this intron (in the genome sense)
        #AND we need to take the reverse complement
        elif strand == '-':
          upstreamintstart = exonstop
          upstreamintend = exonstop + 200
          downstreamintstart = exonstart - 200
          downstreamintend = exonstart
          upstreamseq = seq_dict[chrm].seq[upstreamintstart : upstreamintend].reverse_complement().transcribe() #biopython has a reverse complement function
          downstreamseq = seq_dict[chrm].seq[downstreamintstart : downstreamintend].reverse_complement().transcribe()
          
        #Write sequence in fasta format
        upstreamoutfh.write ('>' + seqid + '\n' + str(upstreamseq) + '\n')
        downstreamoutfh.write ('>' + seqid + '\n' + str(downstreamseq) + '\n')
        
# generate list of kmers
from itertools import product

# make a list of all possible 5mers 
nucleotides = ['A', 'G', 'C', 'U']
allkmers = [''.join(kmer) for kmer in product(nucleotides, repeat=5)]

# build a list of kmers and record each observation of the kmer along the sequence using - countkmers(input fasta, outfile, k is nt in kmer)
def countkmers(fasta, outfile, k):
  k = int(k)
  kmercounts = {} # generate a disctionary {kmer : number of times we observe that kmer}
  with open(fasta, 'r') as fastafh:
    # loop through every sequence in the fasta file
    for record in SeqIO.parse(fastafh, 'fasta'):
      seq = str(record.seq)
      # Count kmers
      for i in range(len(seq) - k + 1):
        kmer = seq[i : i+k]
        # If we haven't seen this kmer before, its count is 1
        if kmer not in kmercounts:
          kmercounts[kmer] = 1
        # If we have seen this kmer before, add one to its count
        elif kmer in kmercounts:
          kmercounts[kmer] +=1
          
  # dictionary contains {kmer : number of times we observe that kmer}, add all other possible kmers to the dictionary with a count of 0
  for kmer in allkmers:
    if kmer not in kmercounts:
      kmercounts[kmer] = 0
      
  # ouptup file will contain each kmer, the counts for that kmer, and counts of all other kmers
  with open(outfile, 'w') as outfh:
    totalkmercounts = sum(kmercounts.values())
    outfh.write('kmer' + '\t' + 'kmercount' + '\t' + 'otherkmerscount' + '\n')
    for kmer in kmercounts:
      kmercount = kmercounts[kmer]
      otherkmerscount = totalkmercounts - kmercount
      outfh.write(kmer + '\t' + str(kmercount) + '\t' + str(otherkmerscount) + '\n')


### RUN FUNCTIONS ###

# run cordstoseq() to get the sequences for each exon type, inclusion, exclusion, and insensitive 

coordstoseq(inclusion_coordinates,inc_upstream_fasta, inc_downstream_fasta)
coordstoseq(exclusion_coordinates,exc_upstream_fasta ,exc_downstream_fasta)
coordstoseq(insensitive_coordinates,insensitive_upstream_fasta, insensitive_downstream_fasta)

# evaluate kmers for each set of sequenes      
countkmers(inc_upstream_fasta, 'inclusionexons.upstream.kmercounts.txt', 5)
countkmers(inc_downstream_fasta, 'inclusionexons.downstream.kmercounts.txt', 5)      
countkmers(exc_upstream_fasta, 'exclusionexons.upstream.kmercounts.txt', 5)
countkmers(exc_downstream_fasta, 'exclusionexons.downstream.kmercounts.txt', 5)
countkmers(insensitive_upstream_fasta, 'insensitives.upstream.kmercounts.txt', 5)
countkmers(insensitive_downstream_fasta, 'insensitives.downstream.kmercounts.txt', 5)
