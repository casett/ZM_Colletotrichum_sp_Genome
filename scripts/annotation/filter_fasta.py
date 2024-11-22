# Python
# To get secreted proteins 
# By Cassie

import sys
import Bio
from Bio import SeqIO
import argparse 

parser = argparse.ArgumentParser(description='Process to obtain secreted proteins (or proteins of interest) based on gene ID. Created to get secreted proteins for use with Effector3-P.')
parser.add_argument('-i', '--input', help='A fasta file with multiple sequences')
parser.add_argument('-g', '--genes', help='A file of gene IDs')
parser.add_argument('-o', '--out', help='A fasta file name to save output too')

	
args = parser.parse_args()

def filter_fasta(input_fasta, output_file, genes): 
											
	seq_records = SeqIO.parse(input_fasta, format='fasta') #parses the fasta file
	
	with open(genes) as f:
			gene_ids_list = f.read().splitlines() #parse the contamination file which is each line as a scaffold id 
	
	OutputFile = open(output_file, 'w') #opens new file to write to
	
	for record in seq_records: 
		#print(record.id)
		if record.id.split('-',1)[0] in gene_ids_list:
			OutputFile.write('>'+ record.id +'\n') #writes the scaffold to the file (or assession) 
			OutputFile.write(str(record.seq)+'\n') #writes the seq to the file
	
	OutputFile.close()


#Run

filter_fasta(args.input, args.out, args.genes)
