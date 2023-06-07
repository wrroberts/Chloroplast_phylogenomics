#! /usr/bin/env python3

from Bio import SeqIO
import argparse
import os,sys

# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
parser = argparse.ArgumentParser(description = "extract features and create fasta files from a directory contining GenBank files" )

# add arguments
parser.add_argument( "ending", help="file ending for summary files" )
parser.add_argument( "directory", help="working directory" )

args = parser.parse_args()

DIR = args.directory + "/"

for i in os.listdir(DIR):
	if i[-len(args.ending):] == args.ending in i:
		if os.path.getsize(i) > 1:
			taxon = i.split('.')[0]
			print(i)
			
			for record in SeqIO.parse(i, 'genbank'):
				seq = str(record.seq)
				for feature in record.features:
					if feature.type == 'CDS' or feature.type == 'rRNA' or feature.type == 'tRNA':
						if 'gene' in feature.qualifiers:
							geneName = feature.qualifiers['gene'][0]
						elif 'product' in feature.qualifiers:
							geneName = feature.qualifiers['product'][0]
						else:
							print('ERROR when parsing feature:')
							print(feature.qualifiers)
							quit()
							
#						if seq_feature.location.strand == 1:
#                          	complement = 'no complement'
#                       if seq_feature.location.strand == -1:
#                        	complement = str(seq_feature.location).reverse_complement()
							
						protFile = open(geneName + '.faa', 'a')
						protFile.write('>')
						protFile.write(taxon + '_' + geneName + '\n')
						if feature.strand == -1:
							protFile.write(str(record.seq[feature.location.start.position: feature.location.end.position].reverse_complement().translate(table=11)))
							protFile.write('\n\n')
						else:
							protFile.write(str(record.seq[feature.location.start.position: feature.location.end.position].translate(table=11)))
							protFile.write('\n\n')
							
							
						geneFile = open(geneName + '.fna', 'a')
#						if os.path.exists(str(geneFile)):
#							with open(geneFile, 'a') as myfile:
#								geneFile.write('>')
#								geneFile.write(taxon + '_' + geneName + '\n')
#								
#								geneFile.write(seq[feature.location.start.position: feature.location.end.position])
#								geneFile.write('\n\n')

						geneFile.write('>')
						geneFile.write(taxon + '_' + geneName + '\n')
						if feature.strand == -1:
							geneFile.write(str(record.seq[feature.location.start.position: feature.location.end.position].reverse_complement()))
							geneFile.write('\n\n')
						else:
							geneFile.write(str(record.seq[feature.location.start.position: feature.location.end.position]))
							geneFile.write('\n\n')
						
			
			
