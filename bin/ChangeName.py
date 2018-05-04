import sys

def ChangeID(input_file,output_file):
	'''
	Change the gene_id with reference_gene_id, transcript_id with reference_transcript_id 
	if the entry has a reference in the gtf file 
	'''
	with open(input_file) as f1, open(output_file,'w') as f2:
		for line in f1.readlines():
			line= line.strip()
			if not line.startswith('#'):
				line = line.split('\t')
				attribute = line[8].strip().split(';')
				reference_gene_id,reference_transcript_id = False,False
				for i in attribute:
					if i.startswith(' reference_id'):
						reference_transcript_id = i.split()[1]
					elif i.startswith(' ref_gene_id'):
						reference_gene_id = i.split()[1]
				if reference_gene_id:
					attribute[0] = 'gene_id '+reference_gene_id
				if reference_transcript_id:
					attribute[1] = ' transcript_id '+reference_transcript_id
				attribute = ';'.join(attribute)
				line[8] = attribute
				line='\t'.join(line)
			f2.write(line+'\n')


def main():
	input_file = sys.argv[1]
	output_file = sys.argv[2]

	ChangeID(input_file,output_file)



if __name__ == '__main__':
	main()