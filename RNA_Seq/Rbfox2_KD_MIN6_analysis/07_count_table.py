from collections import defaultdict
import re

file_list = ['MIN6_RBFOX2_KD/counts/NT1_S5_counts.txt','MIN6_RBFOX2_KD/counts/NT3_S6_counts.txt','MIN6_RBFOX2_KD/counts/NT5_S7_counts.txt','MIN6_RBFOX2_KD/counts/NT6_S8_counts.txt','MIN6_RBFOX2_KD/counts/RbFox2_1_S1_counts.txt','MIN6_RBFOX2_KD/counts/RbFox2_3_S2_counts.txt','MIN6_RBFOX2_KD/counts/RbFox2_5_S3_counts.txt','MIN6_RBFOX2_KD/counts/RbFox2_6_S4_counts.txt'] 
output_file = 'MIN6_RBFOX2_KD/Rbfox2_MIN6_count.txt' 

gene_dict = defaultdict(list)
sample_list = list()

for i in file_list:
	sample_name = re.sub(r'_counts.txt', '', i)
	sample_list.append(sample_name)
	with open(i, "r") as countFile:
		for line in countFile:
			line = line.strip().split("\t")
			if "ENSMUSG" in line[0]:
				ens_id = line[0]
				gene = line[6]
				count = line[8]
				gene_ens = gene + "_" + ens_id
				gene_dict[gene_ens].append(count)

with open(output_file, "w") as count_file:
	count_file.write("gene" + "\t" + "\t".join(sample_list) + "\n")
	for gene in gene_dict:
		count_file.write(gene + "\t" + "\t".join(gene_dict[gene]) + "\n")
