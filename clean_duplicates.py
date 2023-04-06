from Bio import SeqIO
import random
#script to clean the bins from duplicates (both na and ha are the same)

def remove_dups(path, bin1, len_dictionary):

	ha_protein_file = "{}/ha_proteins/ha_{}.fa".format(path,bin1)
	na_protein_file = "{}/na_proteins/na_{}.fa".format(path, bin1)
		
	#concatenating the strings
	concat_dictionary = {}
	for seq_record in SeqIO.parse(ha_protein_file, "fasta"):
		concat_dictionary[seq_record.id] = str(seq_record.seq)
		
	for seq_record in SeqIO.parse(na_protein_file, "fasta"):
		concat_dictionary[seq_record.id] += str(seq_record.seq)
	
	
	#finding the duplicates
	rev_multidict = {}
	for key, value in concat_dictionary.items():
		rev_multidict.setdefault(value, set()).add(key)


	cleaning = [values for key, values in rev_multidict.items() if len(values) > 1]
	
	unique_sequences = [random.choice(list(elem)) for elem in cleaning]
	#print([len(seq) for seq in cleaning])
	final_data = [values.pop() for key, values in rev_multidict.items() if len(values) == 1] + unique_sequences
	len_dictionary[bin1] = {}
	len_dictionary[bin1]["complete"] = len(concat_dictionary)
	len_dictionary[bin1]["unique"] = len(final_data)
	
	ha_proteins = []
	na_proteins = []
	for seq_record in SeqIO.parse(ha_protein_file, "fasta"):
		if seq_record.id in final_data:
			ha_proteins.append(seq_record)
	
	for seq_record in SeqIO.parse(na_protein_file, "fasta"):
		if seq_record.id in final_data:
			na_proteins.append(seq_record)
	
	return ha_proteins, na_proteins, len_dictionary
	#SeqIO.write(ha_proteins, ha_protein_file_new , "fasta")
	#SeqIO.write(na_proteins, na_protein_file_new , "fasta")
	
"""
with open("cleaned_sequences.csv", 'w', newline='') as f:
		
	writer = csv.writer(f)
	header = ["bin", "#complete", "#unique"]	
	writer.writerow(header)
	data = []
		
	for bins in len_dictionary:
		data.append([bins, len_dictionary[bins]["complete"], len_dictionary[bins]["unique"]])
	writer.writerows(data)
"""
	