from Bio import SeqIO, AlignIO
from translate_sequences import translate_sequences
from clean_reassortant import compute_outliers
from clean_duplicates import remove_dups
import os
import csv

first_bin = "10_11"
last_bin = "19_20"
path = "../emboss/H1N1"
protein = "ha"


#options: remove_reassortant, translate, remove_duplicates, full_pipeline (not suggested)
option = "remove_duplicates"

def padding_zero(i):
	return str(i).zfill(2)

if option == "remove_reassortant":
	# with this options new files are created into a new folder called "{protein}_filtered"
	# both HA and NA proteins automatically
	
	# create {protein}_filtered if it doesn't exist
	path_ha_filtered = "ha_filtered"
	path_na_filtered = "na_filtered"
		
	if not os.path.exists(os.path.join(path, path_ha_filtered)):
		os.makedirs(os.path.join(path, path_ha_filtered))
			
	if not os.path.exists(os.path.join(path, path_na_filtered)):
		os.makedirs(os.path.join(path, path_na_filtered))
	
	for i in range(int(first_bin[0:2]), int(last_bin[0:2])):
		bin1 = "{}_{}".format(padding_zero(i), padding_zero(i+1))
		bin2 = "{}_{}".format(padding_zero(i+1), padding_zero(i+2))

			
		outliers_first_bin, outliers_second_bin = compute_outliers(bin1, bin2, path)
		
		#print (outliers_second_bin)
		if bin1 == "10_11" or bin1 == "15_16": # bin 14_15 really small, had to compare 15_16 with the following one
		#if True: # bin 14_15 really small, had to compare 15_16 with the following one
			
			ha_filtered = []
			na_filtered = []
			
			ha_cleaned_file = "{}/ha_cleaned/ha_{}.fa".format(path, bin1)
			na_cleaned_file = "{}/na_cleaned/na_{}.fa".format(path, bin1)
	
			ha_filtered_file = "{}/ha_filtered/ha_{}.fa".format(path, bin1)
			na_filtered_file = "{}/na_filtered/na_{}.fa".format(path, bin1)
			
			
			for seq_record in SeqIO.parse(ha_cleaned_file, "fasta"):
				if seq_record.id.replace('/','_') not in outliers_first_bin:
					 ha_filtered.append(seq_record)
	
			for seq_record in SeqIO.parse(na_cleaned_file, "fasta"):
				if seq_record.id.replace('/','_') not in outliers_first_bin:
					na_filtered.append(seq_record)	
	
			SeqIO.write(ha_filtered, ha_filtered_file , "fasta")
			SeqIO.write(na_filtered, na_filtered_file , "fasta")
			
		ha_filtered = []
		na_filtered = []
			
		ha_cleaned_file = "{}/ha_cleaned/ha_{}.fa".format(path, bin2)
		na_cleaned_file = "{}/na_cleaned/na_{}.fa".format(path, bin2)
	
		ha_filtered_file = "{}/ha_filtered/ha_{}.fa".format(path, bin2)
		na_filtered_file = "{}/na_filtered/na_{}.fa".format(path, bin2)
			
			
		for seq_record in SeqIO.parse(ha_cleaned_file, "fasta"):
			if seq_record.id.replace('/','_') not in outliers_second_bin:
				ha_filtered.append(seq_record)
	
		for seq_record in SeqIO.parse(na_cleaned_file, "fasta"):
			if seq_record.id.replace('/','_') not in outliers_second_bin:
				na_filtered.append(seq_record)	
		
		
		#SeqIO.write(ha_filtered, ha_filtered_file , "fasta")
		#SeqIO.write(na_filtered, na_filtered_file , "fasta")
		

elif option == "translate":
	# Translate a file with multiple sequence alignment obtained with clustalo and stored inside
	# {protein}_alignment. 

	# create {protein}_filtered if it doesn't exist
	path_ha_proteins = "ha_proteins"
	path_na_proteins = "na_proteins"
		
	if not os.path.exists(os.path.join(path, path_ha_proteins)):
		os.makedirs(os.path.join(path, path_ha_proteins))
			
	if not os.path.exists(os.path.join(path, path_na_proteins)):
		os.makedirs(os.path.join(path, path_na_proteins))
	
	for i in range(int(first_bin[0:2]), int(last_bin[0:2])): #but better to do one at the time!!!
		bin1 = "{}_{}".format(padding_zero(i), padding_zero(i+1))
		bin2 = "{}_{}".format(padding_zero(i+1), padding_zero(i+2))

		alignments = list(AlignIO.parse("{}/{}_alignment/{}_{}_{}_align.fa".format(path, protein,protein, bin1, bin2), "fasta"))
		
		bin1_seq, bin2_seq = translate_sequences(alignments[0], path, bin1, bin2, protein, method = 1)
		
		if bin1 == "10_11":
			SeqIO.write(bin1_seq, "{}/{}_proteins/{}_{}.fa".format(path, protein, protein, bin1), "fasta")

		#SeqIO.write(bin2_seq, "{}/{}_proteins/{}_{}.fa".format(path, protein, protein, bin2), "fasta")
		
elif option == "remove_duplicates":
	
	# create removed_duplicatesif it doesn't exist
	
	if not os.path.exists(os.path.join(path, "removed_duplicates")):
		os.makedirs(os.path.join(path, "removed_duplicates"))
	
	new_path = os.path.join(path, "removed_duplicates")
	# create {protein}_filtered inside removed duplicates if it doesn't exist
	path_ha_proteins = "ha_proteins"
	path_na_proteins = "na_proteins"
			
	if not os.path.exists(os.path.join(new_path, path_ha_proteins)):
		os.makedirs(os.path.join(new_path, path_ha_proteins))
				
	if not os.path.exists(os.path.join(new_path, path_na_proteins)):
		os.makedirs(os.path.join(new_path, path_na_proteins))
	
	len_dictionary = {}
	
	for i in range(int(first_bin[0:2]), int(last_bin[0:2]) + 1):
		bin1 = "{}_{}".format(padding_zero(i), padding_zero(i+1))
		
		
		ha_proteins, na_proteins, len_dictionary = remove_dups(path, bin1, len_dictionary)
		
		ha_protein_file_new = "{}/removed_duplicates/ha_proteins/ha_{}.fa".format(path,bin1)
		na_protein_file_new = "{}/removed_duplicates/na_proteins/na_{}.fa".format(path,bin1)
		#SeqIO.write(ha_proteins, ha_protein_file_new, "fasta")
		#SeqIO.write(na_proteins, na_protein_file_new, "fasta")
	
	with open("{}/metadata/dups.csv".format(path), 'w', newline='') as f:
			
		writer = csv.writer(f)
		header = ["bin", "#complete", "#unique"]	
		writer.writerow(header)
		data = []
			
		for bins in len_dictionary:
			data.append([bins, len_dictionary[bins]["complete"], len_dictionary[bins]["unique"]])
		writer.writerows(data)
	
	
	
