import pandas as pd


def translate_sequences(msa, path_metadata, bin1, bin2, protein, method = 1):
	# This function translates all the sequences that are part of an msa. 

	for col in range(msa.get_alignment_length()): #trim at ATG
	
		if (method == 1):
			if msa[:, col] == "A"*len(msa) and msa[:, col + 1] == "T"*len(msa) and msa[:, col + 2] == "G"*len(msa):
				position_start = col
				print("First full column is {}".format(col))
				break 

		elif (method == 2):
	
			if msa[2, col] == "A" and msa[2, col + 1] == "T" and msa[2, col + 2] == "G": # sometimes some sequences start two nucleotide later!
				position_start = col
				print("First full column is {}".format(col))
				break
	
	
	for col in reversed(range(msa.get_alignment_length())): #trim at the last full column
		
		if (method == 1):
		    if not "-" in msa[:, col]:
		        position_end = col
		        print("Last full column is {}".format(col))
		        break
		
		elif (method == 2): #inserting a stop codon check and one sequence here
			if not "-" in msa[1, col] and msa[1, col] == "A" and msa[1, col + 1] == "T" and msa[2, col + 2] == "T":
				position_end = col
				print("Last full column is {}".format(col))
				break
	

	trimmed_msa = msa[:, position_start:position_end - 2] # trimming between ATG and removing the stop codon (position_end - 2)
	

	#sequence = trimmed_msa[2].seq.ungap("-")
	#print (sequence.translate())
	#return 
	name_strain_b1 = list(pd.read_csv('{}/metadata/all_cleaned_bins/{}.csv'.format(path_metadata,bin1))["name_strain"])
	name_strain_b2 = list(pd.read_csv('{}/metadata/all_cleaned_bins/{}.csv'.format(path_metadata,bin2))["name_strain"])
	
	bin1_seq = []
	bin2_seq = []
	
	for seq in trimmed_msa:
		
		if protein == "na" and len(seq.seq.ungap("-")) == 1407:
			seq.seq = seq.seq.ungap("-")
		elif protein == "ha" and len(seq.seq.ungap("-")) == 1698:
			seq.seq = seq.seq.ungap("-")
		
		if seq.id in name_strain_b1:
			try:
				seq.seq = seq.seq.translate()
			except:
				print("An exception occurred bin 1")
				print(seq.id)
				#return 
			bin1_seq.append(seq)
			
			#len check
			if protein == "ha":
				if len(seq.seq) != 566:
					print(seq.id + " from bin " + bin1 + " has length " + str(len(seq.seq)))
			if protein == "na":
				if len(seq.seq) != 469:
					print(seq.id + " from bin " + bin1 + " has length " + str(len(seq.seq)))
			
		elif seq.id in name_strain_b2:
			try:
				seq.seq = seq.seq.translate()
			except:
				print("An exception occurred bin 2 in " + seq.id)
			bin2_seq.append(seq)
			
			#len check
			if protein == "ha":
				if len(seq.seq) != 566:
					print(seq.id + " from bin " + bin2 + " has length " + str(len(seq.seq)))
			if protein == "na":
				if len(seq.seq) != 469: 
					print(seq.id + " from bin " + bin2 + " has length " + str(len(seq.seq)))
		
				
	return bin1_seq, bin2_seq
	
		
	
	
	
