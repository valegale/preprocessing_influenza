import csv
import pandas as pd
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')
import seaborn as sns; sns.set(color_codes=True)


def create_ha_na_matrices(path, bin1, bin2):
	
	distmat_path_ha = "{}/ha_alignment/emboss_ha_{}_{}.dist".format(path, bin1, bin2)
	distmat_path_na = "{}/na_alignment/emboss_na_{}_{}.dist".format(path, bin1, bin2)
	
	with open(distmat_path_ha) as file:
		tsv_file = csv.reader(file, delimiter="\t")

		for i in range(8):
			next(tsv_file)

		distance_matrix_ha = pd.DataFrame(tsv_file)
	with open(distmat_path_na) as file:
		tsv_file = csv.reader(file, delimiter="\t")

		for i in range(8):
			next(tsv_file)

		distance_matrix_na = pd.DataFrame(tsv_file)

    #from triangular matrix (output of distmat) it creates a symmetric matrix
	def upp2sym(a):
		return np.where(a,a,a.T)


	new_names = distance_matrix_ha.iloc[: , -1]
	new_names = pd.Series([name for name, index in new_names.str.split(' ', 1)])

	distance_matrix_ha = distance_matrix_ha.iloc[:,:-2]
	distance_matrix_ha = distance_matrix_ha.iloc[: , 1:]

	distance_matrix_ha = pd.DataFrame(upp2sym(distance_matrix_ha)).apply(pd.to_numeric)
	distance_matrix_ha.rename(columns=new_names, index= new_names, inplace=True)


	new_names = distance_matrix_na.iloc[: , -1]
	new_names = pd.Series([name for name, index in new_names.str.split(' ', 1)])

	distance_matrix_na = distance_matrix_na.iloc[:,:-2]
	distance_matrix_na = distance_matrix_na.iloc[: , 1:]

	distance_matrix_na = pd.DataFrame(upp2sym(distance_matrix_na)).apply(pd.to_numeric)
	distance_matrix_na.rename(columns=new_names, index= new_names, inplace=True)

	name_strain_b1 = pd.read_csv('{}/metadata/all_cleaned_bins/{}.csv'.format(path,bin1))["name_strain"]
	name_strain_b2 = pd.read_csv('{}/metadata/all_cleaned_bins/{}.csv'.format(path,bin2))["name_strain"]
	name_strain_b1 = name_strain_b1.str.replace('/','_')
	name_strain_b2 = name_strain_b2.str.replace('/','_')

	between_distance_matrix_ha = distance_matrix_ha.loc[name_strain_b1]
	between_distance_matrix_ha = between_distance_matrix_ha[name_strain_b2]

	between_distance_matrix_na = distance_matrix_na.loc[name_strain_b1]
	between_distance_matrix_na = between_distance_matrix_na[name_strain_b2]

	return between_distance_matrix_ha, between_distance_matrix_na, name_strain_b1, name_strain_b2


def compute_outliers(bin1, bin2, path, plot = False):
	# The emboss matrices are supposed to be inside the {protein}_alignment folders.
	
	
	between_distance_matrix_ha, between_distance_matrix_na, name_strain_b1, name_strain_b2 = create_ha_na_matrices(path, bin1, bin2)
	
	
	rows_sum = np.mean(between_distance_matrix_ha, axis=0)
	cols_sum = np.mean(between_distance_matrix_ha, axis=1)
	
	if plot:
		# matplotlib histogram
		plt.hist(cols_sum, color = 'blue', edgecolor = 'black',
		         bins = 1000)
		
		# Add labels
		plt.title('Density plot bin1')
		plt.xlabel('average distance')
		plt.ylabel('# sequences')
		
		plt.hist(rows_sum, color = 'blue', edgecolor = 'black',
		         bins = 1000)
		
		# Add labels
		plt.title('Density plot bin2')
		plt.xlabel('average distance')
		plt.ylabel('# sequences')
		
	outliers_first_bin_ha = set(cols_sum[cols_sum > 5].index)
	outliers_second_bin_ha = set(rows_sum[rows_sum > 5].index)
	
	rows_sum = np.mean(between_distance_matrix_na, axis=0)
	cols_sum = np.mean(between_distance_matrix_na, axis=1)
	
	if plot:
		# matplotlib histogram
		plt.hist(cols_sum, color = 'blue', edgecolor = 'black',
		         bins = 1000)
		
		# Add labels
		plt.title('Density plot bin1')
		plt.xlabel('average distance')
		plt.ylabel('# sequences')

		plt.hist(rows_sum, color = 'blue', edgecolor = 'black',
		         bins = 1000)
		
		# Add labels
		plt.title('Density plot bin2')
		plt.xlabel('average distance')
		plt.ylabel('# sequences')
	
	outliers_first_bin_na = set(cols_sum[cols_sum > 5].index)
	outliers_second_bin_na = set(rows_sum[rows_sum > 5].index)
	print (outliers_first_bin_ha)
	print (outliers_first_bin_na)
	print (outliers_second_bin_ha)
	print (outliers_second_bin_na)
	outliers_first_bin = set.union(outliers_first_bin_ha, outliers_first_bin_na)
	outliers_second_bin = set.union(outliers_second_bin_ha, outliers_second_bin_na)
	return outliers_first_bin, outliers_second_bin