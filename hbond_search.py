import os
import sys
import numpy as np
import pandas as pd
from ccdc import io
from ccdc.molecule import Molecule

"""
hbond_search

takes in a CSV file of results from a substructure search, using the identifiers
to obtain crystal structures and hydrogen bonding data

the input CSV must at least have data with columns Identifier, TOR0, TOR1, N1, N2

"""

def hbond_csv_search(CSV_arg,name):
	filename = CSV_arg
	#add trycatch for wrong columns
	data = pd.read_csv(filename, header=0, index_col=False, usecols=['Identifier','TOR0','TOR1', 'N1', 'N2'])
	
	crystals = []
	crystalhbonds = []
	crystalReader = io.EntryReader('CSD')
	#initial search
	print('Obtaining hydrogen bond data...')
	for i in range(len(data)):
		crystals.append(crystalReader.crystal(data['Identifier'][i]))
		crystalhbonds.append(crystals[i].hbonds(path_length_range=(-1,999)))
		
	#search for hbonds on crystal structures with missing hydrogens
	criterion = Molecule.HBondCriterion()
	for i in range(len(crystalhbonds)):
		if crystalhbonds[i] == ():
			crystalhbonds[i] = crystals[i].hbonds(path_length_range=(-1,999), require_hydrogens=False, hbond_criterion = criterion)
	
	#count the number of hbonds a substructure is making
	num_urea_hbonds = np.zeros(len(crystals), dtype=int)
	for i in range(len(crystalhbonds)):
		for j in range(len(crystalhbonds[i])):
			for k in range(len(crystalhbonds[i][j].atoms)):
				if data['N1'][i] == crystalhbonds[i][j].atoms[k].label or data['N2'][i] == crystalhbonds[i][j].atoms[k].label:
					num_urea_hbonds[i] +=1
	data['# Hbonds'] = num_urea_hbonds
	data['Hbonds'] = crystalhbonds
	filename = os.getcwd() + '/hbonds_' + name + '.CSV'
	data.to_csv(filename)
	print('File saved to: ' + filename) 
	
def hbond_df_search(DF):
	data = DF
	
	crystals = []
	crystalhbonds = []
	crystalReader = io.EntryReader('CSD')

	#initial search
	print('Obtaining hydrogen bond data...')
	for i in range(len(data)):
		crystals.append(crystalReader.crystal(data['Identifier'][i]))
		crystalhbonds.append(crystals[i].hbonds(path_length_range=(-1,999)))
		
	#search for hbonds on crystal structures with missing hydrogens
	criterion = Molecule.HBondCriterion()
	for i in range(len(crystalhbonds)):
		if crystalhbonds[i] == ():
			crystalhbonds[i] = crystals[i].hbonds(path_length_range=(-1,999), require_hydrogens=False, hbond_criterion = criterion)
	
	strength = [[] for i in range(len(crystals))]
	#count the number of hbonds a substructure is making
	num_urea_hbonds = np.zeros(len(crystals), dtype=int)
	for i in range(len(crystalhbonds)):
		for j in range(len(crystalhbonds[i])):
			for k in range(len(crystalhbonds[i][j].atoms)):
				if data['N1'][i] == crystalhbonds[i][j].atoms[k].label or data['N2'][i] == crystalhbonds[i][j].atoms[k].label:
					num_urea_hbonds[i] +=1
					strength[i].append(crystalhbonds[i][j].strength)

	data['# Hbonds'] = num_urea_hbonds
	data['Strength'] = strength
	data['Hbonds'] = crystalhbonds
	
	return data
