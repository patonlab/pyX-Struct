from __future__ import print_function
from builtins import input
import datetime
import os
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from ccdc.search import SMARTSSubstructure, SubstructureSearch 
from ccdc.io import csd_version


#for testing
search = True 			#True: the program will perform a search on the input SMILES
graph = True			#True: the program will display angle or torsion data
searchLimit = True 		#True: the search is limited by the maxHits variable
printData = False		#True: the data will print in the terminal for the user
saveData = False		#True: Export data
maxHits = 1000			#Max number of hits until the search ends, default at 1000

#number of measurements for naming purposes d, a, t respectively
measurement = np.zeros(shape=(3,2))

def addDist(atom = None, askIndex = False):
	global substructure, substructure_search, sub_id, measurement
	name = 'D'+str(int(measurement[0][1]))
	
	#the atom indices have been specified, add a search element 
	if not askIndex: 
		substructure_search.add_distance_measurement(name,
			sub_id, atom[0], sub_id, atom[1])
	#obtain indices for distance	
	else: 
		#print molecule to user so they can see indices
		for i in substructure.atoms:
			print(i.index, i)
		#obtain values from user and add measurement search
		atom_i = input("Enter two indices to measure a distance: ")
		atom_i = list(map(int, atom_i.split()))
		substructure_search.add_distance_measurement(name,
			sub_id, atom_i[0], sub_id, atom_i[1])
	
	print('Distance', name, "added to the search.")
	measurement[0][1] += 1
				
def addAngle(atom = None, askIndex = False):
	global substructure, substructure_search, sub_id, measurement
	name = 'A'+str(int(measurement[1][1]))
	
	#the atom indices have been specified, add a search element 
	if not askIndex: 
		substructure_search.add_angle_measurement(name,
			sub_id, atom[0], sub_id, atom[1], sub_id, atom[2])
	#obtain indices for angle	
	else: 
		#print molecule to user so they can see indices
		for i in substructure.atoms:
			print(i.index, i)
		#obtain values from user and add measurement search
		atom_i = input("Enter three indices to measure an angle: ")
		atom_i = list(map(int, atom_i.split()))
		substructure_search.add_angle_measurement(name,
			sub_id, atom_i[0], sub_id, atom_i[1], sub_id, atom_i[2])
	
	print('Angle', name, "added to the search.")		
	measurement[1][1] += 1
		
def addTor(atom = None, askIndex = False):
	global substructure, substructure_search, sub_id, measurement
	name = 'TOR'+str(int(measurement[2][1]))
	
	#the atom indices have been specified, add a search element 
	if not askIndex: 
		substructure_search.add_torsion_angle_measurement(name,
			sub_id, atom[0], sub_id, atom[1], sub_id, atom[2], sub_id, atom[3])
	#obtain indices for angle
	else: 
		#print molecule to user so they can see indices
		for i in substructure.atoms:
			print(i.index, i)
		#obtain values from user and add measurement search
		atom_i = input("Enter four indices to measure a dihedral: ")
		atom_i = list(map(int, atom_i.split()))
		substructure_search.add_torsion_angle_measurement(name,
			sub_id, atom_i[0], sub_id, atom_i[1], sub_id, atom_i[2],sub_id, atom_i[3])
	
	print('Torsion', name, "added to the search.")			
	measurement[2][1] += 1

def checkMeasurement(arg):
	global graph, printData, maxHits, saveData
	#check for distance measurement
	if 'd' in arg:
		dIndices = [j for j, k in enumerate(arg) if k == 'd']
		#keep track of how many specified distances we want to measure
		measurement[0][0] = len(dIndices)
		for j in np.arange(len(dIndices)):
			atom = []
			#check if user specified atom indices in argument
			try:
				for k in np.arange(1,3):
					atom.append(int(arg[dIndices[j]+k]))
				if all(isinstance(x,int) for x in atom):
					addDist(atom=atom)
				else: 
					addDist(askIndex=True)
			except (ValueError, IndexError):
				addDist(askIndex=True)
			
					
	#check for angle measurements
	if 'a' in arg:
		aIndices = [j for j, k in enumerate(arg) if k == 'a']
		#keep track of how many specified angles we want to measure
		measurement[1][0] = len(aIndices)
		for j in np.arange(len(aIndices)):
			atom = []
			#check if user specified atom indices in argument
			try:
				for k in np.arange(1,4):
					atom.append(int(arg[aIndices[j]+k]))
				if all(isinstance(x,int) for x in atom):
					addAngle(atom=atom)
				else: 
					addAngle(askIndex=True)
			except (ValueError, IndexError):
				addAngle(askIndex=True)
						
						
	#check for torsion angle measurement
	if 't' in arg:
		tIndices = [j for j, k in enumerate(arg) if k == 't']
		#keep track of how many specified torsion angles we want to measure
		measurement[2][0] = len(tIndices)
		for j in np.arange(len(tIndices)):
			atom = []
			#check if user specified atom indices in argument
			try:
				for k in np.arange(1,5):
					atom.append(int(arg[tIndices[j]+k]))
				if all(isinstance(x,int) for x in atom):
					addTor(atom=atom)
				else:
					addTor(askIndex=True)
			except (ValueError, IndexError):
				addTor(askIndex=True)
	
	if 'lim' in arg:
		limIndex = arg.index('lim')
		lim= None
		try:
			lim = int(arg[limIndex+1])
			if isinstance(lim, int):
				maxHits = lim
			else:
				lim = int(input("Enter a number for max number of hits for the search: "))
				maxHits = lim
		except (ValueError, IndexError):
			#ask for lim
			lim = int(input("Enter a number for max number of hits for the search: "))
			maxHits = lim
			
	if 's' in arg:
		saveData = True
			
#if the user didn't pass any arguments, ask for a SMILES string
if len(sys.argv) == 1:
	smiles = input("Enter a SMILES string: ")
	#if user put quotations, remove them
	if smiles[0] == ("'" or '"'):
		smiles = smiles[1:]
	if smiles[-1] == ("'" or '"'):
		smiles = smiles[:-1]
else:
	smiles = sys.argv[1]
	
#Prepare a search
substructure = SMARTSSubstructure(smiles)
substructure_search = SubstructureSearch()
sub_id = substructure_search.add_substructure(substructure)

#determine if the user has specified other measurements
if len(sys.argv) > 2:
	arguments = sys.argv
	del arguments[0:2]
	checkMeasurement(arguments)
	#lim or save specified but measurements arent 
	if ('d' and 'a' and 't') not in sys.argv:
		#check if the user wants to search for measurements
		ans = input('Search for any specific measurements on this molecule? (y/n): ')
		if ans == 'y':
			#ask which one, run method
			ans = input('Which measurements? \nEnter d for distance, a for angle, t for torsion angle: ')
			ans = list(map(str, ans.split()))
			checkMeasurement(ans)
		else: 
			#user has no measurement parameters, dont graph anything and print the IDs
			graph = False
			printData = True
else:
	#check if the user wants to search for measurements
	ans = input('Search for any specific measurements on this molecule? (y/n): ')
	if ans == 'y':
		#ask which one, run method
		ans = input('Which measurements? \nEnter d for distance, a for angle, t for torsion angle: ')
		ans = list(map(str, ans.split()))
		checkMeasurement(ans)
	else: 
		#user has no measurement parameters, dont graph anything and print the IDs
		graph = False
		printData = True

#for testing..		
if search == False:
	exit()

#search the CDS for our specified query
if maxHits == 0:
	searchLimit = False
if searchLimit:
	print('Searching for substructures with a limit of ' + str(maxHits) + ' max structures...')
	hits = substructure_search.search(max_hit_structures=maxHits)
else:
	print('Searching for substructures...')
	hits = substructure_search.search()

#if no hits, end program
if len(hits) == 0:
		print("No substructures found.")
		exit()

#determine names of columns for the DataFrame
columns = ['Identifier']
mType = ''
for i in np.arange(3):
	for j in np.arange(measurement[i][0], dtype=int):
		if i == 0:
			mType = 'D'
		elif i == 1:
			mType = 'A'
		else:
			mType = 'TOR'
			
		columns.append(mType+str(j))

#prepare data to be placed into dataframe
ids = []
measureHits = [[] for _ in range(len(columns)-1)]
for h in hits:
	ids.append(h.identifier)
	if len(columns) > 1:
		for i in np.arange(len(columns)-1):
			measureHits[i].append(float('%.2f' % h.measurements[columns[i+1]]))

#create a dataframe and put the data into it
hitData = pd.DataFrame(columns=columns)
hitData['Identifier'] = ids
for i in np.arange(len(measureHits)):
	hitData[columns[i+1]] = np.absolute(measureHits[i])
	
#optionally print data to user
if printData:
	print(hitData.to_string())
		
#print data info to user 		
structData = str(len(hitData.Identifier)) + ' matching substructures in '+ str(len(hitData.Identifier.unique())) + ' different molecules.'
print('Found ' + structData)

#Export data 
if saveData:
	filename = os.getcwd() + '/search_' + datetime.datetime.now().strftime ("%H:%M:%S") + '.CSV'
	hitData.to_csv(filename)
	print('file saved to: ' + filename) 

#ask user which data to display 
if len(hitData.columns) == 3:
	xData = hitData.columns[1]
	yData = hitData.columns[2]
elif len(hitData.columns) > 3:
	print('Which data would you like to display? (enter two indices, # #):')
	print('i | Measurement')
	#only show angular data
	#decide title of graph
	for i in np.arange(len(columns[1:])):
		if columns[i+1][0] != 'D':
			print(str(i) + " | " + str(columns[i+1]))
	ans = input()
	ans = list(map(int, ans.split()))
	xData = hitData.columns[ans[0]+1]
	yData = hitData.columns[ans[1]+1]
	print('Graphing \'' + xData + '\' vs \'' + yData + '\'...')
elif len(hitData.columns) < 3:
	graph = False

	
if graph:
	if xData[0] == 'T' and yData[0] == 'T':
		title = 'Dihedral Angle Magnitude Comparison'
	elif xData[0] == 'A' and yData[0] == 'A':
		title = 'Angle Magnitude Comparison'
	elif xData[0] == ('T' or 'A') and yData[0] == ('T' or 'A') :
		title = 'Dihedral Angle vs. Angle Magnitude Comparison'
	else: 
		title = 'Angular Magnitude Comparison'
	#Joint plot Creator
		#gridsize: higher = smaller hexagons
		#joint_kws bins: lower = darker hexagons
		#marginal_kws bins: lower = thicker histograms
	sns.set(style='white', color_codes=True)
	g = sns.jointplot(x=xData, y=yData,data=hitData, kind="hex",color='b', gridsize=30,
		stat_func=None, space=0, ratio=5 ,marginal_kws=dict(bins=50, color='r')) 
		 #xlim=(-15,195), ylim=(-15,195))
	g.set_axis_labels(xData + ' / degrees',yData + ' / degrees')
	anText = "Data from " + structData # add timestamp/version
	g.fig.text(0.5, 0.05, anText, horizontalalignment='center', verticalalignment='bottom', fontsize = 10)	
	versionText = "CDS version " + str(int(csd_version())*0.01) + ". Date: " + datetime.datetime.now().strftime ("%m/%d/%Y")
	g.fig.text(0.5, 0, versionText, horizontalalignment='center', verticalalignment='bottom', fontsize = 10)	
	g.fig.suptitle(title)
	#set these different for t vs a 
	if xData[0] == 'T':	
		g.ax_joint.xaxis.set_major_locator(ticker.MultipleLocator(30))
	if yData[0] == 'T':	
		g.ax_joint.yaxis.set_major_locator(ticker.MultipleLocator(30))
	plt.setp(g.ax_marg_y.patches, color='g')
	plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.2) 
	cax = g.fig.add_axes([.9, .4, .01, .25])  # x, y, width, height
	plt.colorbar(cax=cax, format='%.0f')
	plt.show()