from __future__ import print_function
try:
   from builtins import input
except ImportError:
   from __builtin__ import input
import datetime
import os
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from ccdc.search import SMARTSSubstructure, SubstructureSearch
from ccdc import io
from ccdc.molecule import Molecule
from hbond_search import hbond_df_search

#initial variables
search = True             #True: the program will perform a search on the input SMILES
graph = True              #True: the program will display angle or torsion data
searchLimit = True        #True: the search is limited by the maxHits variable
printData = False         #True: the data will print in the terminal for the user
saveData = False          #True: Export data
hbond = False             #False: Do not search for urea/thiourea hydrogen bonds by default
maxHits = 1000            #Max number of hits until the search ends, default at 1000

#number of measurements for naming purposes d, a, t respectively
measurement = np.zeros(shape=(3,2))

#for adding a distance measurement
def addDist(atom = None, askIndex = False,):
    global substructure, substructure_search, sub_id, measurement
    name = 'D'+str(int(measurement[0][1]))
    addM = True
    if not askIndex:
        substructure_search.add_distance_measurement(name,
            sub_id, atom[0], sub_id, atom[1])
        print('Distance', name, 'added to the search.')
        measurement[0][1] += 1
    else:
        for i in substructure.atoms:
            print(i.index, i)
        while True:
            try:
                atom_i = input("Enter two indices to measure a distance (# #): ")
                atom_i = list(map(int, atom_i.split()))
            except ValueError:
                if atom_i == 'q':
                    addM = False
                    measurement[0][0]-=1
                    break
                else:
                    print('Try Again with two numbers')
                    continue
            if len(atom_i) != 2:
                print('Try Again with two numbers')
                continue
            elif not all(i < len(substructure.atoms) for i in atom_i):
                print('Enter indices within the range of given atoms')
                continue
            elif not sorted(atom_i) == list(sorted(set(atom_i))):
                print('Enter two different numbers')
                continue
            else:
                break
        if addM:
            substructure_search.add_distance_measurement(name,
                sub_id, atom_i[0], sub_id, atom_i[1])
            print('Distance', name, 'added to the search.')
            measurement[0][1] += 1


#for adding an angle measurement
def addAngle(atom = None, askIndex = False):
    global substructure, substructure_search, sub_id, measurement
    name = 'A'+str(int(measurement[1][1]))
    addM = True
    if not askIndex:
        substructure_search.add_angle_measurement(name,
            sub_id, atom[0], sub_id, atom[1], sub_id, atom[2])
        print('Angle', name, "added to the search.")
        measurement[1][1] += 1
    else:
        for i in substructure.atoms:
            print(i.index, i)
        while True:
            try:
                atom_i = input("Enter three indices to measure an angle (# # #): ")
                atom_i = list(map(int, atom_i.split()))
            except ValueError:
                if atom_i == 'q':
                    addM = False
                    measurement[1][0]-=1
                    break
                else:
                    print('Try Again with three numbers')
                    continue
            if len(atom_i) != 3:
                print('Try Again with three numbers')
                continue
            elif not all(i < len(substructure.atoms) for i in atom_i):
                print('Enter indices within the range of given atoms')
                continue
            elif not sorted(atom_i) == list(sorted(set(atom_i))):
                print('Enter three different numbers')
                continue
            else:
                break
        if addM:
            substructure_search.add_angle_measurement(name,
                sub_id, atom_i[0], sub_id, atom_i[1], sub_id, atom_i[2])
            print('Angle', name, "added to the search.")
            measurement[1][1] += 1

#for adding a torsion angle measurement
def addTor(atom = None, askIndex = False):
    global substructure, substructure_search, sub_id, measurement
    name = 'TOR'+str(int(measurement[2][1]))
    addM = True
    if not askIndex:
        substructure_search.add_torsion_angle_measurement(name,
            sub_id, atom[0], sub_id, atom[1], sub_id, atom[2], sub_id, atom[3])
        print('Torsion', name, "added to the search.")
        measurement[2][1] += 1
    else:
        for i in substructure.atoms:
            print(i.index, i)
        while True:
            try:
                atom_i = input("Enter four indices to measure a torsion angle (# # # #): ")
                atom_i = list(map(int, atom_i.split()))
            except ValueError:
                if atom_i == 'q':
                    addM = False
                    measurement[2][0]-=1
                    break
                else:
                    print('Try Again with four numbers')
                    continue
            if len(atom_i) != 4:
                print('Try Again with four numbers')
                continue
            elif not all(i < len(substructure.atoms) for i in atom_i):
                print('Enter indices within the range of given atoms')
                continue
            elif not sorted(atom_i) == list(sorted(set(atom_i))):
                print('Enter four different numbers')
                continue
            else:
                break
        if addM:
            substructure_search.add_torsion_angle_measurement(name,
                sub_id, atom_i[0], sub_id, atom_i[1], sub_id, atom_i[2],sub_id, atom_i[3])
            print('Torsion', name, "added to the search.")
            measurement[2][1] += 1

#parse out user given arguments
def parseArgs(arg):
    global graph, printData, maxHits, saveData
    #check for distance measurement
    if 'd' in arg:
        dIndices = [j for j, k in enumerate(arg) if k == 'd']
        measurement[0][0] = len(dIndices)
        for j in np.arange(len(dIndices)):
            atom = []
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
        measurement[1][0] = len(aIndices)
        for j in np.arange(len(aIndices)):
            atom = []
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
        measurement[2][0] = len(tIndices)
        for j in np.arange(len(tIndices)):
            atom = []
            try:
                for k in np.arange(1,5):
                    atom.append(int(arg[tIndices[j]+k]))
                if all(isinstance(x,int) for x in atom):
                    addTor(atom=atom)
                else:
                    addTor(askIndex=True)
            except (ValueError, IndexError):
                addTor(askIndex=True)
    #check if search limit was specified
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
            lim = int(input("Enter a number for max number of hits for the search: "))
            maxHits = lim
    #check if user wants to save data
    if 's' in arg:
        saveData = True
        print('Data will be saved.')
    #check if user wants to print data to terminal
    if 'p' in arg:
        printData = True
        print('Data will be printed.')
    #check if user does not want to display a graph
    if 'g' in arg:
        graph = False
        print('Data will not be graphed.')

    if 'h' in arg:
        hbond = True
        print('Hbond data will be searched for.')

def main():
    #if the user didn't pass any arguments, ask for a SMILES string
    if len(sys.argv) == 1:
        smiles = input("Enter a SMILES string: ")
        if smiles[0] == ("'" or '"'):
            smiles = smiles[1:]
        if smiles[-1] == ("'" or '"'):
            smiles = smiles[:-1]
    else:
        smiles = sys.argv[1]

    #Prepare a search, checking for invalid submissions
    while True:
        try:
            substructure = SMARTSSubstructure(smiles)
            if smiles[0] == ("'" or '"'):
                smiles = smiles[1:]
            if smiles[-1] == ("'" or '"'):
                smiles = smiles[:-1]
        except Exception as e:
            print(e)
            smiles = input('Enter another SMILES string:')
            continue
        else:
            break

    substructure_search = SubstructureSearch()
    sub_id = substructure_search.add_substructure(substructure)

    #determine if the user has specified other measurements to keep track of
    if len(sys.argv) > 2:
        arguments = sys.argv
        del arguments[0:2]
        parseArgs(arguments)
        if not any(i in sys.argv for i in ['d', 'a', 't']):
            ans = input('Search for any specific measurements on this molecule? (y/n): ')
            while ans == 'y':
                ans = input('Which measurements? \nEnter d for distance, a for angle, t for torsion angle: ')
                ans = list(map(str, ans.split()))
                parseArgs(ans)
                ans = input('Add any additional measurements? (y/n): ')
    else:
        ans = input('Search for any specific measurements on this molecule? (y/n): ')
        while ans == 'y':
            #ask which one, run method
            ans = input('Which measurements? \nEnter d for distance, a for angle, t for torsion angle: ')
            ans = list(map(str, ans.split()))
            parseArgs(ans)
            ans = input('Add any additional measurements? (y/n): ')


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

    #format data into lists and place it into a dataframe
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

    ids = []
    measureHits = [[] for _ in range(len(columns)-1)]
    for h in hits:
        ids.append(h.identifier)
        if len(columns) > 1:
            for i in np.arange(len(columns)-1):
                measureHits[i].append(float('%.2f' % h.measurements[columns[i+1]]))
    hitData = pd.DataFrame(columns=columns)
    hitData['Identifier'] = ids
    for i in np.arange(len(measureHits)):
        hitData[columns[i+1]] = np.absolute(measureHits[i])

    #print useful search data info to user
    structData = str(len(hitData.Identifier)) + ' matching substructures in '+ str(len(hitData.Identifier.unique())) + ' different molecules.'
    print('Found ' + structData)

    #optionally export data to user
    if saveData:
        filename = os.getcwd() + '/search_' + datetime.datetime.now().strftime ("%H:%M:%S") + '.CSV'
        hitData.to_csv(filename)
        print('File saved to: ' + filename)

    if len(hitData.columns) < 3:
        graph = False

    #optionally print data to user
    if printData:
        print(hitData.to_string())

    #optionally graph, asking user which data to display, formatting graph and axes titles
    if graph:
        if len(hitData.columns) == 3:
            xData = hitData.columns[1]
            yData = hitData.columns[2]
        elif len(hitData.columns) > 3:
            print('Which data would you like to display? (enter two indices, # #):')
            print('i | Measurement')
            for i in np.arange(len(columns[1:])):
                print(str(i) + " | " + str(columns[i+1]))
            while True:
                try:
                    ans = input()
                    ans = list(map(int, ans.split()))
                except ValueError:
                    print('Try Again with two numbers')
                    continue
                if len(ans) != 2:
                    print('Try Again with two numbers')
                    continue
                elif not all(i < len(columns)-1 for i in ans):
                    print('Enter indices within the range of given measurements')
                    continue
                elif not sorted(ans) == list(sorted(set(ans))):
                    print('Enter two different numbers')
                    continue
                else:
                    break
            xData = hitData.columns[ans[0]+1]
            yData = hitData.columns[ans[1]+1]
        if xData[0] == 'T' and yData[0] == 'T':
            title = 'Dihedral Angle Magnitude Comparison'
            xax = xData + ' / degrees'
            yax = yData + ' / degrees'
        elif xData[0] == 'A' and yData[0] == 'A':
            title = 'Angle Magnitude Comparison'
            xax = xData + ' / degrees'
            yax = yData + ' / degrees'
        elif xData[0] == 'D' and yData[0] == 'D':
            title = 'Distance Comparison'
            xax = xData + r' / $\AA$'
            yax = yData + r' / $\AA$'
        elif any(i in xData for i in ['A', 'T']) and any(i in yData for i in ['A', 'T']):
            title = 'Dihedral Angle vs. Angle Magnitude Comparison'
            xax = xData + ' / degrees'
            yax = yData + ' / degrees'
        elif any(i in xData for i in ['D']) and any(i in yData for i in ['A', 'T']):
            title = 'Distance vs. Angle Magnitude Comparison'
            xax = xData + r' / $\AA$'
            yax = yData + ' / degrees'
        elif any(i in xData for i in ['A', 'T']) and any(i in yData for i in ['D']):
            title = 'Angle Magnitude Comparison vs Distance'
            xax = xData + ' / degrees'
            yax = yData + r' / $\AA$'
        else:
            title = 'Angular Magnitude Comparison'


        print('Graphing \'' + xData + '\' vs \'' + yData + '\'...')

        #Joint plot Creator
        #gridsize: larger # = smaller hexagons
        #marginal_kws bins: smaller # = thicker histograms
        sns.set(style='white', color_codes=True)
        g = sns.jointplot(x=xData, y=yData,data=hitData, kind="hex",color='b', gridsize=30,
            stat_func=None, space=0, ratio=5 ,marginal_kws=dict(bins=50, color='r'))
        g.set_axis_labels(xax,yax)
        anText = "Data from " + structData # add timestamp/version
        g.fig.text(0.5, 0.05, anText, horizontalalignment='center', verticalalignment='bottom', fontsize = 10)
        versionText = "CDS version " + str(int(io.csd_version())*0.01) + ". Date: " + datetime.datetime.now().strftime ("%m/%d/%Y")
        g.fig.text(0.5, 0, versionText, horizontalalignment='center', verticalalignment='bottom', fontsize = 10)
        g.fig.suptitle(title)
        #add specific ticks for torsion angles
        if xData[0] == 'T':
            g.ax_joint.xaxis.set_major_locator(ticker.MultipleLocator(30))
        if yData[0] == 'T':
            g.ax_joint.yaxis.set_major_locator(ticker.MultipleLocator(30))
        plt.setp(g.ax_marg_y.patches, color='g')
        plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.2)
        cax = g.fig.add_axes([.9, .4, .01, .25])
        plt.colorbar(cax=cax, format='%.0f')
        plt.show()

    if hbond:
        #obtain urea nitrogens for searching
        n1_label = []
        n2_label = []
        for i in range(len(hits)):
            new_molecule = True
            for j in range(len(hits[i].match_atoms())):
                if hits[i].match_atoms()[j].atomic_number == 7:
                    if new_molecule:
                        n1_label.append(hits[i].match_atoms()[j].label)
                        new_molecule = False
                    else:
                        n2_label.append(hits[i].match_atoms()[j].label)

        hitData['N1'] = n1_label
        hitData['N2'] = n2_label

        hitData = hbond_df_search(hitData)

        filename = os.getcwd() + '/search_' + datetime.datetime.now().strftime ("%H:%M:%S") + '.CSV'
        hitData.to_csv(filename)
        print('File saved to: ' + filename)

if __name__ == "__main__":
    main()
