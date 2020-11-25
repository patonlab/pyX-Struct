[![DOI](https://zenodo.org/badge/137907550.svg)](https://zenodo.org/badge/latestdoi/137907550)


=========

**pyX-Struct** provides a python interface to scrape geometric X-ray Data from the [Cambridge Structural Database](https://www.ccdc.cam.ac.uk/solutions/csd-system/components/csd/).<sup>1</sup> This code has been used to carry out a quantitative comparison of the conformational preferences of diarylureas and diarylthioureas in the solid state.<sup>2</sup> 

## Getting Started 

This Python program is run from Terminal or Unix Shell in a Python environment that
contains the Cambridge Crystallographic Data Centre library (CCDC). 

The program reads a SMILES string as a substructure to search the database with, along with additional optional arguments that allow the user to request measurements such as the distance between two atoms, an angle between three atoms, a
torsion angle between four atoms, or any combination or number of these three measurements from crystallographic data containing the input substructure. Measurements (distances, angles, torsion angles) may be compared graphically. 

If the measurement type and indices of the involved atoms are known, they
may be included in the initial command argument, however, the program will still ask if any additional measurements should be added to the structure.

The user has the option to export the data as a .CSV file which will save the resulting molecule's CSD identifier along with any specified substructure measurements. 

The program may optionally search for hydrogen bonding from urea / thiourea nitrogens and 
save the count results as a CSV.
	
#### Package Dependencies
csd-python-api, future, matplotlib, numpy, pandas, seaborn

#### Optional Arguments
* A SMILES string (encased in quotations if illegal characters are involved) of a molecule will search for crystal structures including the substructure.
* The `d` argument followed by two atom indices measures the distance between the two given atoms.
* The `a` argument followed by three atom indices will measure an angle between the three given atoms.
* The `t` argument followed by four atom indices will measure a torsion angle of the four given atoms.
* The `s` argument will save the crystal identifiers and specified measurement search data as a .CSV file in the current directory.
* The `lim` argument allows the user to specify a limit to the number of search results obtained, default limit is 1000 crystal structures.
* The `p` argument will print search result data to the command line as the found crystal structure identifiers and specified measurements.
* The `g` argument turns graphing of two measurements off, default behavior displays graph.
* The `h` argument permits the search of urea or thiourea hydrogen bonding activity.
	

## Sample Inputs/Outputs

#### Example 1: Search the CSD for a porphyrin ring substructure removing the search limit and exporting the results as a .CSV file
```
python structure_search.py 'C1=CC2=NC1=CC3=NC(=CC4=NC(=CC5=NC(=C2)C=C5)C=C4)C=C3' lim 0 s
Data will be saved.
Search for any specific measurements on this molecule? (y/n): n
Searching for substructures...
Found 18 matching substructures in 12 different molecules.
File saved to: ./search_16:38:39.CSV
```
Output file: [search_16:38:39.CSV](https://github.com/bobbypaton/pyX-Struct/blob/master/search_16:38:39.CSV)

#### Example 2: Search the CSD for a 2-chlorobut-2-ene substructure and measure a torsion angle
```
python structure_search.py 'CC=C(C)Cl' t 0 1 2 4 s
Torsion TOR0 added to the search.
Data will be saved.
Searching for substructures with a limit of 1000 max structures...
Found 2238 matching substructures in 1000 different molecules.
File saved to: ./search_16:42:05.CSV
```
Output file: [search_16:42:05.CSV](https://github.com/bobbypaton/pyX-Struct/blob/master/search_16:42:05.CSV)

#### Example 3: Search the CSD for an ethanol substructure measuring C-O distance and C-C-O angle
```
python structure_search.py 'CCO' d a lim 
0 QueryAtom(C)[atom aromaticity: equal to 0]
1 QueryAtom(C)[atom aromaticity: equal to 0]
2 QueryAtom(O)[atom aromaticity: equal to 0]
Enter two indices to measure a distance (# #): 0 1
Distance D0 added to the search.
0 QueryAtom(C)[atom aromaticity: equal to 0]
1 QueryAtom(C)[atom aromaticity: equal to 0]
2 QueryAtom(O)[atom aromaticity: equal to 0]
Enter three indices to measure an angle (# # #): 2 1 0
Angle A0 added to the search.
Enter a number for max number of hits for the search: 400
Searching for substructures with a limit of 400 max structures...
Found 1804 matching substructures in 400 different molecules.
Graphing 'D0' vs 'A0'...
```
<p align="center">
<img src=Example2Pic.png alt="Example 2 Picture" width=400 height=400 >
</p>


  
## Known Bugs
* Unable to access indices of hydrogen atoms given a SMILES string.
* UserWarning occurs when generating the graph since 'normed' kwarg is depricated.
	
## Other Notes
* If a SMILES string or measurement are not specified, or if indices of a measurement are not known or given,
	interaction with the terminal is required.
* The default search limit is 1000 structures. This can be changed with the `lim` argument.
* If a measurement was added by accident and the program is asking for atom indices for a specific measurement, entering `q` will quit adding the measurement.

#### References
1. “The Cambridge Structural Database” Groom, C. R.; Bruno, I. J.; Lightfoot, M. P.; Ward, S. C.; *Acta Cryst. B*, **2016**, *B72*, 171-179
[**DOI:** 10.1107/S2052520616003954](http://dx.doi.org/10.1107/S2052520616003954)
2. "Data-mining the Diaryl(thio)urea Conformational Landscape: Understanding the Contrasting Behavior of Ureas and Thioureas with Quantum Chemistry" Luchini, G.; Ascough, D. H. M.; Alegre-Requena, J. V.; Gouverneur, V.; Paton, R. S. *Tetrahedron* **2019**, *75*, 697-702 [**DOI:** 10.1016/j.tet.2018.12.033](https://doi.org/10.1016/j.tet.2018.12.033)
