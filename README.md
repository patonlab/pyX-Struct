# pyX-Struct

Scrape Geometric X-ray Data from the Cambridge Structural Database 

## Getting Started 

The program is run from Terminal or Unix Shell in a python environment that
contains the Cambridge Crystallographic Data Centre library (CCDC). The program
may take in in a SMILES string along with additional optional arguments that 
can measure the distance between two atoms, an angle between three atoms, or a 
torsion angle between four atoms by searching the Cambridge Structural
Database of crystal structures and may be used to graphically compare results.
If the measurement type and indices of the involved atoms are known, they
may be included in the command argument. If the indices are not known, the
types of measurement may still be specified. If no measurements are 
specified, the user will later be given the option to add measurements.
A structure search may still take place with no parameters. User has the 
option to export the data as a .CSV file.
	
## Package Dependencies
| Name           |
| :------------: |
| csd-python-api | 
| future         |
| matplotlib     | 
| numpy          | 
| pandas         | 
| seaborn        | 


## Optional Arguments
* SMILES string (encased in quotations if illegal characters are involved)
* d (distance measurement) followed by two atom indices (also optional)
* a (angle measurement) followed by three atom indices (also optional)
* t (torsion angle measurement) followed by three atom indices (also optional)
* s (save data) exports data as a c2m file 
* lim (search hit limit) followed by a number to limit the number of search hits
* p (print) prints the data to the command line
* g (graph) turns graphing of two measurements off (default behavior displays graph)
	
### Sample Arguments
Measures two torsion angles with specified indices:
```
'S=C([NH]c1ccccc1)[NH]c1ccccc1' t 0 1 2 3 t 0 1 9 10
```
Measures two torsion angles, user needs to later specify indices:
```
'O=C([NH]c1ccccc1)[NH]c1ccccc1' t t lim 5000
```
Measures distance between CC and angle between CCO:
```
CCO d 0 1 a 0 1 2 
```
Removes a search limit and exports the resulting structural data:
```
'C1=CC2=NC1=CC3=NC(=CC4=NC(=CC5=NC(=C2)C=C5)C=C4)C=C3' lim 0 s 
```
  
## Known Bugs
* Unable to access indices of hydrogen atoms given a SMILES string.
* UserWarning occurs when generating the graph since 'normed' kwarg is depricated.
	
## Other Notes
* If a measurement is not specified or indices of a measurement are not given,
	interaction with the terminal is required.
* The default search limit is 1000 structures.

## Authors
* **Guilian Luchini** - *Initial work*
	
