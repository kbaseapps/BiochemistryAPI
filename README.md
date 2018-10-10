
# BiochemistryAPI
---

This module serves Biochemistry data for KBase

## Release Notes
#### 0.3.0 (10/04/18)
- Updated the search methods to approximate a poor man's search engine. Terms are tokenized on _-; and space, and matching results returned in the same order as they appear in the source. 
- Added a optional limit on search results
- Added a search_reactions method

#### 0.2.0 (06/26/18)
- convert to Python3 module

#### 0.1.4 (06/25/18)
- Add basic compound search methods
- Add 3D compound generation method

#### 0.1.3 (05/04/18)
- Add the ability to perform similarity and substructure searches

#### 0.1.2 (10/20/17)
Adds the ability to generate SVG depictions of chemical structures submitted in
SMILES and InChI format.