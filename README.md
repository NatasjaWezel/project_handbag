# Finding and quantifying the directionality of intermolecular interactions
This is a project that I do for 't van het Hoff Institute for Molecular Sciences at the UvA.

Supervisor: dhr. dr. T.J. (Tiddo) Mooibroek

Auteur: Natasja Wezel (BSc)


### Installation
pip install requirements.txt

### Preparation
Search query in conquest, export fragments xyz coordinates.

### What it does
It takes as input orthogonal data from fragments. You can specify the central and contact groups.

### Running
Align the fragments by running `python load_from_coords.py <inputfile.xyz> <outputfile.csv>`. 

Plot the overlapping fragments by running `python plot_result.py <inputfile.csv>`.

Plot the 4D graph by running `python plot_density.py <inputfile.csv> <resolution>`.


[![BCH compliance](https://bettercodehub.com/edge/badge/NatasjaWezel/MasterProject?branch=master)](https://bettercodehub.com/)

