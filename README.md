# Finding and quantifying the directionality of intermolecular interactions
This is a project that I do for 't van het Hoff Institute for Molecular Sciences at the UvA.

Supervisor: dhr. dr. T.J. (Tiddo) Mooibroek

Auteur: Natasja Wezel (BSc)


### Installation - Scripting
```pip install -r docs/requirements.txt```

### Preparation
Search query in conquest, export fragments xyz coordinates.

### What it does
It takes as input orthogonal data from fragments. You can specify the central and contact groups.

### Running

#### Step 1: Getting data from the CSD
ConQuest search. Definitions of the groups used.

#### Step 2: Aligning fragments
Align the fragments by running `python load_from_coords.py <inputfile.xyz> <outputfile.csv>`. 

Choosing atoms to align on

Plot the overlapping fragments by running `python plot_result.py <inputfile.csv>`.

#### Step 3: Plot density
Plot the 4D graph by running `python plot_density.py <inputfile.csv> <resolution>`.
Plot the density

#### Step 4: Analyzing the density - are there any clusters?
Get the directionality coefficient

[![BCH compliance](https://bettercodehub.com/edge/badge/NatasjaWezel/MasterProject?branch=master)](https://bettercodehub.com/)


## Contributing
Pull request templates are provided. Issues and pull requests are always welcome.
For other questions, you can always contact me!
