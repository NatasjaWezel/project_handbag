# Installation
Install python 3.8.8 [here](https://www.python.org/downloads/release/python-388/). Scroll down and click the 
Windows 64-bit installer.
Install Windows terminal from the Microsoft Marketplace or click [here](https://www.microsoft.com/en-us/p/windows-terminal/9n0dx20hk701?activetab=pivot:overviewtab).

Follow the detailed tutorial [here](docs/TUTORIAL.md) if you're not used to working with terminal.

Download the source code from GitHub, open a terminal in the src folder

# Getting data from the CSD
Search query in conquest, export fragments xyz coordinates. ConQuest search. Definitions of the groups used.

# Running
Start your virtual environment ```python venv handbag```
Activate the environment
Install dependencies of the project: ```pip install -r docs/requirements.txt```

Specify the labels from your central groups in the central groups [file](../src/files/central_groups.csv). You can open this file and adjust it with any code editor you like, notepad, excel, etc.

If your files have the nameformat CENTRAL_CONTACT.cor and CENTRAL_CONTACT.csv, you can now run the program with: ```python quantify.py --input INPUTFILE --crp CONTACT_RP```
Else, you can specify the names with the extra options.
* --labels : the csv file containing the parameters, from conquest
* --central : the name of your central group
* --contact : the name of your contact grou

## Specifying argument files 
If you do not want to type all the input arguments each time, you can also create a file containing the specifications. Examples for the test data are given [here](./src/arg_files/).


# Using single scripts

## Step 2: Aligning fragments
Align the fragments by running `python load_from_coords.py <inputfile.xyz> <outputfile.csv>`. 

Choosing atoms to align on

Plot the overlapping fragments by running `python plot_result.py <inputfile.csv>`.

## Step 3: Plot density
Plot the 4D graph by running `python plot_density.py <inputfile.csv> <resolution>`.
Plot the density

## Step 4: Analyzing the density - are there any clusters?
Get the directionality coefficient
