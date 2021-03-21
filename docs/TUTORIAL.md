# How to use the CLI of project: Handbag
Install python 3.8.8 [here](https://www.python.org/downloads/release/python-388/). Scroll down and click the 
Windows 64-bit installer.

Install Windows terminal from the Microsoft Marketplace or click [here](https://www.microsoft.com/en-us/p/windows-terminal/9n0dx20hk701?activetab=pivot:overviewtab).

Download the source code of this page by using the green download button on the github page.
[download button](../figures/tutorial/download_source_code.png)

## First time: Set up a virtual environment and install the dependencies of this project
Navigate to the src folder in file explorer and open a windows terminal by right clicking the white space in the folder and choosing open Windows Terminal.
[open windows terminal](../figures/tutorial/open_terminal.png)
Create your virtual environment by typing ```python -m venv handbag-env``` and pressing enter. It might take some time.
Activate the virtual environment with: ```./handbag-env/Scripts/Activate.ps1```. While typing folder names, you can always use 'tab' to let the terminal auto-complete the filepath you're typing. When the virtual environment is activated, this should be visible. 
[activated venv](../figures/tutorial/activated_venv.png)
Install dependencies of the project: ```pip install -r requirements.txt```. This can take a bit of time.
You can deactivate the virtual environment with saying ```deactivate```, but you don't have to: you can always just close the terminal.

## After the first time:
Open the terminal in the src folder and activate the virtual environment with ```./handbag-env/Scripts/Activate.ps1```

## Testing the installation
Test if everything works with running: ```python .\quantify.py -f .\arg_files\test_h2o_xh_o.txt```. If everything worked correctly, the program should now start to run and you should be able to explore all the plots.
[working program](../figures/tutorial/working_program.png)

## Tips for using the terminal
```cd```: changes directory. Change to directory scripts by typing ```cd scripts```, and back by ```cd ../```, where ```../``` means óne folder "up", so: move to the parent directory of the folder you're currently in.
```ls```: list. List all files and directories that are in the current directory. This can help you navigate through the file structure on your computer.

## Running other data
>! This is a spoiler
Test

