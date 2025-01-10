Welcome to Navground Blofun, standing for Navigation playground block functions. 
This repository is for using topological data analysis techniques for the study of 
autonomous robots navigation simulations performed using Navground. 
Our analysis is based on the induced partial matching, a technique to match clusters of 
two datasets based on underlying assignment, such as when a dataset is a subset of another. 

Installation Instructions
=========================
Dependencies
------------
- Python 3
- Numpy
- Matplotlib
- Navground
- Scipy
- tslearn
- seaborn

In order to install Navground, go to a directory where you want to save the repository. 
Next, open the command line and execute the following instructions to clone the repository.
```sh
    git clone git@github.com:Cimagroup/navground-blofun.git
    cd navground-blofun
```
Optionally, you might wish to create and activate a virtual environment as follows.
```sh
    python3 -m venv venv
    . venv/bin/activate
```
Finally, install the remaining dependencies and install the perdiver module, which is contained in this repository.
```sh
    python3 -m pip install numpy matplotlib scipy navground jupyter notebook tslearn seaborn
    python3 -m pip install .
```
Now you are ready to execute the notebooks. Run the following commands and you are ready to go!
```sh
    cd notebooks
    python3 -m jupyter notebook
```
