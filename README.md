# ehull_electrochemical_window.py

This script is used to calculate the energy above hull and eletrochemical window (reduction/oxidation potential) for common inorganic solid-state electrolytes.

### environmental requirements
1. pymatgen
2. mp-api
   
### building the environment
1. conda create -n pymatgen
2. conda install pymatgen -c conda-forge -y
3. conda activate pymatgen
4. pip install mp-api

### Usage
1. replace you mp_api in this code "your_api_key_here"
2. go to vasp calculation folder and using the command "pmg analyze ." to get the vasp_data.gz
3. keep the vasp_data.gz and ehull_electrochemical_window.py in the same folder
4. python ehull_electrochemical_window.py

