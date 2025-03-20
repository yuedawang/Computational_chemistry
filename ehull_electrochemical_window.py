import os
import pandas as pd
from pymatgen.core import Composition, Element
from pymatgen.apps.borg.queen import BorgQueen
from pymatgen.apps.borg.hive import VaspToComputedEntryDrone
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility, MaterialsProjectCompatibility
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter
from pymatgen.ext.matproj import MPRester

# Get the directory of the current script and change the working directory
script_dir = os.path.dirname(os.path.realpath(__file__))
os.chdir(script_dir)

# Load entries from VASP output
queen = BorgQueen(VaspToComputedEntryDrone(), number_of_drones=1)
queen.load_data("vasp_data.gz")
vasp_entries = queen.get_data()

if not vasp_entries:
    raise ValueError("Error: No entries found in vasp_data.gz!")

entry = vasp_entries[0]  # Single entry analysis

# Fix run type and apply energy correction
def correct_energy(entry):
    if hasattr(entry, "run_type") and entry.run_type == "PBE":
        entry.run_type = "GGA"
    if "run_type" in entry.parameters and entry.parameters["run_type"] == "PBE":
        entry.parameters["run_type"] = "GGA"
    entry.parameters.setdefault("is_hubbard", False)
    entry.parameters.setdefault("hubbards", {})
    entry = MaterialsProjectCompatibility().process_entry(entry)
    a = MaterialsProject2020Compatibility().get_adjustments(entry)
    entry.energy_adjustments = a
    return entry

corrected_entry = correct_energy(entry)

# Retrieve MP entries for the same chemical system
API_KEY = "your_api_key_here"
with MPRester(API_KEY) as m:
    chem_sys = list({el.symbol for el in entry.composition.elements} | {'Li'})
    mp_entries = m.get_entries_in_chemsys(chem_sys)

# Construct phase diagram
all_entries = mp_entries + [corrected_entry]
pd = PhaseDiagram(all_entries)
e_above_hull = pd.get_e_above_hull(corrected_entry)

# Write the outputs to stability.txt
with open("stability.txt", "w") as f:
    f.write("Chemical system:\n")
    f.write(str(chem_sys) + "\n")
    f.write("E above hull (eV/atom): %s\n" % e_above_hull)
    
    li_entries = [e for e in all_entries if e.composition.reduced_formula == "Li"]
    uli0 = min(li_entries, key=lambda e: e.energy_per_atom).energy_per_atom
    
    el_profile = pd.get_element_profile(Element("Li"), corrected_entry.composition)
    for i, d in enumerate(el_profile):
        voltage = -(d["chempot"] - uli0)
        f.write("Voltage: %s V\n" % voltage)
        f.write(str(d["reaction"]) + "\n")
        f.write("\n")
