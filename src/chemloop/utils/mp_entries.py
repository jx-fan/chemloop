import json
import os
from pathlib import Path
from typing import List, Set
from dotenv import load_dotenv

from chemloop.data import CL_AMMONIA_PATH, CL_COMBUSTION_PATH
from mp_api.client import MPRester
from pymatgen.entries.computed_entries import ComputedStructureEntry


def get_entries_from_api(chemsys) -> List[ComputedStructureEntry]:
    load_dotenv()
    mp_api_key = os.getenv("MP_API_KEY") or ""
    with MPRester(mp_api_key) as mpr:
        entries = mpr.get_entries_in_chemsys(chemsys)
        print("Materials Project database version: ", mpr.get_database_version())
    return entries


def get_entries_from_json(cations_in_redox_materials: List[str],
                          chemsys_net_rxn: Set[str]) -> List[ComputedStructureEntry]:
    entries = []
    entry_paths = _get_file_paths(cations_in_redox_materials, chemsys_net_rxn)
    for path in entry_paths:
        with open(path, "r") as f:
            entries.extend(json.load(f))
    return [ComputedStructureEntry.from_dict(e) for e in entries]


def _get_file_paths(cations_in_redox_materials: List[str],
                    chemsys_net_rxn: Set[str]
                    ) -> List[Path]:
    entry_paths = []
    if chemsys_net_rxn == {"O", "N", "H"}:  # chemical looping ammonia. Included the intermediate H2O.
        base_entries_path = CL_AMMONIA_PATH / "computed_str_entries_O_N_H_2021_11_10.json"
        entry_paths.append(base_entries_path)
        if len(cations_in_redox_materials) == 1:
            binary_fn = "computed_str_entries_" + cations_in_redox_materials[0] + "_O_N_H_2021_11_10.json"
            entry_paths.append(CL_AMMONIA_PATH / binary_fn)
            return entry_paths
        elif len(cations_in_redox_materials) == 2:
            for c in cations_in_redox_materials:
                binary_fn = "computed_str_entries_" + c + "_O_N_H_2021_11_10.json"
                entry_paths.append(CL_AMMONIA_PATH / binary_fn)
            ternary_fn = "computed_str_entries_" + "_".join(cations_in_redox_materials) + "_O_N_H_2021_11_10.json"
            entry_paths.append(CL_AMMONIA_PATH / ternary_fn)
            return entry_paths
    elif chemsys_net_rxn == {"O", "C", "H"}:
        base_entries_path = CL_COMBUSTION_PATH / "computed_str_entries_O_H_C_2021_11_10.json"
        entry_paths.append(base_entries_path)
        binary_fn = "computed_str_entries_" + cations_in_redox_materials[0] + "_O_H_C_2021_11_10.json"
        entry_paths.append(CL_COMBUSTION_PATH / binary_fn)
        return entry_paths
