from pathlib import Path

import numpy as np
import pytest
from chemloop.analysis.clas import AnalyseHydroPathwaySet, ammonia_yield_energy
from monty.serialization import loadfn
from pymatgen.core import Composition, Element
from rxn_network.reactions.basic import BasicReaction

TEST_FILES_PATH = Path(__file__).parent / "test_files"
PATHWAY_FILE = "NaMnO2_Mn2N_paths.json.gz"
NET_REACTION_DICT = {"oxide": "NaMnO2",
                     "nitride": "Mn2N",
                     "net_rxn": "H2O + 0.5 Mn2N + NaHO -> NaMnO2 + 0.5 H3N + 0.75 H2",
                     "net_rxn_energy": -0.0792372224581654}


@pytest.fixture(scope="module")
def pathway_set():
    return loadfn(TEST_FILES_PATH / PATHWAY_FILE)


@pytest.fixture(scope="module")
def analyser(pathway_set):
    return AnalyseHydroPathwaySet(pathway_set=pathway_set,
                                  net_rxn=BasicReaction.from_string(NET_REACTION_DICT["net_rxn"]),
                                  net_rxn_energy=NET_REACTION_DICT["net_rxn_energy"],
                                  nitride=Composition(NET_REACTION_DICT["nitride"]),
                                  oxide=Composition(NET_REACTION_DICT["oxide"])
                                  )


def test_paths(analyser, pathway_set):
    assert analyser.paths == sorted(pathway_set.get_paths(), key=lambda p: float(np.mean(p.costs)))


def test_lowest_cost_pathway(analyser, pathway_set):
    assert analyser.lowest_cost_pathway == sorted(pathway_set.get_paths(), key=lambda p: float(np.mean(p.costs)))[0]


def test_lowest_cost(analyser):
    assert analyser.lowest_cost == pytest.approx(0.28960340748)


def test_cations(analyser):
    assert analyser.cations == [Element("Mn")]


def test_from_file():
    analyser_from_file = AnalyseHydroPathwaySet.from_file(file_pathway=TEST_FILES_PATH / PATHWAY_FILE,
                                                          file_energy=TEST_FILES_PATH / "results_tot_NaOH_773K.csv",
                                                          rxn_column_name="sub_rxn",
                                                          e_column_name="e",
                                                          normalize_to="NH3")
    assert analyser_from_file.oxide == Composition("NaMnO2")


def test_net_rxn(analyser):
    assert analyser.net_rxn == BasicReaction.from_string("H2O + 0.5 Mn2N + NaHO -> NaMnO2 + 0.5 H3N + 0.75 H2")


def test_net_rxn_energy(analyser):
    assert analyser.net_rxn_energy == -0.0792372224581654


def test_nitride(analyser):
    assert analyser.nitride == Composition("Mn2N")


def test_oxide(analyser):
    assert analyser.oxide == Composition("NaMnO2")


def test_ammonia_yield_energy(analyser):
    energy = ammonia_yield_energy(pathway=analyser.lowest_cost_pathway)
    assert energy == pytest.approx(-0.0889902456607885)
