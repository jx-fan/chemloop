from pathlib import Path

import numpy as np
import pytest
from chemloop.analysis.clas import AnalyseHydroPathwaySet
from monty.serialization import loadfn
from pymatgen.core import Composition, Element
from rxn_network.reactions.basic import BasicReaction

from chemloop.analysis.filter import ReactionFilter
from src.chemloop.analysis.clas import limiting_step, ammonia_yield_steps

TEST_FILES_PATH = Path(__file__).parent / "test_files"
PATHWAY_FILE = "MnO_Mn2N_paths_1.json.gz"
NET_REACTION_DICT = {"oxide": Composition("MnO"),
                     "nitride": Composition("Mn2N"),
                     "net_rxn": BasicReaction.from_string("0.5 Mn2N + H2O -> MnO + 0.5 H3N + 0.25 H2"),
                     "net_rxn_energy": -0.1420334859928212}


@pytest.fixture(scope="module")
def pathway_set():
    return loadfn(TEST_FILES_PATH / PATHWAY_FILE)


@pytest.fixture(scope="module")
def analyser(pathway_set):
    return AnalyseHydroPathwaySet(pathway_set=pathway_set, **NET_REACTION_DICT)


@pytest.fixture(scope="module")
def reaction_to_include():
    return BasicReaction.balance(reactants=[Composition("H2O"), Composition("Mn2N")],
                                 products=[Composition("H2"), Composition("MnN"), Composition("MnO")])


@pytest.fixture(scope="module")
def analyser_with_filter(pathway_set, reaction_to_include):
    return AnalyseHydroPathwaySet(pathway_set=pathway_set,
                                  pathway_filter=ReactionFilter(rxns_to_include={reaction_to_include}),
                                  **NET_REACTION_DICT
                                  )


@pytest.fixture(scope="module")
def lowest_cost_pathway_default():
    return loadfn(TEST_FILES_PATH / "lowest_arithmetic_cost_pathway_default.json.gz")


@pytest.fixture(scope="module")
def lowest_cost_pathway_with_filter():
    return loadfn(TEST_FILES_PATH / "lowest_arithmetic_cost_pathway_with_filter.json.gz")


def test_analyser_path_filter(pathway_set, analyser_with_filter, reaction_to_include):
    analyser = AnalyseHydroPathwaySet(pathway_set=pathway_set, **NET_REACTION_DICT)
    assert len(analyser.paths) != len(analyser_with_filter.paths)
    analyser.pathway_filter = ReactionFilter(rxns_to_include={reaction_to_include})
    assert len(analyser.paths) == len(analyser_with_filter.paths)
    assert analyser.lowest_cost_pathway == analyser_with_filter.lowest_cost_pathway


def test_pathway_set(analyser, pathway_set):
    assert analyser.pathway_set == pathway_set


def test_paths(analyser, pathway_set):
    assert analyser.paths == sorted(pathway_set.get_paths(), key=lambda p: float(np.mean(p.costs)))


def test_lowest_cost_pathway(analyser, analyser_with_filter,
                             lowest_cost_pathway_default,
                             lowest_cost_pathway_with_filter):
    assert analyser.lowest_cost_pathway == lowest_cost_pathway_default
    assert analyser_with_filter.lowest_cost_pathway == lowest_cost_pathway_with_filter


def test_lowest_cost(analyser, analyser_with_filter):
    assert analyser.lowest_cost == pytest.approx(0.2727760806494658)
    assert analyser_with_filter.lowest_cost == pytest.approx(0.2784896832452982)


def test_cations(analyser):
    assert analyser.cations == [Element("Mn")]


def test_from_file():
    analyser_from_file = AnalyseHydroPathwaySet.from_file(file_pathway=str(TEST_FILES_PATH / PATHWAY_FILE),
                                                          file_energy=str(TEST_FILES_PATH / "results_tot_773K.csv"),
                                                          rxn_column_name="sub_rxn1",
                                                          e_column_name="e1"
                                                          )
    assert analyser_from_file.oxide == Composition("MnO")


def test_net_rxn(analyser):
    assert analyser.net_rxn == BasicReaction.from_string("0.5 Mn2N + H2O -> MnO + 0.5 H3N + 0.25 H2")


def test_net_rxn_energy(analyser):
    assert analyser.net_rxn_energy == -0.1420334859928212


def test_nitride(analyser):
    assert analyser.nitride == Composition("Mn2N")


def test_oxide(analyser):
    assert analyser.oxide == Composition("MnO")


def test_ammonia_yield_steps(analyser):
    reactions, energy = ammonia_yield_steps(pathway=analyser.lowest_cost_pathway,
                                            net_rxn=analyser.net_rxn)
    assert reactions == analyser.lowest_cost_pathway.reactions[-2:-1]
    assert energy == pytest.approx(0.6686309157599197)


def test_softplus(analyser):
    assert analyser.softplus(900, 0) == pytest.approx(0.26492508532916464)


def test_net_rxn_cost(analyser):
    assert analyser.net_rxn_cost == pytest.approx(0.26728068094549096)


def test_temperature(analyser):
    assert analyser.temperature == 773


def test_limiting_step(analyser):
    reaction, energy = limiting_step(analyser.lowest_cost_pathway,
                                     analyser.net_rxn)
    assert reaction == analyser.lowest_cost_pathway.reactions[0]
    assert energy == pytest.approx(11.301137676060156)
