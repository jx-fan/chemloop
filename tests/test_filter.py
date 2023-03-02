from pathlib import Path

import pytest
from chemloop.analysis.filter import ReactionFilter
from monty.serialization import loadfn
from pymatgen.core import Composition
from rxn_network.reactions.basic import BasicReaction

TEST_FILES_PATH = Path(__file__).parent / "test_files"
PATHWAY_FILE = "MnO_Mn2N_paths_1.json.gz"


@pytest.fixture(scope="module")
def pathway_set():
    return loadfn(TEST_FILES_PATH / PATHWAY_FILE)


@pytest.fixture(scope="module")
def reaction_filter(pathway_set):
    rxn_filter = ReactionFilter(rxns_to_include={BasicReaction.balance(reactants=[Composition("H2O"),
                                                                                  Composition("Mn2N")],
                                                                       products=[Composition("H2"),
                                                                                 Composition("MnN"),
                                                                                 Composition("MnO")]
                                                                       )})
    return rxn_filter.filter(pathway_set)


def test_filter(reaction_filter):
    assert len(reaction_filter) == 16
