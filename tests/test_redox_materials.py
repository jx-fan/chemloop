import pytest
from pymatgen.core import Composition, Element

from src.chemloop.core.redox_materials import RedoxMaterialsSet


@pytest.fixture
def manganese_set() -> RedoxMaterialsSet:
    return RedoxMaterialsSet(materials=[Composition("Mn2N"), Composition("MnO")])


def test_add(manganese_set):
    original_len = len(manganese_set)
    manganese_set.add(Composition("Mn"))
    assert len(manganese_set) == original_len + 1
    assert Composition("Mn") in manganese_set


def test_discard(manganese_set):
    original_len = len(manganese_set)
    manganese_set.discard(Composition("Mn2N"))
    assert len(manganese_set) == original_len - 1
    assert Composition("Mn") not in manganese_set


def test_num_materials(manganese_set):
    assert manganese_set.num_materials == 2


def test_cations(manganese_set):
    assert manganese_set.cations == {Element("Mn")}


def test_anions(manganese_set):
    assert manganese_set.anions == {Element("N"), Element("O")}


def test_chemical_system(manganese_set):
    assert manganese_set.chemical_system == {"Mn", "N", "O"}


def test_add_charges():
    assert False
