import pytest
from pymatgen.core import Composition

from src.chemloop.core.net_reactions import NetReaction


@pytest.fixture
def ammonia_synthesis() -> NetReaction:
    return NetReaction(oxidant=Composition("N2"),
                       reducing_agent=Composition("H2"),
                       products=[Composition("NH3")],
                       coefficients=[-1, -3, 2]
                       )


def test_compositions(ammonia_synthesis):
    assert ammonia_synthesis.compositions == [Composition("N2"),
                                              Composition("H2"),
                                              Composition("NH3")
                                              ]


def test_chemical_system(ammonia_synthesis):
    assert ammonia_synthesis.chemical_system == {"H", "N"}


def test_equation(ammonia_synthesis):
    assert ammonia_synthesis.equation == "1 N2 + 3 H2 -> 2 H3N"


