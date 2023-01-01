import warnings
from typing import List, Union, Optional, Dict, Set

import numpy as np
from pymatgen.core import Composition
from rxn_network.reactions.basic import BasicReaction


class NetReaction(BasicReaction):
    """
    Define the net reaction for a chemical looping process. For instance, in the case of chemical looping ammonia
    synthesis, the net reaction is:
        N2 + 3 H2 --> 2 NH3
    """

    def __init__(self,
                 oxidant: Composition,
                 reducing_agent: Composition,
                 products: List[Composition],
                 coefficients: Union[np.ndarray, List[float]],
                 lowest_num_errors: Union[int, float] = 0,
                 ):
        """
        A NetReaction object is defined as a redox reaction containing oxidant, reducing agent, and products. The
        chemical reaction coefficients can be automatically determined or manually supplied.
        Args:
            oxidant: The Composition of the oxidant.
            reducing_agent: The Composition of the reducing agent.
            products: List of Composition objects in the final products
            coefficients: List of coefficients, where negative coeff distinguishes a
                reactant.
            lowest_num_errors: the minimum number of errors reported by the reaction
                balancing algorithm (see the balance() method). A number of errors
                >= 1 means that the reaction may be different than intended (some
                phases may be shuffled or removed entirely).
        """
        self._oxidant = oxidant
        self._reducing_agent = reducing_agent
        self._products = products
        self._coefficients = coefficients
        self._compositions = [self._oxidant, self._reducing_agent, *self._products]
        super().__init__(compositions=self._compositions,
                         coefficients=self._coefficients,
                         lowest_num_errors=lowest_num_errors
                         )

    @classmethod
    def balance(
            cls,
            oxidant: Composition,
            reducing_agent: Composition,
            products: List[Composition],
            data: Optional[Dict] = None
    ) -> "NetReaction":
        """
        Balance the reaction. Oxidant and reducing agent to be specified as pymatgen.core.Composition. Products to be
        input as list of Composition object.
        Args:
            oxidant: The Composition of the oxidant.
            reducing_agent: The Composition of the reducing agent.
            products: List of Composition objects in the final products
            data: Optional dictionary containing extra data about the reaction.
        """
        coefficients, lowest_num_errors, _ = cls._balance_coeffs([oxidant, reducing_agent], products)
        if any([abs(c) < 1e-10 for c in coefficients]):
            raise ValueError("Reaction can not be balanced! Please check the input of reactants and products.")
        if any([c > 0 for c in coefficients[:2]]) or any([c < 0 for c in coefficients[-2:]]):
            warnings.warn("Some of the compounds in the reaction may have swapped side.")
        return cls(
            oxidant=oxidant,
            reducing_agent=reducing_agent,
            products=products,
            coefficients=coefficients,
            lowest_num_errors=lowest_num_errors
        )

    @property
    def oxidant(self):
        """
        The oxidant in this redox reaction.
        """
        return self._oxidant

    @property
    def reducing_agent(self):
        """
        The reducing agent in this redox reaction.
        """
        return self._reducing_agent

    @property
    def reactants(self) -> List[Composition]:
        """
        List of reactants in this redox reaction.
        """
        return [self.oxidant, self._reducing_agent]

    @property
    def products(self) -> List[Composition]:
        """
        List of products in this redox reaction.
        """
        return self._products

    @property
    def coefficients(self) -> np.ndarray:
        """
        Array of reaction coefficients.
        """
        return self._coefficients

    @property
    def energy(self) -> float:
        """The energy of this reaction"""
        raise ValueError("No energy for a net reaction!")

    @property
    def chemical_system(self) -> Set[str]:
        """Returns the chemical system as set of elements"""
        return {e.symbol for e in self.elements}
