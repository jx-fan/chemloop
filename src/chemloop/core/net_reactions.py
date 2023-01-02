from collections import Set
from dataclasses import dataclass

from pymatgen.core import Composition


@dataclass(frozen=True)
class NetReaction:
    oxidant: Composition
    reducing_agent: Composition
    products: list[Composition]
    coefficients: list  # coefficients for reactants are negative. e.g. N2 + 3 H2 -> 2 NH3, coefficients: [-1, -3, 2]

    @property
    def compositions(self) -> list[Composition]:
        return [self.oxidant, self.reducing_agent, *self.products]

    @property
    def chemical_system(self) -> Set[str]:
        """Returns the chemical system as set of elements"""
        chemical_system = set()
        for c in self.compositions:
            chemical_system.update([e.symbol for e in c.elements])
        return chemical_system

    @property
    def equation(self) -> str:
        reactant_str = [self.oxidant.reduced_formula, self.reducing_agent.reduced_formula]
        product_str = [p.reduced_formula for p in self.products]
        equation_str = []
        for coefficient, formula in zip(self.coefficients, reactant_str + product_str):
            equation_str.append(str(abs(coefficient)) + " " + formula)
        return " + ".join(equation_str[:2]) + " -> " + " + ".join(equation_str[2:])
