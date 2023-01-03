from dataclasses import dataclass

from pymatgen.core import Composition


@dataclass(frozen=True)
class NetReaction:
    oxidant: Composition
    reducing_agent: Composition
    products: list[Composition]
    coefficients: list  # coefficients should be provided in following order: [oxidant, reducing_agent, products].

    @property
    def compositions(self) -> list[Composition]:
        return [self.oxidant, self.reducing_agent, *self.products]

    @property
    def chemical_system(self) -> set[str]:
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
            if abs(coefficient) != 1:
                equation_str.append(str(abs(coefficient)) + " " + formula)
            else:
                equation_str.append(formula)
        return " + ".join(equation_str[:2]) + " -> " + " + ".join(equation_str[2:])

    @property
    def balanced(self) -> bool:
        def tot_element_number(element, materials_list, coefficients):
            return sum(m.get(element) * abs(c) for m, c in zip(materials_list, coefficients))

        if all(tot_element_number(element, [self.oxidant, self.reducing_agent],
                                  self.coefficients[:2]) == tot_element_number(element,
                                                                               self.products, self.coefficients[2:])
               for
               element in self.chemical_system):
            return True
        else:
            return False
