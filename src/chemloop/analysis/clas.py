"""
Analysis of reaction pathways in Chemical Looping Ammonia Synthesis (CLAS)
"""
from pathlib import Path

import pandas as pd
from monty.serialization import loadfn
from pymatgen.core.composition import Element, Composition
from rxn_network.pathways.balanced import BalancedPathway
from rxn_network.pathways.pathway_set import PathwaySet
from rxn_network.reactions.basic import BasicReaction


class AnalyseHydroPathwaySet:
    def __init__(self,
                 pathway_set: PathwaySet,
                 net_rxn: BasicReaction,
                 net_rxn_energy: float,
                 nitride: Composition,
                 oxide: Composition,
                 cost_method: str = "arithmetic",
                 max_combo: int = 5
                 ):
        """

        Args:
            pathway_set:
            net_rxn:
            net_rxn_energy:
            nitride:
            oxide:
            cost_method:
            max_combo:
        """
        self.pathway_set = pathway_set
        self.net_reaction = net_rxn
        self.net_reaction_energy = net_rxn_energy
        self.nitride = nitride
        self.oxide = oxide
        self.cost_method = cost_method
        self.max_combo = max_combo

    @staticmethod
    def arithmetic_cost(costs) -> float:
        """
        Calculate the arithmetic average of the costs
        Args:
            costs:

        Returns:

        """
        return round(sum(costs) / len(costs), 3)

    @property
    def paths(self) -> list[BalancedPathway]:
        """

        Returns:

        """
        default_paths = [path for path in self.pathway_set.get_paths() if len(path.costs) <= self.max_combo]
        if self.cost_method == "mcdermott":
            return default_paths
        elif self.cost_method == "arithmetic":
            return sorted(default_paths, key=lambda p: self.arithmetic_cost(p.costs))

    @property
    def lowest_cost_pathway(self) -> BalancedPathway:
        """

        Returns:

        """
        return self.paths[0]

    @property
    def lowest_cost(self) -> float:
        """

        Returns:

        """
        if self.cost_method == "mcdermott":
            return self.lowest_cost_pathway.average_cost
        elif self.cost_method == "arithmetic":
            return self.arithmetic_cost(self.lowest_cost_pathway.costs)

    @property
    def cations(self) -> list[Element]:
        """

        Returns:

        """
        try:
            cations = [e for e in self.nitride if e.symbol not in ["O", "N"]]
        except:
            try:
                cations = [e for e in self.oxide if e.symbol not in ["O", "N"]]
            except:
                raise ValueError("something is wrong here.")
        return cations

    @classmethod
    def from_file(cls,
                  file_pathway: Path,
                  file_energy: Path,
                  rxn_column_name: str,
                  e_column_name: str,
                  normalize_to: str,
                  cost_method: str = "arithmetic",
                  max_combo: int = 5
                  ) -> "AnalyseHydroPathwaySet":
        """
        Customised method for loading pre-calculated data.
        Args:
            max_combo:
            cost_method:
            normalize_to:
            file_pathway:
            file_energy:
            rxn_column_name:
            e_column_name:

        Returns:

        """
        oxide, nitride = [Composition(formula) for formula in file_pathway.name.split("_")[:2]]
        df = pd.read_csv(file_energy, sep=";", index_col=[0, 1],
                         na_values='', keep_default_na=False,  # avoid filtering out sodium nitride (NaN)
                         )
        net_rxn = BasicReaction.from_string(df.loc[(oxide.reduced_formula, nitride.reduced_formula),
                                                   rxn_column_name])
        energy = df.loc[(oxide.reduced_formula, nitride.reduced_formula), e_column_name]
        normalized_net_rxn = net_rxn.normalize_to(Composition(normalize_to))
        factor = normalized_net_rxn.num_atoms / net_rxn.num_atoms
        return cls(pathway_set=loadfn(file_pathway),
                   net_rxn=normalized_net_rxn,
                   net_rxn_energy=energy * factor,
                   nitride=nitride,
                   oxide=oxide,
                   cost_method=cost_method,
                   max_combo=max_combo
                   )


def ammonia_yield_energy(pathway: BalancedPathway,
                         normalise_to_per_ammonia: bool = True
                         ) -> float:
    """

    Args:
        pathway:
        normalise_to_per_ammonia:

    Returns:

    """
    ammonia_steps = []
    for r in pathway.reactions:
        if Composition("NH3") in r.products:
            ammonia_steps.append(r)
    if normalise_to_per_ammonia:
        ammonia_steps = [step.normalize_to(Composition("NH3")) for step in ammonia_steps]
    return sum([r.energy for r in ammonia_steps]) / len(ammonia_steps)
