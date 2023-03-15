"""
Analysis of reaction pathways in Chemical Looping Ammonia Synthesis (CLAS)
"""
from pathlib import Path

import numpy as np
import pandas as pd
from monty.serialization import loadfn
from pymatgen.core.composition import Element, Composition
from rxn_network.pathways.balanced import BalancedPathway
from rxn_network.pathways.pathway_set import PathwaySet
from rxn_network.reactions.basic import BasicReaction
from rxn_network.reactions.computed import ComputedReaction

from chemloop.analysis.filter import AbstractPathwayFilter


class AnalyseHydroPathwaySet:
    def __init__(self,
                 pathway_set: PathwaySet,
                 nitride: Composition,
                 oxide: Composition,
                 net_rxn: BasicReaction = None,
                 net_rxn_energy: float = None,
                 pathway_filter: AbstractPathwayFilter = None,
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
        self._net_rxn = net_rxn
        self._net_rxn_energy = net_rxn_energy
        self._nitride = nitride
        self._oxide = oxide
        self._pathway_set = pathway_set
        self.pathway_filter = pathway_filter
        self.cost_method = cost_method
        self.max_combo = max_combo

    @staticmethod
    def softplus(t: float,
                 e: float) -> float:
        """
        Simple Softplus function that only takes temperature and reaction energy as input.
        Args:
            t:
            e:

        Returns:

        """
        return np.log(1 + (273 / t) * np.exp(e))

    @property
    def pathway_set(self):
        if self.pathway_filter:
            return self.pathway_filter.filter(self._pathway_set)
        else:
            return self._pathway_set

    @property
    def net_rxn_cost(self) -> float:
        return self.softplus(self.temperature, self.net_rxn_energy)

    @property
    def temperature(self) -> float:
        return list(self.lowest_cost_pathway.entries)[0].temperature

    @property
    def net_rxn(self) -> BasicReaction:
        return self._net_rxn

    @property
    def net_rxn_energy(self) -> float:
        return self._net_rxn_energy

    @property
    def nitride(self) -> Composition:
        return self._nitride

    @property
    def oxide(self) -> Composition:
        return self._oxide

    @property
    def paths(self) -> list[BalancedPathway]:
        """

        Returns:

        """
        default_paths = [path for path in self.pathway_set.get_paths() if len(path.costs) <= self.max_combo]
        if self.cost_method == "mcdermott":
            return default_paths
        elif self.cost_method == "arithmetic":
            return sorted(default_paths, key=lambda p: float(np.mean(p.costs)))

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
            return float(np.mean(self.lowest_cost_pathway.costs))

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
                  file_pathway: str,
                  file_energy: str,
                  rxn_column_name: str,
                  e_column_name: str,
                  cost_method: str = "arithmetic",
                  max_combo: int = 5,
                  pathway_filter: AbstractPathwayFilter = None,
                  ) -> "AnalyseHydroPathwaySet":
        """
        Customised method for loading pre-calculated data.
        Args:
            pathway_filter:
            max_combo:
            cost_method:
            normalize_to:
            file_pathway:
            file_energy:
            rxn_column_name:
            e_column_name:

        Returns:

        """
        file_pathway = Path(file_pathway)
        file_energy = Path(file_energy)
        oxide, nitride = [Composition(formula) for formula in file_pathway.name.split("_")[:2]]
        df = pd.read_csv(file_energy, sep=";", index_col=[0, 1],
                         na_values='', keep_default_na=False,  # avoid filtering out sodium nitride (NaN)
                         )
        net_rxn = BasicReaction.from_string(df.loc[(oxide.reduced_formula, nitride.reduced_formula), rxn_column_name])
        energy = df.loc[(oxide.reduced_formula, nitride.reduced_formula), e_column_name]  # eV/atom
        return cls(pathway_set=loadfn(file_pathway),
                   net_rxn=net_rxn,
                   net_rxn_energy=energy,
                   nitride=nitride,
                   oxide=oxide,
                   cost_method=cost_method,
                   max_combo=max_combo,
                   pathway_filter=pathway_filter
                   )


def balanced_reactions(pathway: BalancedPathway) -> list[ComputedReaction]:
    """
    
    Args:
        pathway:

    Returns:

    """
    return [ComputedReaction(r.entries, r.coefficients * c)
            for r, c in zip(pathway.reactions, pathway.coefficients)]


def ammonia_yield_steps(pathway: BalancedPathway,
                        net_rxn: BasicReaction,
                        normalise_to_per_ammonia: bool = True,
                        ) -> tuple[list[ComputedReaction], float]:
    """

    Args:
        net_rxn:
        pathway:
        normalise_to_per_ammonia:

    Returns:

    """
    ammonia_steps = []
    for r in balanced_reactions(pathway):
        if Composition("NH3") in r.products:
            ammonia_steps.append(r)
    if normalise_to_per_ammonia:
        factor = net_rxn.normalize_to(Composition("NH3")).coefficients[0] / net_rxn.coefficients[0]
        ammonia_steps = [ComputedReaction(r.entries, r.coefficients * factor) for r in ammonia_steps]
        return ammonia_steps, float(np.mean([r.energy for r in ammonia_steps]) * 96)  # kJ/mol NH3
    else:
        return ammonia_steps, float(np.mean([r.energy_per_atom for r in ammonia_steps]))  # eV/atom


def limiting_step(pathway: BalancedPathway,
                  net_rxn: BasicReaction,
                  normalise_to_per_ammonia: bool = True,
                  ) -> tuple[ComputedReaction, float]:
    """

    Args:
        pathway:
        net_rxn:
        normalise_to_per_ammonia:

    Returns:

    """
    sorted_reactions = sorted(balanced_reactions(pathway), key=lambda x: x.energy)
    if not normalise_to_per_ammonia:
        return sorted_reactions[-1], sorted_reactions[-1].energy  # eV/atom
    else:
        factor = net_rxn.normalize_to(Composition("NH3")).coefficients[0] / net_rxn.coefficients[0]
        limiting_rxn = ComputedReaction(sorted_reactions[-1].entries,
                                        sorted_reactions[-1].coefficients * factor)
        return limiting_rxn, limiting_rxn.energy * 96  # kJ/mol NH3
