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
        self._pathway_set = pathway_set
        self._net_rxn = net_rxn
        self._net_rxn_energy = net_rxn_energy
        self._nitride = nitride
        self._oxide = oxide
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
    def net_rxn_cost(self, normalised=False) -> float:
        if normalised:
            return self.softplus(self.temperature, self.net_rxn_energy/self.net_rxn.num_atoms/96)
        else:
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
        default_paths = [path for path in self._pathway_set.get_paths() if len(path.costs) <= self.max_combo]
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
                  normalize_to: str = None,
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
        file_pathway = Path(file_pathway)
        file_energy = Path(file_energy)
        oxide, nitride = [Composition(formula) for formula in file_pathway.name.split("_")[:2]]
        df = pd.read_csv(file_energy, sep=";", index_col=[0, 1],
                         na_values='', keep_default_na=False,  # avoid filtering out sodium nitride (NaN)
                         )
        net_rxn = BasicReaction.from_string(df.loc[(oxide.reduced_formula, nitride.reduced_formula),
                                                   rxn_column_name])
        energy = df.loc[(oxide.reduced_formula, nitride.reduced_formula), e_column_name]  # eV/atom
        if normalize_to:
            normalized_net_rxn = net_rxn.normalize_to(Composition(normalize_to))
            return cls(pathway_set=loadfn(file_pathway),
                       net_rxn=normalized_net_rxn,
                       net_rxn_energy=energy * normalized_net_rxn.num_atoms * 96,  # kJ/mol molecule
                       nitride=nitride,
                       oxide=oxide,
                       cost_method=cost_method,
                       max_combo=max_combo
                       )
        else:
            return cls(pathway_set=loadfn(file_pathway),
                       net_rxn=net_rxn,
                       net_rxn_energy=energy,
                       nitride=nitride,
                       oxide=oxide,
                       cost_method=cost_method,
                       max_combo=max_combo
                       )


def ammonia_yield_energy(pathway: BalancedPathway,
                         normalise_to_per_ammonia: bool = True,
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
        return float(np.mean([r.energy for r in ammonia_steps]) * 96)  # kJ/mol NH3
    else:
        return float(np.mean([r.energy_per_atom for r in ammonia_steps]))  # eV/atom
