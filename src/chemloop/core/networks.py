from abc import ABCMeta, abstractmethod
from typing import List, Tuple, Union

from chemloop.core.chemical_loops import AbstractChemicalLoop, ChemicalLoopTwoStep, ChemicalLoopThreeStep
from chemloop.utils.mp_entries import get_entries_from_json, get_entries_from_api
from rxn_network.core.enumerator import Enumerator
from rxn_network.costs.softplus import Softplus
from rxn_network.entries.entry_set import GibbsEntrySet
from rxn_network.network.network import ReactionNetwork
from rxn_network.pathways.balanced import BalancedPathway
from rxn_network.pathways.pathway_set import PathwaySet
from rxn_network.pathways.solver import PathwaySolver
from rxn_network.reactions.computed import ComputedReaction


class AbstractNetwork(metaclass=ABCMeta):
    _chemical_loop = AbstractChemicalLoop

    @abstractmethod
    def build_networks(self) -> List[ReactionNetwork]:
        """
        Build reaction network for the chemical looping process
        """

    @abstractmethod
    def solve_pathways(self) -> List[BalancedPathway]:
        """
        Find the pathways for each subreactions
        """

    @property
    @abstractmethod
    def chemical_system(self):
        """
        Returns the chemical system.
        """
        return self._chemical_loop.chemical_system


class ChemicalLoopNetwork(AbstractNetwork):
    """
    The object analysing the reaction networks in each chemical loop process steps.
    """

    def __init__(self,
                 chemical_loop: Union[ChemicalLoopTwoStep, ChemicalLoopThreeStep],
                 temps: List[float]
                 ):
        """
        Supply the type of chemical looping process and the temperatures each subreaction is carried out. Ensure the
        number of temperatures imported is the same as the number of steps in the chemical looping process.
        Args:
            chemical_loop: Chemical loop object defining the redox materials and the reactions.
            temps: List of float numbers defining the reaction temperatures in Kelvin.
        """
        self.chemical_loop = chemical_loop
        self.temps = temps
        if self.chemical_loop.number_of_steps != len(self.temps):
            raise ValueError("The input of temperatures do not match the number of steps. Please double check")
        self._networks = None

    def _get_gibbs_entry_set(self,
                             temp: float,
                             e_above_hull: float,
                             from_json=True,
                             include_compounds=True,
                             ) -> GibbsEntrySet:
        """
        Internal method to get the gibbs entry set.
        Args:
            temp: Reaction temperature in K.
            e_above_hull: Energy above hull used to filter out unstable materials in the GibbsEntrySet
            from_json: Whether or not to read the materials data from json files.

        Returns:
            GibbsEntrySet object for each subreaction.

        """
        if from_json:
            cations = [e.symbol for e in self.chemical_loop.redox_pair.cations]
            entries = get_entries_from_json(cations_in_redox_materials=sorted(cations),
                                            chemsys_net_rxn={cation for cation in self.chemical_loop.chemical_system
                                                             if cation not in cations})
        else:
            entries = get_entries_from_api()
        entry_set = GibbsEntrySet.from_computed_entries(entries, temp, include_nist_data=True)
        if e_above_hull >= 0:
            filtered_entry_set = entry_set.filter_by_stability(e_above_hull)
            if include_compounds:
                for c in list(self.chemical_loop.all_materials_set):
                    filtered_entry_set.add(
                        entry_set.get_min_entry_by_formula(c.reduced_formula))  # make sure materials in ReactionNetwork
        else:
            raise ValueError("Wrong hull energy: ", e_above_hull)
        return filtered_entry_set

    def entry_sets(self, e_above_hull: float) -> List[GibbsEntrySet]:
        """
        Get a list of materials entries for all subreactions in the chemical looping process.
        Args:
            e_above_hull:

        Returns:
            A list of GibbsEntrySet for constructing the reaction networks.

        """
        return [self._get_gibbs_entry_set(temp, e_above_hull=e_above_hull) for temp in self.temps]

    def build_networks(self,
                       enumerator: Enumerator,
                       e_above_hull: float = 0.01) -> List[ReactionNetwork]:
        """
        Construct the list of reaction network graph objects for the chemical looping process.

        Args:
            e_above_hull:
            enumerator:

        Returns:
            None

        """
        networks = []
        for entries, temp, subrxn in zip(self.entry_sets(e_above_hull), self.temps, self.chemical_loop.subreactions):
            rxns = enumerator.enumerate(entries)
            rn = ReactionNetwork(rxns, Softplus(temp))
            print("Building network for: ", subrxn)
            rn.build()
            networks.append(rn)
        self._networks = networks

    def solve_pathways(self, k=5, max_num_combos=5) -> Tuple[List[ComputedReaction], List[PathwaySet]]:
        networks = self._networks
        subreactions = self.chemical_loop.subreactions
        computed_subrxns = []
        balanced_paths_cl = []
        for network, subrxn, temp in zip(networks,
                                         subreactions,
                                         self.temps):
            network.set_precursors(subrxn.reactants)
            for product in subrxn.products:
                if product in self.chemical_loop.redox_pair:
                    network.set_target(product)
            paths = network.find_pathways(subrxn.products, k=k)
            cf = Softplus(temp)
            ps = PathwaySolver(paths, network.entries, cf)
            entries = network.entries
            computed_subrxn = ComputedReaction.balance([entries.get_min_entry_by_formula(c.reduced_formula)
                                                        for c in subrxn.reactants],
                                                       [entries.get_min_entry_by_formula(c.reduced_formula)
                                                        for c in subrxn.products])
            print("Solving pathways for: ", computed_subrxn)
            balanced_paths = ps.solve(computed_subrxn, max_num_combos=max_num_combos,
                                      intermediate_rxn_energy_cutoff=0.0,
                                      use_minimize_enumerator=True,
                                      filter_interdependent=True)
            computed_subrxns.append(computed_subrxn)
            balanced_paths_cl.append(balanced_paths)
        return computed_subrxns, balanced_paths_cl

    @property
    def chemical_system(self):
        pass


def save_balanced_paths(subrxns, paths_list, filename="results.log", all_paths=True):
    with open(filename, "w") as f:
        for subrxn, paths in zip(subrxns, paths_list):
            f.write("=" * 30 + "\n")
            f.write(str(subrxn) + " (dG = " + str(round(subrxn.energy_per_atom, 3)) + " eV/atom)\n")
            f.write("=" * 30 + "\n")
            if all_paths:
                for idx, path in enumerate(paths):
                    f.write(f"Path {idx + 1}\n")
                    f.write(str(path))
                    f.write("\n\n")
            else:
                f.write("Path 1\n")
                f.write(str(path))
                f.write("\n\n")
