from chemloop.core.redox_materials import RedoxMaterialsSet
from chemloop.core.net_reactions import NetReaction
from rxn_network.reactions.basic import BasicReaction
from rxn_network.reactions.computed import ComputedReaction
from rxn_network.core.reaction import Reaction
from typing import List, Union, Dict, Optional, Tuple, Set
from pymatgen.core import Composition
from abc import ABCMeta, abstractmethod


class AbstractChemicalLoop(metaclass=ABCMeta):
    """
    Base definition for a chemical loop process.
    """
    @property
    @abstractmethod
    def subreactions(self) -> Tuple[Union[BasicReaction, ComputedReaction]]:
        """
        subreactions of a chemical loop process.
        """

    @property
    @abstractmethod
    def chemical_system(self) -> Set[str]:
        """
        chemical system for the looping process.
        """


class ChemicalLoopTwoStep(AbstractChemicalLoop):
    """
    Define a 2-step chemical looping process.
    MX1 + Red -> MX2 + ...
    MX2 + Oxd -> MX1 + ...
    """
    def __init__(self,
                 redox_pair: RedoxMaterialsSet,
                 net_rxn: NetReaction,
                 rxn_materials: Optional[List] = None
                 ):
        """
        A chemical looping process is defined by two redox materials with shared cation(s) and the net reaction.
        Args:
            redox_pair: Set of redox materials. Has to contain 2 materials with the same cation group.
            net_rxn: The net reaction after completing one loop.
        """
        if len(redox_pair) != 2:
            raise ValueError("Number of redox materials is not equal to 2. Current size: ", len(rp))
        self._redox_pair = redox_pair
        self._net_rxn = net_rxn
        self._subreactions = self.set_subreactions(rxn_materials)
        # self.temp = temp

    def set_subreactions(self,
                         rxn_materials: Optional[List] = None
                         ) -> Tuple[BasicReaction, BasicReaction]:
        """

        Args:
            rxn_materials: Dict contains the oxidised and the reduced redox materials.

        Returns:
            Subreactions in the chemical looping process.
        """
        redox_materials = {}
        if not rxn_materials:
            m1 = list(self._redox_pair)[0]
            m2 = list(self._redox_pair)[1]
            m1_oxi_states = m1.oxi_state_guesses()
            m2_oxi_states = m2.oxi_state_guesses()
            if not m1_oxi_states or not m2_oxi_states:
                raise ValueError("No oxidation information for redox materials. Check input.")
            if all([m1.oxi_state_guesses()[0][c.symbol] >= m2.oxi_state_guesses()[0][c.symbol]
                    for c in self._redox_pair.cations]):
                redox_materials = {"oxidised": m1, "reduced": m2}
            elif all([m1.oxi_state_guesses()[0][c.symbol] <= m2.oxi_state_guesses()[0][c.symbol]
                      for c in self._redox_pair.cations]):
                redox_materials = {"oxidised": m2, "reduced": m1}
            rxn1 = BasicReaction.balance(reactants=[self._net_rxn.reducing_agent, redox_materials["oxidised"]],
                                         products=[*self._net_rxn.products, redox_materials["reduced"]])
            rxn2 = BasicReaction.balance(reactants=[self._net_rxn.oxidant, redox_materials["reduced"]],
                                         products=[redox_materials["oxidised"]])
        else:
            rxn1 = BasicReaction.balance(reactants=[Composition(c) for c in rxn_materials[0]["reactants"]],
                                         products=[Composition(c) for c in rxn_materials[0]["products"]])
            rxn2 = BasicReaction.balance(reactants=[Composition(c) for c in rxn_materials[1]["reactants"]],
                                         products=[Composition(c) for c in rxn_materials[1]["products"]])
        return rxn1, rxn2

    @property
    def subreactions(self) -> Tuple[BasicReaction, BasicReaction]:
        return self._subreactions

    @property
    def number_of_steps(self) -> int:
        return len(self._subreactions)

    @property
    def chemical_system(self) -> Set[str]:
        """
        Returns the chemical system of all elements involved in this chemical looping process.
        """
        return self._redox_pair.chemical_system | self._net_rxn.chemical_system

    @property
    def redox_pair(self) -> RedoxMaterialsSet:
        """
        Returns the redox pair.
        """
        return self._redox_pair

    @property
    def all_materials_set(self) -> Set[str]:
        mat_set = set()
        for subrxn in self.subreactions:
            mat_set.update(set(subrxn.reactants+subrxn.products))
        return mat_set

    @property
    def net_rxn(self) -> NetReaction:
        """
        Returns the net reaction.
        """
        return self._net_rxn

    @staticmethod
    def get_avrg_valences(material: Composition):
        return material.oxi_state_guesses()

    def __str__(self):
        rxns = [rxn.__str__() for rxn in self.subreactions]
        rxns.append("="*len(max(rxns, key=len)))
        rxns.append("Net rxn: " + str(self._net_rxn))
        return '\n'.join(rxns)


class ChemicalLoopThreeStep(AbstractChemicalLoop):
    # TODO: implement this class for three step chemical looping process
    @property
    def subreactions(self) -> Tuple[Union[BasicReaction, ComputedReaction]]:
        pass

    @property
    def chemical_system(self) -> List[str]:
        pass
