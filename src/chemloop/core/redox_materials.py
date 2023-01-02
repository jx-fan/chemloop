from collections.abc import MutableSet
from pathlib import Path
from typing import Set, Iterable

from pymatgen.core import Composition
from pymatgen.core.periodic_table import Element

ENTRIES_DATA_PATH = Path("entries_data")


class RedoxMaterialsSet(MutableSet):
    """
    A container for manipulating redox materials.
    """

    def __init__(self, materials: Iterable[Composition]):
        """
        The supplied collection of Composition objects will automatically be converted to a set of unique items.
        Args:
            materials: A collection of Composition entries that will make up the redox materials set
        """
        self.materials = set(materials)
        if not set.intersection(*[self._get_element_set(m, "cation") for m in self.materials]):
            raise ValueError("Wrong inputs for the set of redox materials: cations are not shared in all materials.")

    def __contains__(self, item):
        return item in self.materials

    def __iter__(self):
        return self.materials.__iter__()

    def add(self, material: Composition):
        """
        Add a material to the set.

        Args:
            material: material composition.

        Returns:

        """
        self.materials.add(material)
        self._clear_cache()

    def discard(self, material: Composition):
        """
        Discard a material.

        Args:
            material: material composition.

        Returns:

        """
        self.materials.discard(material)
        self._clear__cache()

    @staticmethod
    def _get_element_set(material: Composition, ele_type: str) -> Set[Element]:
        """
        Helper method to return a set of Element objects based on the type (cation or anion) from a material.
        Args:
            material: The composition of a material
            ele_type: Returns either the cations or the anions of the material

        Returns:
            A set of Elements based on the element type.

        """
        cations: Set[Element] = set()
        anions: Set[Element] = set()
        for e in material.elements:
            if e.symbol in ["N", "O", "H", "F", "S"]:  # common CL processes
                anions.add(e)
            else:
                cations.add(e)
        if ele_type == "cation":
            return cations
        elif ele_type == "anion":
            return anions
        else:
            raise ValueError("Wrong element type: ", ele_type)

    def __len__(self):
        return len(self.materials)

    @property
    def num_materials(self) -> int:
        """
        Returns the number of materials in the set of redox materials.
        """
        return self.__len__()

    @property
    def cations(self) -> Set[Element]:
        """
        Returns the shared cations in the provided redox materials set. e.g. {"Li", "Fe"}
        """
        return self._get_element_set(list(self.materials)[0], ele_type="cation")

    @property
    def anions(self) -> Set[Element]:
        """
        Returns the set of anions in the provided redox materials set. e.g. {"O", "N"}
        """
        anions = set()
        for material in self.materials:
            anions.update(self._get_element_set(material, ele_type="anion"))
        return anions

    @property
    def chemical_system(self) -> Set[str]:
        """
        Returns a set representing the chemical system, e.g. {"Co", "Mo", "O", "N"}
        """
        return {ele.symbol for ele in self.cations | self.anions}

    def add_charges(self) -> "RedoxMaterialsSet":
        """
        Returns a RedoxMaterialsSet that consists of Composition objects containing guessed oxidation states.
        """
        return RedoxMaterialsSet([m.add_charges_from_oxi_state_guesses() for m in self.materials])

    def _clear_cache(self):
        """
        Clears cached properties.
        """
        try:
            del self.entries_list
        except AttributeError:
            pass

        try:
            del self.pd_dict
        except AttributeError:
            pass
