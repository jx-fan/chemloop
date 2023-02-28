from abc import ABCMeta, abstractmethod

from rxn_network.pathways.pathway_set import PathwaySet
from rxn_network.reactions.basic import BasicReaction


class AbstractPathwayFilter(metaclass=ABCMeta):
    """
    Base definition for a pathway filter.
    """

    @abstractmethod
    def filter(self,
               pathway: PathwaySet
               ) -> PathwaySet:
        """
        Apply the filter
        Returns:

        """


class ReactionFilter(AbstractPathwayFilter):
    """
    Allow additional criteria on whether certain reactions can exist.
    """

    def __init__(self,
                 rxns_to_include: set[BasicReaction],
                 rxns_to_exclude: set[BasicReaction] = None,
                 filter_reactants_only: bool = True
                 ):
        """

        Args:
            rxns_to_include:
            rxns_to_exclude:
            filter_reactants_only:
        """
        self.rxns_to_include = rxns_to_include
        self.rxns_to_exclude = rxns_to_exclude if rxns_to_exclude else set()
        self.filter_reactants_only = filter_reactants_only

    def filter(self,
               pathway: PathwaySet
               ) -> PathwaySet:
        """

        Args:
            pathway:

        Returns:

        """
        total_paths = pathway.get_paths()
        filtered_paths = []
        for p in total_paths:
            # only compare if reactants and products are matching to avoid: A + B -> C != 2 A + 2 B -> 2 C
            total_compositions = [(set(r.reactants), set(r.products)) for r in p.reactions]
            if self.filter_reactants_only:
                if any([set(r.reactants) in [c[0] for c in total_compositions]
                        for r in self.rxns_to_exclude]):
                    continue
                elif any([set(r.reactants) in [c[0] for c in total_compositions]
                          for r in self.rxns_to_include]):
                    filtered_paths.append(p)
            else:
                if any([(set(r.reactants), set(r.products)) in total_compositions
                        for r in self.rxns_to_exclude]):
                    continue
                if any([(set(r.reactants), set(r.products)) in total_compositions
                        for r in self.rxns_to_include]):
                    filtered_paths.append(p)
        return PathwaySet.from_paths(filtered_paths)
