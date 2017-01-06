"""This module contains classes for filtering matches against the genome. This
is a key step in determing the mapping of a sequence to a genome. All programs
will produce some alignment, but different RNA's from different organisms will
require different criteria for being a good match. For example, in humans long
RNA's will have to allow for splicing events and thus gaps in the genomic
sequence, while in bacteria the RNA's should match exactly with no gaps in
either sequence.

The base filter, which does not implement the selection logic is MappingFilter.
All other filters are expected to inherit from it.
"""

import abc


class MappingFilter(object):
    """This is the base class that all other mappers should inherit from. It
    does not contain any logic to detect if a match is valid or not, but
    provides other functionality.
    """
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def is_valid_match(self, mapping):
        """This is the key method all inheriting classes will have to
        implement, it should check if the given mapping is 'vaild' for some
        criteria.

        Parameters
        ----------
        mapping : Mapping
            The mapping object to filter.

        Returns
        -------
        valid : bool
            True if the given mapping is valid.
        """
        pass

    def filter_matches(self, mappings):
        """Filter all matches to only those that are valid.

        Parameters
        ----------
        mappings : list
            List of Mapping objects to filter with self.is_valid_match.

        Returns
        -------
        valid_mappings : list
            List of mappings which are valid.
        """

        return [m for m in mappings if self.is_valid_match(m)]


class ExactMappingFilter(MappingFilter):
    """This is a simple filter which requires that the mapping have the same
    length in the target and query sequences and there be no gaps in either
    sequence. This is an 'exact' match.
    """

    def is_valid_match(self, mapping):
        """Check if the mapping is exact.

        Parameters
        ----------
        mapping : Mapping
        
        Returns
        -------
        is_exact : bool
            True if the mapping is exact.
        """
        return mapping.stats.hit_len == mapping.stats.query_len and \
            mapping.identical == mapping.stats.hit_len
