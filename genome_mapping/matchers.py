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
import sys
import operator as op
import itertools as it

from genome_mapping import utils as ut


def known():
    return ut.names_of_children(sys.modules[__name__], Base)


def fetch(name):
    return ut.get_child(sys.modules[__name__], Base, name)


class Base(object):
    """This is the base class that all other mappers should inherit from. It
    does not contain any logic to detect if a match is valid or not, but
    provides other functionality.
    """
    __metaclass__ = abc.ABCMeta

    @abc.abstractproperty
    def name(self):
        pass

    @abc.abstractmethod
    def is_valid_hit(self, hit):
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

    def filter_matches(self, hits):
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

        for hit in hits:
            if self.is_valid_hit(hit):
                yield hit


class ExactMappingFilter(Base):
    """This is a simple filter which requires that the mapping have the same
    length in the target and query sequences and there be no gaps in either
    sequence. This is an 'exact' match.
    """

    name = 'exact'

    def is_valid_hit(self, hit):
        return hit.stats.completeness.query == 1 and \
            hit.stats.identical == hit.stats.length.query


class PercentIdentityFilter(Base):
    name = 'identity'

    def __init__(self, min=100.0, max=100.0, completeness=100.0):
        self.min = float(min)
        self.max = float(max)
        self.completeness = float(completeness)

    def identity(self, hit):
        return 100 * float(hit.stats.identical) / hit.stats.length.query

    def is_valid_hit(self, hit):
        if (100 * hit.stats.completeness.query) < self.completeness:
            return False
        query_identity = self.identity(hit)
        assert 0 <= query_identity <= 100, "Query percent %f too large" % query_identity
        return self.min <= query_identity <= self.max


class HighestIdentityFilter(PercentIdentityFilter):
    name = 'best-match'

    def __init__(self, max_hits=1, **kwargs):
        self.max_hits = max_hits
        super(HighestIdentityFilter, self).__init__(**kwargs)

    def best_in_group(self, group):
        ordered = sorted(group, key=self.identity, reverse=True)
        best = ordered[0]
        top_hits = it.takewhile(lambda h: self.identity(h) == best, ordered)
        top_hits = list(top_hits)
        if len(top_hits) < self.max_hits:
            offset = len(top_hits)
            count = self.max_hits - len(top_hits)
            top_hits.extend(ordered[offset:count])
        return top_hits

    def filter_matches(self, hits):
        grouped = it.ifilter(self.is_valid_hit, hits)
        grouped = sorted(hits, key=op.attrgetter('urs'))
        grouped = it.groupby(hits, op.attrgetter('urs'))
        for _, group in grouped:
            for hit in self.best_in_group(group):
                yield hit


class PassThroughFilter(Base):
    """Always accept all matches.
    """
    name = 'passthrough'

    def is_valid_hit(self, mapping):
        return True
