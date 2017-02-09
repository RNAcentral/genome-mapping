import re
import operator as op
import itertools as it
import collections as coll

import gffutils as gff
from intervaltree import Interval
from intervaltree import IntervalTree

from genome_mapping.data import Hit
from genome_mapping.data import Comparision


def feature_upi(feature):
    return feature.attributes['Name'][0]


def upi_of(data):
    if isinstance(data, Interval):
        return upi_of(data.data)
    if isinstance(data, Hit):
        return data.upi
    if isinstance(data, gff.Feature):
        return feature_upi(data)
    raise ValueError("No way to get UPI")


class Tree(object):
    def __init__(self, filename):
        self.db = gff.create_db(filename, ':memory:')
        self.trees = coll.defaultdict(list)
        self.locations = coll.defaultdict(list)
        for (start, stop), feature in self.intervals():
            self.trees[feature.seqid].append(Interval(start, stop, feature))
            for name in feature.attributes['Name']:
                name = re.sub('_\d+$', '', name)
                self.locations[name].append(feature)

        for chromosome, intervals in self.trees.items():
            self.trees[chromosome] = IntervalTree(intervals)

    def intervals(self):
        def as_key(feature):
            getter = op.attrgetter('seqid', 'start', 'end')
            return tuple([feature_upi(feature), getter(feature)])

        seen = set()
        for feature in self.db.all_features():
            # Simplify intervals to only the most interesting ones.
            # In cases where a transcript has a *single* exon that has the same
            # endpoints we don't need the exon by itself, we only want the
            # transcript then. The ranges are the same and the sequence the
            # same so using the data isn't helpful.
            key = as_key(feature)
            if key in seen:
                continue
            seen.add(key)
            yield (feature.start, feature.stop), feature

    def search(self, start, stop):
        return {i.data for i in self.tree.search(start, stop)}

    def find(self, upi):
        upi = re.sub('_\d+$', '', upi)
        return self.locations[upi]

    def compare_to_known(self, hits, reduce_duplicates=True,
                         ignore_missing_chromosome=True):
        seen = set()
        compared = []
        for hit in hits:
            tree = self.trees[hit.chromosome]
            if not tree:
                if ignore_missing_chromosome:
                    continue
                raise ValueError("No tree for chromosome %s" % hit.chromosome)
            intervals = tree.search(hit.start, hit.stop)
            if not intervals:
                compared.append(Comparision.build(hit, None))
                continue

            if reduce_duplicates and len(intervals) > 1:
                upi = upi_of(hit)
                seen.update(intervals)
                limited = [i for i in intervals if upi_of(i) == upi]
                if limited:
                    intervals = limited

            for interval in intervals:
                feature = interval.data
                compared.append(Comparision.build(hit, feature))
                seen.add(interval)

        rest = it.chain.from_iterable(self.trees.itervalues())
        rest = it.ifilter(lambda f: f not in seen, rest)
        rest = it.imap(lambda i: Comparision.build(None, i.data), rest)
        compared.extend(rest)
        return compared

    def best_hits_within(self, hits, max_range):
        # For each feature search for all hits within max_range of the feature.
        # Then select only the 'best' hits, ie max similarity/completeness
        comparisions = []
        hit_trees = self.__hits_to_tree__(hits)
        found = set()
        for chromosome, tree in self.trees.iteritems():
            hit_tree = hit_trees[chromosome]
            for interval in tree:
                feature = interval.data
                start = interval.begin - max_range
                end = interval.end + max_range
                hits = [i.data for i in hit_tree.search(start, end)]
                if not hits:
                    continue

                similarity = max(h.stats.identity for h in hits)
                best = (h for h in hits if h.stats.identity == similarity)
                best = [Comparision.build(h, feature) for h in best]
                if best:
                    found.add(feature)
                comparisions.extend(best)

        for tree in self.trees.itervalues():
            for interval in tree:
                if interval.data not in found:
                    comparisions.append(Comparision.build(None, feature))

        return comparisions

    def __hits_to_tree__(self, hits):
        as_tree = IntervalTree.from_tuples
        by_chromosome = sorted(hits, key=op.attrgetter('chromosome'))
        by_chromosome = it.imap(lambda m: (m.start, m.stop, m), by_chromosome)
        by_chromosome = it.groupby(by_chromosome, lambda i: i[2].chromosome)
        return {k: as_tree(v) for (k, v) in by_chromosome}
