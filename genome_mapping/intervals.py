import operator as op
import itertools as it
import collections as coll

import gffutils as gff
from intervaltree import IntervalTree

from genome_mapping.data import urs_of
from genome_mapping.data import Comparision
from genome_mapping.data import FeatureData


class Tree(object):
    def __init__(self, filename):
        self.db = gff.create_db(filename, ':memory:')
        self.trees = self.__build_tree__(self.intervals())

    def intervals(self):
        def as_key(feature):
            getter = op.attrgetter('seqid', 'start', 'end')
            return tuple([urs_of(feature), getter(feature)])

        grouped = coll.defaultdict(list)
        for feature in self.db.all_features(featuretype='noncoding_exon'):
            # Group the features by their parent. This is required to merge the
            # gff3 exons into a single unified feature. This way we can
            # correclty detect if something is spliced or not. Without this we
            # end up with many single exons.
            parent = feature.attributes['Parent'][0]
            grouped[parent].append(feature)

        for subfeatures in grouped.itervalues():
            yield FeatureData.build(subfeatures)

    def search(self, start, stop):
        return {i.data for i in self.tree.search(start, stop)}

    def compare_to_known(self, hits, reduce_duplicates=True,
                         ignore_missing_chromosome=True):
        seen = set()
        compared = []
        for hit in hits:
            if hit.chromosome not in self.trees:
                if ignore_missing_chromosome:
                    continue
                raise ValueError("No tree for chromosome %s" % hit.chromosome)

            tree = self.trees[hit.chromosome]
            intervals = tree.search(hit.start, hit.stop)
            if not intervals:
                compared.append(Comparision.build(hit, None))
                continue

            if reduce_duplicates and len(intervals) > 1:
                urs = urs_of(hit)
                seen.update(intervals)
                limited = [i for i in intervals if urs_of(i) == urs]
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
        hit_trees = self.__build_tree__(hits)
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

    def __build_tree__(self, data):
        as_tree = IntervalTree.from_tuples
        by_chromosome = sorted(data, key=op.attrgetter('chromosome'))
        by_chromosome = it.imap(lambda m: (m.start, m.stop, m), by_chromosome)
        by_chromosome = it.groupby(by_chromosome, lambda i: i[2].chromosome)
        return {k: as_tree(v) for (k, v) in by_chromosome}
