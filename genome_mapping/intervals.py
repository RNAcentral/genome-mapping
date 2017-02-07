import re
import operator as op
import itertools as it
import collections as coll

import gffutils as gff
from intervaltree import Interval
from intervaltree import IntervalTree

from genome_mapping.data import Comparision


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
        seen = set()
        for feature in self.db.all_features():
            if feature.featuretype != 'transcript':
                continue

            name = feature.attributes['Name'][0]
            key = (name, feature.start, feature.end)
            if key in seen:
                continue
            seen.add(key)
            yield (feature.start, feature.stop), feature

    def search(self, start, stop):
        return {i.data for i in self.tree.search(start, stop)}

    def find(self, upi):
        upi = re.sub('_\d+$', '', upi)
        return self.locations[upi]

    def compare_to_known(self, matches):
        as_tree = IntervalTree.from_tuples
        by_chromosome = sorted(matches, key=op.attrgetter('chromosome'))
        by_chromosome = it.imap(lambda m: (m.start, m.stop, m), by_chromosome)
        by_chromosome = it.groupby(by_chromosome, lambda i: i[2].chromosome)
        by_chromosome = {k: as_tree(v) for (k, v) in by_chromosome}

        seen = set()
        compared = []
        for name, features in self.locations.items():
            for feature in features:
                chromosome = feature.seqid
                tree = by_chromosome[chromosome]
                start, end = sorted([feature.start, feature.end])
                intervals = tree[start:end]
                if not intervals:
                    intervals = tree.search(start, end)

                if not intervals:
                    compared.append(Comparision.build(None, feature))
                    continue

                for interval in intervals:
                    match = interval.data
                    compared.append(Comparision.build(match, feature))
                    seen.add(interval)

        for chromosome, tree in by_chromosome.items():
            for interval in tree:
                if interval not in seen:
                    compared.append(Comparision.build(interval.data, None))

        return compared
