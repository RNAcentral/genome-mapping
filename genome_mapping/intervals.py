import re
import operator as op
import collections as coll

import attr
import gffutils as gff
from intervaltree import Interval
from intervaltree import IntervalTree


@attr.s()
class Shift(object):
    start = attr.ib()
    stop = attr.ib()

    @classmethod
    def from_interval(cls, start, stop, interval):
        return cls(start=interval.begin - start,
                   stop=interval.end - stop)

    @classmethod
    def from_match(cls, match, feature):
        return cls(start=match.start - feature.start,
                   stop=match.stop - feature.end)

    @classmethod
    def cross_chromosome(cls):
        return cls(start=float('-inf'), stop=float('inf'))

    @property
    def total(self):
        return abs(self.start) + abs(self.stop)

    def is_exact(self):
        return not self.start and not self.stop


@attr.s()
class Overlap(object):
    feature = attr.ib(hash=False)
    location = attr.ib()
    shift = attr.ib()

    @classmethod
    def from_interval(cls, start, stop, interval):
        return cls(feature=interval.data,
                   location=None,
                   shift=Shift.from_interval(start, stop, interval))

    @classmethod
    def from_match(cls, match, feature):
        shift = Shift.from_match(match, feature)
        if match.chromosome != feature.seqid:
            shift = Shift.cross_chromosome()

        return cls(feature=feature,
                   location=match,
                   shift=shift)

    def is_exact(self):
        return self.shift.is_exact()


class OneBased(object):
    def normalize(self, start, stop):
        return (start - 1, stop)


class Tree(object):
    def __init__(self, filename, coordinates=OneBased()):
        self.db = gff.create_db(filename, ':memory:')
        self.trees = coll.defaultdict(list)
        self.locations = coll.defaultdict(list)
        for (start, stop), feature in self.features(coordinates):
            self.trees[feature.seq].append(Interval(start, stop, feature))
            for name in feature.attributes['Name']:
                name = re.sub('_\d+$', '', name)
                self.locations[name].append(feature)

        for chromosome, intervals in self.trees.items():
            self.trees[chromosome] = IntervalTree(intervals)

        # intervals = (Interval(f.start, f.end, f) for f in self.features())
        # self.tree = IntervalTree(intervals)

    def intervals(self, coordinates):
        seen = set()
        for feature in self.db.all_features():
            name = feature.attributes['Name'][0]
            key = (name, feature.start, feature.end)
            if key in seen:
                continue
            seen.add(key)

            yield coordinates.normalize(feature.start, feature.stop), feature

    def search(self, start, stop):
        matches = set()
        for interval in self.tree.search(start, stop):
            matches.add(Overlap.from_interval(start, stop, interval))
        return matches

    def find(self, upi):
        upi = re.sub('_\d+$', '', upi)
        return self.locations[upi]

    def find_correct_overlaps(self, matches):
        overlaps = []
        for match in matches:
            correct = self.find(match.name)
            if not correct:
                raise ValueError("No correct locations for %s", match.name)
            correct = [Overlap.from_match(match, c) for c in correct]
            overlaps.append(min(correct, key=lambda o: o.shift.total))
        return overlaps

    def find_all_overlaps(self, matches):
        overlaps = []
        shift = op.attrgetter('shift')
        for match in matches:
            correct = self.find(match.name)
            choosable = [Overlap.from_match(match, c) for c in correct]
            choosable.extend(self.search(match.start, match.stop))
            nearest = min(choosable, key=shift)
            overlaps.append(nearest)
        return overlaps

    def find_other_overlaps(self, matches):
        pass
