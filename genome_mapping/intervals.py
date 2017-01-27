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
        return cls(start=match.start - feature.begin,
                   stop=match.start - feature.end)

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
        return cls(feature=feature,
                   location=match,
                   shift=Shift.from_match(match, feature))

    def is_exact(self):
        return self.shift.is_exact()


class Tree(object):
    def __init__(self, filename):
        self.db = gff.create_db(filename, ':memory:')
        intervals = (Interval(f.start, f.end, f) for f in self.features())
        self.tree = IntervalTree(intervals)
        self.locations = coll.defaultdict(list)
        for feature in self.features():
            for name in feature.attributes['Name']:
                self.locations[name].append(feature)

    def features(self):
        seen = set()
        for feature in self.db.all_features():
            name = feature.attributes['Name'][0]
            key = (name, feature.start, feature.end)
            if key in seen:
                continue
            seen.add(key)
            yield feature

    def search(self, start, stop):
        matches = set()
        for interval in self.tree.search(start, stop):
            matches.add(Overlap.from_interval(start, stop, interval))
        return matches

    def find(self, upi):
        return self.locations[upi]

    def find_overlaps(self, matches):
        overlaps = []
        for match in matches:
            seen = set()
            for overlap in self.search(match.start, match.stop):
                seen.add(overlap.feature.id)
                overlaps.append(overlap)

            for correct in self.find(match.name):
                if correct.id in seen:
                    continue
                overlaps.append(Overlap.from_match(match, correct))

        return overlaps
