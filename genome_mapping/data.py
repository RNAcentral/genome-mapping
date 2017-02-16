import attr

import gffutils as gff
import intervaltree as itt

from attr.validators import optional
from attr.validators import instance_of as is_a

UNKNOWN = 'UNKNOWN'

RESULT_TYPE = frozenset([
    'correct_exact',
    'correct_within',
    'correct_enclose',
    'correct_5p_shift',
    'correct_3p_shift',
    'correct_5p_disjoint',
    'correct_3p_disjoint',
    'incorrect_exact',
    'incorrect_within',
    'incorrect_enclose',
    'incorrect_5p_shift',
    'incorrect_3p_shift',
    'incorrect_5p_disjoint',
    'incorrect_3p_disjoint',
    'correct_UNKNOWN',
    'incorrect_UNKNOWN',
    'novel',
    'missing',
])


IS_INT = is_a(int)
IS_FLOAT = is_a(float)
IS_BOOL = is_a(bool)
IS_STR = is_a(basestring)
IS_SET = is_a(set)
IS_NUM = is_a((int, float))
IS_TUPLE = is_a(tuple)
IS_DICT = is_a(dict)
IS_LIST = is_a(list)


def urs_of(data):
    if hasattr(data, 'urs'):
        return data.urs
    if isinstance(data, itt.Interval):
        return urs_of(data.data)
    if isinstance(data, gff.Feature):
        return data.attributes['Name'][0]
    raise ValueError("No way to get urs of: %s" % data)


@attr.s(frozen=True, slots=True)
class SequenceSummary(object):
    urs = attr.ib(validator=IS_STR)
    id = attr.ib(validator=IS_STR)
    header = attr.ib(validator=IS_STR)
    length = attr.ib(validator=IS_INT)


@attr.s(frozen=True, slots=True)
class PairStat(object):
    query = attr.ib(validator=IS_NUM)
    hit = attr.ib(validator=IS_NUM)

    @property
    def total(self):
        return self.query + self.hit

    @classmethod
    def from_summation(cls, pairs):
        return cls(
            query=sum(p.query for p in pairs),
            hit=sum(p.hit for p in pairs)
        )


@attr.s(frozen=True, slots=True)
class Stats(object):
    total_gaps = attr.ib(validator=IS_INT)
    identical = attr.ib(validator=IS_INT)
    gaps = attr.ib(validator=is_a(PairStat))
    length = attr.ib(validator=is_a(PairStat))
    completeness = attr.ib(validator=is_a(PairStat))

    @classmethod
    def from_summation(cls, stats):
        return cls(
            total_gaps=sum(s.total_gaps for s in stats),
            identical=sum(s.identical for s in stats),
            gaps=PairStat.from_summation(s.gaps for s in stats),
            length=PairStat.from_summation(s.length for s in stats),
            completeness=PairStat.from_summation(s.completeness for s in stats),
        )


@attr.s(frozen=True, slots=True)
class Hit(object):
    urs = attr.ib(validator=IS_STR)
    chromosome = attr.ib(validator=IS_STR)
    start = attr.ib(validator=IS_INT)
    stop = attr.ib(validator=IS_INT)
    fragments = attr.ib(validator=IS_LIST, hash=False)
    is_forward = attr.ib(validator=IS_BOOL)
    input_sequence = attr.ib(validator=is_a(SequenceSummary))
    stats = attr.ib(validator=is_a(Stats))

    @property
    def sequence_type(self):
        if len(self.fragments) == 1:
            return 'unspliced'
        if len(self.fragments) > 1:
            return 'spliced'
        return UNKNOWN


@attr.s(frozen=True, slots=True)
class Fragment(object):
    name = attr.ib(validator=IS_STR)
    chromosome = attr.ib(validator=IS_STR)
    start = attr.ib(validator=IS_INT)
    stop = attr.ib(validator=IS_INT)
    is_forward = attr.ib(validator=IS_BOOL)
    stats = attr.ib(validator=is_a(Stats))


@attr.s(frozen=True, slots=True)
class ComparisionType(object):
    pretty = attr.ib(validator=IS_STR)
    location = attr.ib(validator=optional(IS_STR))
    match = attr.ib(validator=optional(IS_STR))
    feature_type = attr.ib(validator=optional(IS_STR))
    hit_type = attr.ib(validator=optional(IS_STR))

    @classmethod
    def build(cls, shift, hit, feature):
        shift_type = '{match_type}_{location_type}_{feature_type}_{hit_type}'

        if not hit and not feature:
            raise ValueError("Atleast one of hit and feature must be given")

        if not feature:
            return cls(
                location='novel',
                match=None,
                pretty='novel',
                feature_type=None,
                hit_type=hit.sequence_type
            )

        if not hit:
            return cls(
                location=None,
                match='missing',
                pretty='missing',
                feature_type=feature.sequence_type,
                hit_type=None,
            )

        if hit.urs == feature.urs:
            match_type = 'correct'
        else:
            match_type = 'incorrect'

        if shift.is_exact():
            location_type = 'exact'
        else:
            if hit.stop < feature.start:
                location_type = '5p_disjoint'
            elif hit.start > feature.stop:
                location_type = '3p_disjoint'
            elif shift.start >= 0 and shift.stop <= 0:
                location_type = 'within'
            elif shift.start < 0 and shift.stop > 0:
                location_type = 'enclose'
            elif shift.start > 0 and shift.stop >= 0:
                location_type = "3p_shift"
            elif shift.start < 0 and shift.stop <= 0:
                location_type = "5p_shift"
            else:
                location_type = UNKNOWN

        return cls(
            pretty=shift_type.format(match_type=match_type,
                                     location_type=location_type,
                                     feature_type=feature.sequence_type,
                                     hit_type=hit.sequence_type,
                                     ),
            location=location_type,
            match=match_type,
            feature_type=feature.sequence_type,
            hit_type=hit.sequence_type,
        )


@attr.s(frozen=True, slots=True)
class Shift(object):
    start = attr.ib(validator=IS_NUM)
    stop = attr.ib(validator=IS_NUM)

    @classmethod
    def cross_chromosome(cls):
        return cls(start=float('-inf'), stop=float('inf'))

    @classmethod
    def build(cls, match, feature):
        if not match:
            return cls.cross_chromosome()
        if not feature:
            return cls.cross_chromosome()
        return cls(start=match.start - feature.start,
                   stop=match.stop - feature.stop)

    @property
    def total(self):
        return abs(self.start) + abs(self.stop)

    def same_chromosome(self):
        return self.start.start != float('-inf') and \
            self.start.stop != float('inf')

    def is_exact(self):
        return (not self.start and not self.stop) or \
            (self.start == -1 and not self.stop)


@attr.s(frozen=True, slots=True)
class FeatureFragment(object):
    seqid = attr.ib(validator=IS_STR)
    source = attr.ib(validator=IS_STR)
    feature_type = attr.ib(validator=IS_STR)
    start = attr.ib(validator=IS_INT)
    stop = attr.ib(validator=IS_INT)
    score = attr.ib(validator=IS_STR)
    strand = attr.ib(validator=IS_STR)
    frame = attr.ib(validator=IS_STR)
    attributes = attr.ib(validator=IS_DICT, hash=False)
    extra = attr.ib(validator=IS_LIST, hash=False)

    @classmethod
    def build(cls, feature):
        return cls(
            seqid=feature.seqid,
            source=feature.source,
            feature_type=feature.featuretype,
            start=feature.start,
            stop=feature.end,
            score=feature.score,
            strand=feature.strand,
            frame=feature.frame,
            attributes=dict(feature.attributes),
            extra=feature.extra,
        )

    @property
    def urs(self):
        return self.attributes['Name'][0]

    def as_gff(self):
        return gff.Feature(
            seqid=self.seqid,
            source=self.source,
            featuretype=self.feature_type,
            start=self.start,
            end=self.stop,
            score=self.score,
            strand=self.strand,
            frame=self.frame,
            attributes=self.attributes,
            extra=self.extra,
        )


@attr.s(frozen=True, slots=True)
class FeatureData(object):
    chromosome = attr.ib(validator=IS_STR)
    source = attr.ib(validator=IS_STR)
    start = attr.ib(validator=IS_INT)
    stop = attr.ib(validator=IS_INT)
    strand = attr.ib(validator=IS_STR)
    fragments = attr.ib(validator=IS_LIST, hash=False)

    @classmethod
    def build(cls, subfeatures):
        def value_of(name):
            possible = {getattr(f, name) for f in subfeatures}
            if len(possible) > 1:
                raise ValueError("All FeatureFragment's in a FeatureData must"
                                 " have singular %s, found %s" %
                                 (name, possible))
            return possible.pop()

        assert subfeatures
        fragments = [FeatureFragment.build(f) for f in subfeatures]
        return cls(
            chromosome=value_of('seqid'),
            source=value_of('source'),
            start=min(f.start for f in subfeatures),
            stop=max(f.stop for f in subfeatures),
            strand=value_of('strand'),
            fragments=fragments
        )

    @property
    def urs(self):
        return self.fragments[0].urs

    def as_gff(self):
        return [f.as_gff() for f in self.fragments]

    @property
    def pretty(self):
        return '\n'.join(str(f) for f in self.as_gff())

    @property
    def sequence_type(self):
        if len(self.fragments) == 1:
            return 'unspliced'
        if len(self.fragments) > 1:
            return 'spliced'
        return UNKNOWN


@attr.s(frozen=True, slots=True)
class Comparision(object):
    hit = attr.ib(validator=optional(is_a(Hit)))
    feature = attr.ib(validator=optional(is_a(FeatureData)))
    shift = attr.ib(validator=is_a(Shift))
    type = attr.ib(validator=is_a(ComparisionType))

    @classmethod
    def build(cls, hit, feature):
        if not hit and not feature:
            raise ValueError("At least one of hit and feature must be given")

        shift = Shift.build(hit, feature)
        type = ComparisionType.build(shift, hit, feature)
        return cls(hit=hit, feature=feature, shift=shift, type=type)
