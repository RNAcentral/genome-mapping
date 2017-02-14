import json

import attr

from attr.validators import optional
from attr.validators import instance_of as is_a

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


@attr.s(frozen=True, slots=True)
class SequenceSummary(object):
    uri = attr.ib(validator=IS_STR)
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


@attr.s(frozen=True, slots=True)
class Stats(object):
    total_gaps = attr.ib(validator=IS_INT)
    gaps = attr.ib(validator=is_a(PairStat))
    length = attr.ib(validator=is_a(PairStat))
    completeness = attr.ib(validator=is_a(PairStat))


@attr.s(frozen=True, slots=True)
class Hit(object):
    urs = attr.ib(validator=IS_STR)
    chromosome = attr.ib(validator=IS_STR)
    start = attr.ib(validator=IS_INT)
    stop = attr.ib(validator=IS_INT)
    fragments = attr.ib(validator=is_a(list))
    is_forward = attr.ib(validator=IS_BOOL)
    input_sequence = attr.ib(validator=is_a(SequenceSummary))
    stats = attr.ib(validator=is_a(Stats))

    @property
    def uri(self):
        return self.input_sequence.uri


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
    match = attr.ib(optional(IS_STR))

    @classmethod
    def build(cls, shift, hit, feature):
        shift_type = '{match_type}_{location_type}'

        if not hit and not feature:
            raise ValueError("Atleast one of hit and feature must be given")

        if not feature:
            return cls(location='novel', match=None, pretty='novel')

        if not hit:
            return cls(location=None, match='missing', pretty='missing')

        if hit.uri == feature.attributes['Name'][0]:
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
                location_type = 'UNKNOWN'

        return cls(
            pretty=shift_type.format(match_type=match_type,
                                     location_type=location_type),
            location=location_type,
            match=match_type,
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
                   stop=match.stop - feature.end)

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
class Feature(object):
    data = attr.ib(validator=optional(IS_TUPLE), hash=False, repr=False)
    pretty = attr.ib(validator=IS_STR, cmp=False)

    @classmethod
    def build(cls, feature):
        pretty = ''
        data = None
        if feature is not None:
            data = list(feature.astuple())
            data[9] = json.loads(data[9])
            data[10] = json.loads(data[10])
            data = tuple(data)
            pretty = str(feature)
        return cls(data=data, pretty=pretty)


@attr.s(frozen=True, slots=True)
class Comparision(object):
    hit = attr.ib(validator=optional(is_a(Hit)))
    feature = attr.ib(validator=optional(is_a(Feature)))
    shift = attr.ib(validator=is_a(Shift))
    type = attr.ib(validator=is_a(ComparisionType))

    @classmethod
    def build(cls, hit, feature):
        if not hit and not feature:
            raise ValueError("At least one of hit and feature must be given")

        shift = Shift.build(hit, feature)
        type = ComparisionType.build(shift, hit, feature)
        feat = Feature.build(feature)
        return cls(hit=hit, feature=feat, shift=shift, type=type)
