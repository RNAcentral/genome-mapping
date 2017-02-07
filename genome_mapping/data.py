import attr

RESULT_TYPE = frozenset([
    'correct/exact',
    'correct/within',
    'correct/enclose',
    'correct/5p-overlap',
    'correct/3p-overlap',
    'incorrect/exact',
    'incorrect/within',
    'incorrect/enclose',
    'incorrect/5p-overlap',
    'incorrect/3p-overlap',
    'novel',
    'missing',
])

IS_INT = attr.validators.instance_of(int)
IS_FLOAT = attr.validators.instance_of(float)
IS_BOOL = attr.validators.instance_of(bool)
IS_STR = attr.validators.instance_of(basestring)
IS_SET = attr.validators.instance_of(set)


@attr.s(frozen=True, slots=True)
class SequenceSummary(object):
    upi = attr.ib(validator=IS_STR)
    id = attr.ib(validator=IS_STR)
    header = attr.ib(validator=IS_STR)


@attr.s(frozen=True, slots=True)
class Stats(object):
    identical = attr.ib(validator=IS_INT)
    identity = attr.ib(validator=IS_FLOAT)
    gaps = attr.ib(validator=IS_INT)
    query_length = attr.ib(validator=IS_INT)
    hit_length = attr.ib(validator=IS_INT)


@attr.s(frozen=True, slots=True)
class Hit(object):
    name = attr.ib(validator=IS_STR)
    chromosome = attr.ib(validator=IS_STR)
    start = attr.ib(validator=IS_INT)
    stop = attr.ib(validator=IS_INT)
    is_forward = attr.ib(validator=IS_BOOL)
    input_sequence = attr.ib()
    stats = attr.ib(validator=attr.validators.instance_of(Stats))


@attr.s(frozen=True, slots=True)
class LocationMatchType(object):
    location = attr.ib()
    match = attr.ib()
    pretty = attr.ib(validator=IS_STR)

    @classmethod
    def build(cls, shift, hit, feature):
        shift_type = '{match_type}/{location_type}'

        if not hit and not feature:
            return cls(location=None, match=None, pretty='IMPOSSIBLE')

        if not feature:
            return cls(location='novel', match=None, pretty='novel')

        if not hit:
            return cls(location=None, match='missing', pretty='missing')

        if hit.input_sequence.upi == feature.attributes['Name'][0]:
            match_type = 'correct'
        else:
            match_type = 'incorrect'

        if shift.is_exact():
            location_type = 'exact'
        else:
            if shift.start >= 0 and shift.stop <= 0:
                location_type = 'within'
            if shift.start < 0 and shift.stop > 0:
                location_type = 'enclose'
            if shift.start > 0 and shift.stop >= 0:
                location_type = "3p-shift"
            if shift.start < 0 and shift.stop <= 0:
                location_type = "5p-shift"
            else:
                # print(hit)
                # print(feature)
                # print(shift)
                location_type = 'UNKNOWN'

        return cls(
            location=location_type,
            match=match_type,
            pretty=shift_type.format(match_type=match_type,
                                     location_type=location_type),
        )


@attr.s(frozen=True, slots=True)
class Shift(object):
    start = attr.ib()
    stop = attr.ib()

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

    def is_exact(self):
        return (not self.start and not self.stop) or \
            (self.start == -1 and not self.stop)


@attr.s(frozen=True, slots=True)
class Comparision(object):
    hit = attr.ib()
    feature = attr.ib()
    shift = attr.ib(validator=attr.validators.instance_of(Shift))
    type = attr.ib()

    @classmethod
    def build(cls, hit, feature):
        feat = feature
        if feat is not None:
            feat = str(feat)
        shift = Shift.build(hit, feature)
        type = LocationMatchType.build(shift, hit, feature)
        return cls(hit=hit, feature=feat, shift=shift, type=type)
