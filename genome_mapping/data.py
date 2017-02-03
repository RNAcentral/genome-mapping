import attr

RESULT_TYPE = frozenset(['exact', 'additional', 'missing', 'inexact'])

IS_INT = attr.validators.instance_of(int)
IS_FLOAT = attr.validators.instance_of(float)
IS_BOOL = attr.validators.instance_of(bool)
IS_STR = attr.validators.instance_of(basestring)
IS_SET = attr.validators.instance_of(set)


@attr.s(frozen=True)
class Stats(object):
    identical = attr.ib(validator=IS_INT)
    identity = attr.ib(validator=IS_FLOAT)
    gaps = attr.ib(validator=IS_INT)
    query_length = attr.ib(validator=IS_INT)
    hit_length = attr.ib(validator=IS_INT)


@attr.s(frozen=True)
class Hit(object):
    name = attr.ib(validator=IS_STR)
    chromosome = attr.ib(validator=IS_STR)
    start = attr.ib(validator=IS_INT)
    stop = attr.ib(validator=IS_INT)
    is_forward = attr.ib(validator=IS_BOOL)
    input_sequence = attr.ib()
    stats = attr.ib(validator=attr.validators.instance_of(Stats))


@attr.s(frozen=True)
class Shift(object):
    start = attr.ib()
    stop = attr.ib()

    @classmethod
    def cross_chromosome(cls):
        return cls(start=float('-inf'), stop=float('inf'))

    @property
    def total(self):
        return abs(self.start) + abs(self.stop)

    @property
    def type(self):
        if not self.start and not self.stop:
            return 'exact'
        # if self.start > 0 and self.stop < 0:
        #     return 'within'
        return 'inexact'


@attr.s(frozen=True)
class Comparision(object):
    hit = attr.ib(validator=attr.validators.instance_of(Hit))
    feature = attr.ib()
    shift = attr.ib(validator=attr.validators.instance_of(Shift))

    @property
    def type(self):
        if not self.feature:
            return 'additional'
        if not self.hit:
            return 'missing'
        if self.hit and self.feature:
            return self.shift.type()
