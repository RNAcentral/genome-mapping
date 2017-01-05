import attr

is_int = attr.validators.instance_of(int)
is_bool = attr.validators.instance_of(bool)
is_str = attr.validators.instance_of(str)
is_set = attr.validators.instance_of(set)


@attr.s()
class SummaryCounts(object):
    valid_matches = attr.ib(default=0)
    invalid_matches = attr.ib(default=0)
    extra_matches = attr.ib(default=0)
    missing_matches = attr.ib(default=0)
    incomplete_matches = attr.ib(default=0)


@attr.s(frozen=True)
class Mapping(object):
    chromosome = attr.ib(validator=is_int)
    start = attr.ib(validator=is_int)
    stop = attr.ib(validator=is_int)
    is_forward = attr.ib(validator=is_bool)

    def __hash__(self):
        return hash((self.chromosome, self.start, self.stop, self.is_forward))


@attr.s()
class SequenceResults(object):
    md5 = attr.ib(validator=is_str)
    mappings = attr.ib(validator=is_set, default=attr.Factory(set()))

    def add(self, mapping):
        self.mappings.add(mapping)

    def __hash__(self):
        return hash(self.name)


@attr.s(frozen=True)
class HitStats(object):
    identical = attr.ib(validator=is_int)
    gaps = attr.ib(validator=is_int)
    query_length = attr.ib(validator=is_int)
    hit_length = attr.ib(validator=is_int)


@attr.s(frozen=True)
class MappingHit(object):
    chromosome = attr.ib(validator=is_int)
    start = attr.ib(validator=is_int)
    stop = attr.ib(validator=is_int)
    is_forward = attr.ib(validator=is_bool)
    stats = attr.ib(validator=attr.validators.instance_of(HitStats))

    def as_mapping(self):
        return Mapping(chromosome=self.chromosome,
                       start=self.start,
                       stop=self.stop,
                       is_forward=self.is_forward
                       )
