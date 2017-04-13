import abc
import csv
import sys
import json
import itertools as it

import attr
from gffutils import Feature

from genome_mapping import data as dat
from genome_mapping import utils as ut


def known():
    return ut.names_of_children(sys.modules[__name__], Base)


def fetch(name):
    return ut.get_child(sys.modules[__name__], Base, name)


def format(data, name, stream):
    klass = fetch(name)
    obj = klass()
    return obj(data, stream)


class Base(object):
    __metaclass__ = abc.ABCMeta

    @abc.abstractproperty
    def name(self):
        pass


class Json(Base):
    name = 'json'

    def format(self, data):
        if attr.has(data.__class__):
            return attr.asdict(data)
        if isinstance(data, (tuple, list)):
            return [self.format(d) for d in data]
        if isinstance(data, dict):
            return {k: self.format(v) for k, v in data.iteritems()}

        raise ValueError("Cannot handle this data type")

    def __call__(self, data, stream):
        return json.dump(self.format(data), stream)


class Gff3(Base):
    name = 'gff3'

    def format_feature(self, feature):
        return feature.pretty

    def format_hit(self, hit, hit_type=None, **extra):
        attributes = {
            'Name': hit.input_sequence.urs,
            'Header': hit.input_sequence.header,
            'HitSize': hit.stats.length.hit,
            'QuerySize': hit.stats.length.query,
        }
        attributes.update(extra or {})

        normalized = {}
        for key, value in attributes.items():
            if not isinstance(value, list):
                normalized[key] = [str(value)]
            else:
                normalized[key] = [str(e) for e in value]

        feature_type = 'hit'
        if hit_type:
            feature_type = '%s-hit' % hit_type

        return str(Feature(
            seqid=hit.chromosome,
            source='map_sequences.py',
            featuretype=feature_type,
            start=hit.start + 1,
            end=hit.stop,
            score='',
            strand='+' if hit.is_forward else '-',
            frame='.',
            attributes=normalized,
        ))

    def format(self, entry):
        if isinstance(entry, dat.Hit):
            yield self.format_hit(entry)
        elif isinstance(entry, dat.FeatureData):
            yield self.format_feature(entry)
        elif isinstance(entry, dat.Comparision):
            if entry.hit:
                yield self.format_hit(entry.hit,
                                      hit_type=entry.type.match,
                                      type=entry.type.pretty)

                if entry.feature:
                    yield self.format_feature(entry.feature)
            raise ValueError('Cannot format all data to gff')

    def __call__(self, data, stream):
        assert isinstance(data, (list, tuple))
        stream.write('##gff-version 3\n')
        for entry in data:
            for row in self.format(data):
                stream.write(row + '\n')


class Insertable(Base):
    name = 'insertable'

    def format(self, data):
        if isinstance(data, (tuple, list)):
            return [self.format(d) for d in data]
        elif isinstance(data, dat.Hit):
            urs, taxid = data.urs.split('_')
            return {
                'urs': urs,
                'taxid': int(taxid),
                'exons': self.format(data.fragments)
            }
        elif isinstance(data, dat.Fragment):
            return {
                'primary_start': data.start + 1,
                'primary_end': data.stop,
                'name': data.chromosome,
                'strand': 1 if data.is_forward else -1
            }
        raise ValueError("Cannot handle given data")

    def __call__(self, data, stream):
        return json.dump(self.format(data), stream)
