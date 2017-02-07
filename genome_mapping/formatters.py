import abc
import sys
import json

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

    def __call__(self, data, stream):
        json.dump([attr.asdict(d) for d in data], stream)


class Gff3(Base):
    name = 'gff'

    def format_feature(self, feature):
        return feature.pretty + '\n'

    def format_hit(self, hit, **extra):
        attributes = {
            'Name': [hit.input_sequence.upi],
            'Header': [hit.input_sequence.header],
        }
        attributes.update(extra or {})

        return str(Feature(
            seqid=hit.chromosome,
            source='map_sequences.py',
            featuretype='hit',
            start=hit.start,
            end=hit.stop,
            score='',
            strand='+' if hit.is_forward else '-',
            frame='.',
            attributes=attributes,
        )) + '\n'

    def __call__(self, data, stream):
        stream.write('##gff-version 3\n')
        for entry in data:
            if isinstance(entry, dat.Hit):
                stream.write(self.format_hit(entry))
            elif isinstance(entry, dat.Feature):
                stream.write(self.format_feature(entry))
            elif isinstance(entry, dat.Comparision):
                if entry.hit:
                    stream.write(self.format_hit(entry.hit,
                                                 type=[entry.type.pretty]))

                if entry.feature.data:
                    stream.write(self.format_feature(entry.feature))
            else:
                raise ValueError('Cannot format all data to gff')
