#!/usr/bin/env python

import csv
import json
import pickle
from pprint import pprint
import collections as coll

import attr
import click

from genome_mapping import mappers
from genome_mapping import matchers
from genome_mapping.intervals import Tree
from genome_mapping.data import RESULT_TYPE


class ReadablePickleFile(click.File):
    name = 'readable-pickle'

    def __init__(self):
        super(ReadablePickleFile, self).__init__(mode='rb')

    def convert(self, value, param, ctx):
        fileobj = super(ReadablePickleFile, self).convert(value, param, ctx)
        return pickle.load(fileobj)


class WritablePickleFile(click.File):
    name = 'writeable-pickle'

    def __init__(self):
        super(WritablePickleFile, self).__init__(mode='wb')

    def convert(self, value, param, ctx):
        fileobj = super(WritablePickleFile, self).convert(value, param, ctx)
        return lambda d: pickle.dump(d, fileobj)


class KeyValue(click.ParamType):
    name = 'key-value'

    def convert(self, value, param, ctx):
        if value is None or value == {}:
            return value
        key, value = value.split('=')
        return {key: value}


@click.group()
def cli():
    pass


@cli.command('find')
@click.argument('genome', type=click.Path(exists=True, readable=True))
@click.argument('targets', type=click.Path(exists=True, readable=True))
@click.argument('save', type=WritablePickleFile())
@click.option('--method', default='blat',
              type=click.Choice(mappers.known()))
def hits(genome, targets, save, method='blat'):
    """
    Search the genome for the given targets using the specified program.

    Parameters
    ----------

    genome :
        The path to the FASTA file of the genome to search in.

    targets :
        The path to the FASTA file of targets to search with.

    save :
        The path to the file to save pickled objects to.
    """
    mapper_class = mappers.fetch(method)
    mapper = mapper_class()
    save(list(mapper(genome, targets)))


@cli.command('select-matches')
@click.argument('hits', type=ReadablePickleFile())
@click.argument('save', type=WritablePickleFile())
@click.option('--matcher', default='exact',
              type=click.Choice(matchers.known()))
@click.option('--define', multiple=True, default={}, type=KeyValue())
def filter_hits(hits, save, matcher='exact', define={}):
    """
    Filter out the hits that only pass a match criteria. A match criteria is
    something like 'exact' or sequence identity.
    """
    matcher_class = matchers.fetch(matcher)
    matcher = matcher_class(**(define or {}))
    save(list(matcher.filter_matches(hits)))


@cli.command('compare')
@click.argument('hits', type=ReadablePickleFile())
@click.argument('correct', type=click.Path(exists=True, readable=True))
@click.argument('save', type=WritablePickleFile())
def compare_matches(hits, correct, save):
    """Compare some hits to known examples to see how well they overlap.
    """
    tree = Tree(correct)
    save(tree.compare_to_known(hits))


@cli.command('select-type')
@click.argument('comparisions', type=ReadablePickleFile())
@click.argument('save', type=WritablePickleFile())
@click.option('--types', multiple=True, type=click.Choice(RESULT_TYPE))
def extract_hit_class(comparisions, save, types):
    """
    Filter comparisions to only those of the given type(s).
    """
    save([comp for comp in comparisions if comp.type in set(types)])


@cli.group('extract')
def extract():
    pass


@extract.command('hit')
@click.argument('comparisions', type=ReadablePickleFile())
@click.argument('save', type=WritablePickleFile())
def extract_hits(comparisions, save):
    """
    Extract only the hits from comparisions.
    """
    save([c.hit for c in comparisions])


@extract.command('feature')
@click.argument('comparisions', type=ReadablePickleFile())
@click.argument('save', type=WritablePickleFile())
def extract_features(comparisions, save):
    """
    Extract the features from comparisions.
    """
    save([c.feature for c in comparisions])


@cli.command('summary')
@click.argument('comparisions', type=ReadablePickleFile())
@click.argument('save', type=click.File(mode='wb'))
def summarize_hits(comparisions, save):
    """
    Compute a summary of the number of each type of comparisions.
    """
    summary = dict(coll.Counter(c.type for c in comparisions))
    writer = csv.DictWriter(save, sorted(summary.keys()))
    writer.writeheader()
    writer.writerow(summary)


@cli.group('as')
def format():
    """
    Format the data into some format.
    """
    pass


@format.command('json')
@click.argument('data', type=ReadablePickleFile())
@click.argument('save', type=click.File(mode='wb'))
def as_json(data, save):
    """
    Convert any data to JSON.
    """
    data = [attr.asdict(d) for d in data]
    json.dump(data, save)


@cli.command('pp')
@click.argument('data', type=ReadablePickleFile())
def display(data):
    """
    Pretty print data to stdout.
    """
    pprint(data)


if __name__ == '__main__':
    cli()
