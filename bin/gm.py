#!/usr/bin/env python -W ignore

import csv
import pickle
import itertools as it
from pprint import pprint
import collections as coll

import click

from genome_mapping import utils
from genome_mapping import mappers
from genome_mapping import matchers
from genome_mapping import formatters
from genome_mapping.intervals import Tree
from genome_mapping.data import Comparision
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
@click.argument('genome', type=click.Path(readable=True))
@click.argument('targets', type=click.Path(exists=True, readable=True))
@click.argument('save', type=WritablePickleFile())
@click.option('--method', default='blat',
              type=click.Choice(mappers.known()))
@click.option('--organism', default='UNKNOWN')
def find(genome, targets, save, method='blat', organism='UNKNOWN'):
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


@cli.group('hits')
def hits():
    """
    Group of commands dealing with the hits.
    """
    pass


@hits.command('select')
@click.argument('hits', type=ReadablePickleFile())
@click.argument('matcher', type=click.Choice(matchers.known()))
@click.argument('save', type=WritablePickleFile())
@click.option('--define', multiple=True, default={}, type=KeyValue())
def hits_select(hits, matcher, save, define={}):
    """
    Filter out the hits that only pass a match criteria. A match criteria is
    something like 'exact' or sequence identity.
    """
    matcher_class = matchers.fetch(matcher)
    matcher = matcher_class(**(define or {}))
    save(list(matcher.filter_matches(hits)))


@hits.command('compare')
@click.argument('hits', type=ReadablePickleFile())
@click.argument('correct', type=click.Path(exists=True, readable=True))
@click.argument('save', type=WritablePickleFile())
def compare_matches(hits, correct, save):
    """
    Compare some hits to known examples to see how well they overlap.
    """
    tree = Tree(correct)
    save(tree.compare_to_known(hits))


@hits.command('best-within')
@click.argument('hits', type=ReadablePickleFile())
@click.argument('features', type=click.Path(exists=True, readable=True))
@click.argument('max-range', type=int)
@click.argument('save', type=WritablePickleFile())
def best_within(hits, features, max_range, save):
    """
    Find the best hits within some distance of the given features.
    """
    tree = Tree(features)
    save(tree.best_hits_within(hits, max_range))


@hits.command('merge')
@click.argument('hits', type=ReadablePickleFile(), nargs=-1)
@click.argument('save', type=WritablePickleFile())
def merge_hits(hits, save):
    """Merge several collections hits. This should not produce any duplicate
    hits.
    """
    save(sorted(set(it.chain.from_iterable(hits))))


@cli.group('comparisions')
def comparisions():
    """Group of commands dealing with maninpulating comparisions.
    """
    pass


@comparisions.command('select')
@click.argument('comparisions', type=ReadablePickleFile())
@click.argument('filter', type=str, nargs=-1)
@click.argument('save', type=WritablePickleFile())
def comparisions_select(comparisions, filter, save):
    """
    Filter comparisions to only those of the given type(s).
    """

    replace = {'is': '=='}
    processed = ' '.join(replace.get(w, w) for w in filter)
    ast = compile(processed, '<string>', mode='eval')

    def checker(obj):
        fields = utils.properities_of(obj.__class__)
        locals = {(f, getattr(obj, f)) for f in fields}
        each = it.chain.from_iterable(t.split('_') for t in RESULT_TYPE)
        locals.update((t, t) for t in it.chain(RESULT_TYPE, each))
        return eval(ast, globals(), dict(locals))

    save([comp for comp in comparisions if checker(comp)])


@comparisions.command('extract')
@click.argument('comparisions', type=ReadablePickleFile())
@click.argument('property')
@click.argument('save', type=WritablePickleFile())
@click.option('--skip-missing', is_flag=True, default=False)
def comparisions_extract(comparisions, property, save, skip_missing=False):
    """Extract parts of the given comparisions.
    """
    extracted = []
    for comparision in comparisions:
        entry = getattr(comparision, property)
        if skip_missing and not entry:
            continue
        extracted.append(entry)
    save(extracted)


@comparisions.command('summary')
@click.argument('comparisions', type=ReadablePickleFile())
@click.argument('save', type=click.File(mode='wb'))
def summarize_hits(comparisions, save):
    """
    Compute a summary of the number of each type of comparisions.
    """
    summary = dict(coll.Counter(c.type.pretty for c in comparisions))
    writer = csv.DictWriter(save, sorted(summary.keys()))
    writer.writeheader()
    writer.writerow(summary)


@cli.command('as')
@click.argument('format', type=click.Choice(formatters.known()))
@click.argument('data', type=ReadablePickleFile())
@click.argument('save', type=click.File(mode='wb'))
def format(format, data, save):
    """
    Format the data into some format.

    Paramaters
    ----------
    format : str
        The name of the format to use.
    data : path
        The path to the data to read, '-' means stdin.
    save : path
         The path of where to save data, '-' means stdout.
    """
    formatters.format(data, format, save)


@cli.command('pp')
@click.argument('data', type=ReadablePickleFile())
def display(data):
    """
    Pretty print data to stdout.
    """
    pprint(data)


if __name__ == '__main__':
    cli()
