#!/usr/bin/env python -W ignore

import os
import csv
import sys
import cPickle
import json
import itertools as it
from pprint import pprint
import collections as coll

import click

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from genome_mapping import utils
from genome_mapping import mappers
from genome_mapping import matchers
from genome_mapping import formatters
from genome_mapping.intervals import Tree
from genome_mapping.data import RESULT_TYPE


class ReadableDataFile(click.File):
    name = 'readable-pickle'

    def __init__(self):
        super(ReadableDataFile, self).__init__(mode='rb')

    def convert(self, value, param, ctx):
        fileobj = super(ReadableDataFile, self).convert(value, param, ctx)
        try:
            while True:
                yield cPickle.load(fileobj)
        except EOFError:
            raise StopIteration


class WritableDataFile(click.File):
    name = 'writeable-pickle'

    def __init__(self):
        super(WritableDataFile, self).__init__(mode='wb')

    def convert(self, value, param, ctx):
        fileobj = super(WritableDataFile, self).convert(value, param, ctx)
        return lambda d: cPickle.dump(d, fileobj)


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
@click.argument('save', type=WritableDataFile())
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
    save(mapper(genome, targets))


@cli.group('hits')
def hits():
    """
    Group of commands dealing with the hits.
    """
    pass


@hits.command('from-format')
@click.argument('data', type=click.Path(exists=True, readable=True))
@click.argument('targets', type=click.Path(exists=True, readable=True))
@click.argument('save', type=WritableDataFile())
@click.option('--format',
              type=click.Choice(mappers.known_formats().keys() + [None]),
              default=None)
def format_to_hits(data, targets, save, format=None):
    if not format:
        _, ext = os.path.splitext(data)
        format = ext[1:]
        if format not in mappers.known_formats():
            raise ValueError("Unknown inferred format %s" % format)

    for hit in mappers.from_format(data, targets, format):
        save(hit)


@hits.command('select')
@click.argument('hits', type=ReadableDataFile())
@click.argument('matcher', type=click.Choice(matchers.known()))
@click.argument('save', type=WritableDataFile())
@click.option('--define', multiple=True, default={}, type=KeyValue())
def hits_select(hits, matcher, save, define={}):
    """
    Filter out the hits that only pass a match criteria. A match criteria is
    something like 'exact' or sequence identity.
    """
    definitions = {}
    for definition in define:
        definitions.update(definition)

    matcher_class = matchers.fetch(matcher)
    matcher = matcher_class(**definitions)
    for filtered in matcher.filter_matches(hits):
        save(filtered)


@hits.command('select-using-spec')
@click.argument('hits', type=ReadableDataFile())
@click.argument('spec-file', type=click.File('rb'))
@click.argument('save', type=WritableDataFile())
def hits_select_spec(hits, spec_file, save):
    """
    Select hits using the specifications in the given file. The file be a json
    object with a matcher entry that is the name of the matcher to use. It may
    also contain a JSON object of definitions to build the matcher with.
    """
    spec = json.load(spec_file)
    matcher = spec['matcher']
    if matcher not in matchers.known():
        raise ValueError("Unknown Matcher")
    definitions = spec.get('definitions', {})
    matcher_class = matchers.fetch(matcher)
    matcher = matcher_class(**definitions)
    for filtered in matcher.filter_matches(hits):
        save(filtered)


@hits.command('compare')
@click.argument('hits', type=ReadableDataFile())
@click.argument('correct', type=click.Path(exists=True, readable=True))
@click.argument('save', type=WritableDataFile())
def compare_matches(hits, correct, save):
    """
    Compare some hits to known examples to see how well they overlap.
    """
    tree = Tree(correct)
    for compared in tree.compare_to_known(hits):
        save(compared)


@hits.command('best-within')
@click.argument('hits', type=ReadableDataFile())
@click.argument('features', type=click.Path(exists=True, readable=True))
@click.argument('max-range', type=int)
@click.argument('save', type=WritableDataFile())
def best_within(hits, features, max_range, save):
    """
    Find the best hits within some distance of the given features.
    """
    tree = Tree(features)
    for best in tree.best_hits_within(hits, max_range):
        save(best)


@cli.group('comparisons')
def comparisons():
    """Group of commands dealing with maninpulating comparisons.
    """
    pass


@comparisons.command('select')
@click.argument('comparisons', type=ReadableDataFile())
@click.argument('filter', type=str, nargs=-1)
@click.argument('save', type=WritableDataFile())
def comparisons_select(comparisons, filter, save):
    """
    Filter comparisons to only those of the given type(s).
    """

    replace = {'is': '=='}
    processed = ' '.join(replace.get(w, w) for w in filter)
    ast = compile(processed, '<string>', mode='eval')

    def checker(obj):
        def hashable(name):
            return isinstance(getattr(obj, name), coll.Hashable)

        fields = utils.properities_of(obj.__class__)
        locals = {(f, getattr(obj, f)) for f in fields if hashable(f)}
        each = it.chain.from_iterable(t.split('_') for t in RESULT_TYPE)
        locals.update((t, t) for t in it.chain(RESULT_TYPE, each))
        return eval(ast, {}, dict(locals))

    save(comp for comp in comparisons if checker(comp))


@comparisons.command('extract')
@click.argument('comparisons', type=ReadableDataFile())
@click.argument('property')
@click.argument('save', type=WritableDataFile())
@click.option('--skip-missing', is_flag=True, default=False)
def comparisons_extract(comparisons, property, save, skip_missing=False):
    """Extract parts of the given comparisons.
    """
    for comparision in comparisons:
        entry = getattr(comparision, property)
        if skip_missing and not entry:
            continue
        save(entry)


@comparisons.command('summary')
@click.argument('comparisons', type=ReadableDataFile())
@click.argument('save', type=click.File(mode='wb'))
def comparisons_summarize(comparisons, save):
    """
    Compute a summary of the number of each type of comparisons.
    """
    pretty_counts = dict(coll.Counter(c.type.pretty for c in comparisons))
    match_counts = dict(coll.Counter(c.type.match for c in comparisons))
    counts = dict(pretty_counts)
    counts.update(match_counts)
    keys = sorted(RESULT_TYPE)
    writer = csv.DictWriter(save, ['type', 'count'])
    writer.writeheader()
    writer.writerows([{'type': k, 'count': counts[k]} for k in keys])
    writer.writerow({'type': 'total', 'count': len(comparisons)})



@comparisons.group('group')
def comparisons_group():
    pass


@comparisons_group.command('by-hit-urs')
@click.argument('comparisons', type=ReadableDataFile())
@click.argument('save', type=WritableDataFile())
def comparisons_group_hit_urs(comparisons, save):
    def key(comparision):
        if comparision.hit:
            return comparision.hit.urs
        return None
    ordered = sorted(comparisons, key=key)
    groups = it.groupby(ordered, key)
    save({urs: list(comps) for urs, comps in groups})


@comparisons_group.command('summarize-types')
@click.argument('grouped', type=ReadableDataFile())
@click.argument('save', type=click.File(mode='wb'))
def comparisons_group_summary_type(grouped, save):
    header = ['urs'] + list(RESULT_TYPE)
    writer = csv.DictWriter(save, header)
    writer.writeheader()
    for urs, hits in grouped.iteritems():
        entry = coll.defaultdict(int)
        entry['urs'] = urs
        entry.update(coll.Counter(h.type.pretty for h in hits))
        writer.writerow(entry)


@cli.command('as')
@click.argument('format', type=click.Choice(formatters.known()))
@click.argument('data', type=ReadableDataFile())
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
@click.argument('data', type=ReadableDataFile())
def display(data):
    """
    Pretty print data to stdout.
    """
    pprint(data)


if __name__ == '__main__':
    cli()
