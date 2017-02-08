#!/usr/bin/env python

import csv
import hashlib

import click
from Bio import SeqIO


class ReadableBioFile(click.File):
    def __init__(self, format):
        super(ReadableBioFile, self).__init__(mode='rb')
        self._format = format

    def convert(self, value, param, ctx):
        fileobj = super(ReadableBioFile, self).convert(value, param, ctx)
        return SeqIO.parse(fileobj, self._format)


class WriteableBioFile(click.File):
    def __init__(self, format):
        super(WriteableBioFile, self).__init__(mode='wb')
        self._format = format

    def convert(self, value, param, ctx):
        fileobj = super(WriteableBioFile, self).convert(value, param, ctx)
        return lambda s: SeqIO.write(s, fileobj, self._format)


@click.group()
def main():
    pass


@main.group()
def md5():
    pass


@md5.command('display')
@click.argument('fasta', type=ReadableBioFile('fasta'))
@click.argument('save', type=click.File(mode='wb'))
def md5_display(fasta, save):
    writer = csv.writer(save, delimiter='\t')
    for record in fasta:
        seq = str(record.seq).strip()
        md5 = hashlib.md5(seq).hexdigest()
        writer.writerow([record.id, md5])


@md5.command('validate')
@click.argument('filename', type=click.File(mode='rb'))
@click.argument('correct', type=click.File(mode='rb'))
@click.argument('save', type=click.File(mode='wb'))
def md5_validate(filename, correct, save):
    correct = {r[0]: r[1] for r in csv.reader(correct, delimiter='\t')}
    writer = csv.writer(save, delimiter='\t')
    for row in csv.reader(filename, delimiter='\t'):
        status = 'OK' if row[1] == correct[row[0]] else 'NOT OK'
        writer.writerow((row[0], status))


@main.command('dna-to-rna')
@click.argument('fasta', type=ReadableBioFile('fasta'))
@click.argument('save', type=WriteableBioFile('fasta'))
def dna_to_rna(fasta, save):
    for record in fasta:
        record.seq = record.seq.transcribe()
        save(record)


@main.command('merge-by-id')
@click.argument('fasta', type=ReadableBioFile('fasta'))
@click.argument('save', type=WriteableBioFile('fasta'))
@click.option('--validate-hash', is_flag=True, default=False)
def merge_by_id(fasta, save, validate_hash):
    unique = {}
    for record in fasta:
        if record.id not in unique:
            unique[record.id] = record
            continue

        rest = record.description.replace(record.id, '')
        unique[record.id].description += ',%s' % rest

    for record in unique.itervalues():
        save(record)


if __name__ == '__main__':
    main()
