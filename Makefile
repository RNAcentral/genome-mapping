gm=bin/gm.py
all_md5=data/md5.tsv

all : pombe

###############################################################################
### Mapping pombe
###############################################################################

pombe : data/pombe/inferred.json data/pombe/summary.csv
	@echo "Found $(jq 'length' data/pombe/inferred.json) coordinates of $(grep -c '^>' data/pombe/unknown.fasta)"

data/pombe/genome.fasta : bin/fetch
	$^ "ftp://ftp.ensemblgenomes.org/pub/current/fungi/fasta/schizosaccharomyces_pombe/dna/Schizosaccharomyces_pombe.ASM294v2.dna.toplevel.fa.gz" > $@

data/pombe/known.gff3 : bin/fetch
	$^ "ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/6.0/genome_coordinates/Schizosaccharomyces_pombe.ASM294v2.gff3.gz" > $@

data/pombe/raw.bed : bin/fetch-bed
	$^ "ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/6.0/genome_coordinates/Schizosaccharomyces_pombe.ASM294v2.bed.gz" > $@

pombe : data/pombe/%-selected.pickle : $(gm) data/pombe/genome.fasta data/pombe/%-hits.pickle
	$(gm) hits select --define min=99 data/pombe/%-hits.pickle best-match  $@

###############################################################################
### Mapping C. elegans
###############################################################################

# celegans: data/celegans/inferred.json data/celegans/summary.csv
# 	@echo "Found $(jq 'length' data/celegans/inferred.json) coordinates of $(grep -c '^>' data/celegans/unknown.fasta)"


###############################################################################
## Here are all the generic definitions that are constant for all mappings
###############################################################################

$(all_md5) : bin/fetch
	$^ "ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral//releases/6.0/md5/md5.tsv.gz" > $@

data/%/known-md5.txt : bin/known-md5s $(all_md5) data/celegans/bad-ids
	$^ > $@

data/%/targets.bed : bin/targets-bed data/%/raw.bed data/%/known.gff3
	$^ > $@

data/%/targets.fasta : bin/targets data/celegans/genome.fasta data/celegans/targets.bed data/celegans/known-md5.txt
	$^ > $@

data/%/targets-hits.pickle : $(gm) data/%/genome.fasta data/%/targets.fasta
	$(gm) find data/celegans/genome.fasta data/celegans/targets.fasta $@

data/%/unknown-hits.pickle : $(gm) data/%/genome.fasta data/%/unknown.fasta
	$(gm) find data/celegans/genome.fasta data/celegans/unknown.fasta $@

data/%s/targets-compared.pickle : $(gm) data/%/targets-selected.pickle data/%/known.gff3
	$(gm) hits compare data/%/selected.pickle data/%/known.gff3 $@

data/%/compared.gff3 : $(gm) data/%/targets-compared.pickle
	$(gm) as gff3 data/%/targets-compared.pickle  $@

data/%/missing.gff3 : $(gm) data/%/targets-compared.pickle
	$(gm) comparisons select data/%/targets-compared.pickle type.pretty is missing - | $(gm).py as gff3 - $@

data/%/novel.gff3 : $(gm) data/%/targets-compared.pickle
	$(gm) comparisons select data/%/targets-compared.pickle type.pretty is novel - | $(gm).py as gff3 - $@

data/%/incorrect.gff3 : $(gm) data/%/targets-compared.pickle
	$(gm) comparisons select data/%/targets-compared.pickle type.pretty is incorrect - | $(gm).py as gff3 - $@

data/%/summary.csv : $(gm) data/%/targets-compared.pickle
	$(gm) comparisons summary data/%/targets-compared.pickle - | tee $@ | xsv table

data/%/inferred.json : $(gm) data/%/unknown-selected.pickle
	$(gm) as insertable data/%/unknown-selected.pickle $@
