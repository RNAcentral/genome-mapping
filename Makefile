gm=bin/gm.py
all_md5=data/md5.tsv

all : pombe

###############################################################################
# pombe
###############################################################################

pombe : data/pombe/inferred.json data/pombe/summary.csv data/pombe/unknown.fasta
	@echo "Found `jq 'length' data/pombe/inferred.json` coordinates of `grep -c '^>' data/pombe/unknown.fasta`"

data/pombe/genome.fasta : bin/fetch
	$^ "ftp://ftp.ensemblgenomes.org/pub/current/fungi/fasta/schizosaccharomyces_pombe/dna/Schizosaccharomyces_pombe.ASM294v2.dna.toplevel.fa.gz" > $@

data/pombe/known.gff3 : bin/fetch
	$^ "ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/6.0/genome_coordinates/Schizosaccharomyces_pombe.ASM294v2.gff3.gz" > $@

data/pombe/raw.bed : bin/fetch-bed
	$^ "ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/6.0/genome_coordinates/Schizosaccharomyces_pombe.ASM294v2.bed.gz" > $@



###############################################################################
# Generic Definitions
###############################################################################

$(all_md5) : bin/fetch
	$^ "ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral//releases/6.0/md5/md5.tsv.gz" > $@
	xsv index $@

data/%/known-md5.txt : bin/known-md5s $(all_md5) data/%/known.gff3 data/%/bad-ids
	$^ > $@

data/%/targets.bed : bin/targets-bed data/%/raw.bed data/%/known.gff3 data/%/bad-ids
	$^ > $@

data/%/targets.fasta : bin/targets data/%/genome.fasta data/%/targets.bed data/%/md5sums.txt data/%/known-md5.txt
	$^ > $@

data/%/targets-hits.pickle : $(gm) data/%/genome.fasta data/%/targets.fasta
	$(gm) find $(word 2,$^) $(word 3,$^) $@

data/%/unknown-hits.pickle : $(gm) data/%/genome.fasta data/%/unknown.fasta
	$(gm) find $(word 2,$^) $(word 3,$^) $@

data/%/targets-compared.pickle : $(gm) data/%/targets-selected.pickle data/%/known.gff3
	$(gm) hits compare $(word 2,$^) $(word 3,$^) $@

data/%/compared.gff3 : $(gm) data/%/targets-compared.pickle
	$(gm) as gff3 $(word 2,$^) $@

data/%/missing.gff3 : $(gm) data/%/targets-compared.pickle
	$(gm) comparisons select $(word 2,$^) type.pretty is missing - | $(gm).py as gff3 - $@

data/%/novel.gff3 : $(gm) data/%/targets-compared.pickle
	$(gm) comparisons select $(word 2,$^) type.pretty is novel - | $(gm).py as gff3 - $@

data/%/incorrect.gff3 : $(gm) data/%/targets-compared.pickle
	$(gm) comparisons select $(word 2,$^) type.pretty is incorrect - | $(gm).py as gff3 - $@

data/%/summary.csv : $(gm) data/%/targets-compared.pickle
	$(gm) comparisons summary $(word 2,$^) - | tee $@ | xsv table

data/%/targets-selected.pickle : $(gm) data/%/targets-hits.pickle data/%/select-spec.json
	$(gm) hits select-using-spec $(word 2,$^) $(word 3,$^) $@

data/%/unknown-selected.pickle : $(gm) data/%/unknown-hits.pickle data/%/select-spec.json
	$(gm) hits select-using-spec $(word 2,$^) $(word 3,$^) $@

data/%/inferred.json : $(gm) data/%/unknown-selected.pickle
	$(gm) as insertable $(word 2,$^) $@

.PHONY : all
