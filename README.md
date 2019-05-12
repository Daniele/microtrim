# Microtrim

Trimmer for microRNA sequences.




## Usage
The scripts require `python 3.6` and it using microconda is suggested.

Install all the requirements with the command:

```
conda install --yes --file requirements.txt
```

Trim and align a microRNA dataset with the following:

```
python3 microtrim.py -i SRR8311267.fastq        \ 
                     -o trimmed.fastq           \
                     -m ndleven                 \      
                     -a TGGAATTCTCGGGTGCCAAGG   \
                     --max-distance .3          \
                     --match-only 13            \
                     --trim-to 23               \
                     --workers 16               \
                     --aligner bowtie_htseq
```

Check all the available option with `python3 microtrim.py --help`.


## Datasets

**MicroRNA**:
[ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR831/007/SRR8311267/SRR8311267.fastq.gz](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR831/007/SRR8311267/SRR8311267.fastq.gz)

**Human genome**: [ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz](ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz)


## License

MIT License