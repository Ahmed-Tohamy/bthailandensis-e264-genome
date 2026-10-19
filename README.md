# Burkholderia thailandensis E264 genome and methylome

Analysis code for the closed genome and base modification profile of
*B. thailandensis* E264.

Ahmed Al-Tohamy and Anne Grove, Louisiana State University.

## Data

Sequencing reads are in the SRA under accession SRXXXXXXX, as two
unaligned BAM files, one per SMRT Cell. The format matters: the
per-base polymerase kinetics that modification calling depends on are
stored as tags that do not survive conversion to FASTQ.

The assembly is in GenBank under accession CPXXXXXX (chromosome 1) and
CPXXXXXX (chromosome 2). The per-site modification calls are archived at
Zenodo, DOI 10.5281/zenodo.XXXXXXX.

Published result tables are in `tables/`, so a rerun can be checked
without downloading anything from the journal.

## Layout

```
scripts/       one script per Methods section, in order
environment/   conda environment files
tables/        result tables as published
verify/        output of the verification runs described below
```

## Running the scripts

Edit the paths at the top of `scripts/config.sh`, then run the scripts
in numerical order. Each reads what the previous ones produced.

```
cd scripts
bash 01_assemble.sh
bash 02_annotate.sh
bash 03_qc.sh
bash 04_call_modifications.sh
bash 05_find_motifs.sh
bash 07_compare_assemblies.sh
bash 09a_select_methyltransferases.sh
bash 09_methyltransferases.sh
```

Four steps take arguments rather than reading the config, because they
act on a region or a file the user chooses. They need an environment
with pysam; `ipd_env` has it.

```
python3 06_oninstrument_calls.py ../work/assembly.fasta \
    ../work/merged_unaligned.bam ../work/aligned.bam \
    ../work/oninstrument_recovery.tsv

python3 08_classify_variants.py ../work/assembly.fasta \
    --snps 2005=../work/compare/variants_GCF_000012365.1.snps \
           2014=../work/compare/variants_GCA_003568605.1.snps \
    --output ../work/variants_classified.tsv

python3 10_copy_number.py ../work/assembly.fasta ../work/depth.txt \
    ../work/aligned.bam --replicon chromosome_1 \
    --start 1905085 --end 2113703 --output ../work/copy_number.tsv

python3 11_locus_crosswalk.py ../work/rebase/hits.tsv \
    ../work/compare/GCF_000012365.1.gff ../work/locus_crosswalk.tsv
```

Expect a few hours on 20 cores, most of it in assembly and
modification calling. Peak memory is around 32 GB.

## Software

Conda environment files are in `environment/`:

```
conda env create -f environment/flye_env.yml
```

Two things are not on conda and have to be fetched by hand.

SMRT Link 10.1.0.119588 comes from PacBio. The yml file for it is
nearly empty, since conda only installed a wrapper and the tools
themselves live inside the tarball. Scripts 04 and 05 call
`ipdSummary`, `motifMaker`, `pbmm2` and `ccs-kinetics-bystrandify`,
all of which come from there.

The Bakta database v6.0 is a 70 GB download (DOI
10.5281/zenodo.14916843). Use that version rather than a newer one if
you want the same annotation, since the ISButh2 calls come from the
ISFinder release bundled with it.

## Verification

The four Python scripts were rerun from this repository against the
original inputs and compared with the published tables.

| Script | Published | Rerun |
|---|---|---|
| 06 on-instrument calls | 1,184,624 calls placed; CACAG 100.0%, GTWWAC 98.9%, GAGHNNNCTC 99.9% | identical |
| 08 variant classification | 30 against 2005, 71 against 2014, none shared | identical |
| 11 locus crosswalk | BTH_I0041, BTH_II1029, both at 100% identity | identical |
| 10 copy number | GC corrected 1.64, raw 1.83, 18 of 49 reads at the junction | 1.63, 1.85, 17 of 48 |

The small differences in script 10 come from three choices made when
reimplementing it, each documented in the script: the background is
taken from flanking windows rather than two fixed flanking regions,
the raw ratio is a ratio of mean depths rather than a median of
per-window ratios, and junction reads are sampled at positions derived
from the segment end rather than at round coordinates. The quantity
the paper reports, the GC corrected copy number, agrees.

Output of these runs is in `verify/`.

The shell scripts were not rerun, since that would mean reassembling
and recalling the genome. Every command and version in them was taken
from the logs of the runs that produced the published results.

## A note on the kinetic model

Modification calling used a kinetic model built for a different
sequencing chemistry. The available version of kineticsTools predates
the Revio platform and ships no model for its chemistry; with
automatic selection it refuses to run. A Sequel II model was supplied
explicitly instead.

The interpulse duration ratios and quality scores in the deposited
calls are therefore not calibrated in absolute terms. The motifs do
not depend on that calibration, having been recovered independently by
two discovery methods and checked against the instrument's own
modification calls, but anyone comparing these kinetic values against
a chemistry-matched run should know where they came from.

## Two coordinate frames

The genome was annotated twice, once as assembled and once after
reorientation to start at *dnaA*. Reorientation rotates both replicons
and reverse complements one, so the two annotations share no
coordinates and carry different locus tags. Everything in the paper
uses the unreoriented frame. Script 11 builds the crosswalk between
our locus tags, the RefSeq tags of the 2005 reference, and the
original BTH tags the literature uses; the result is
`tables/locus_crosswalk.tsv`.

## Citation

Al-Tohamy A, Grove A. PacBio HiFi genome assembly and methylome of the
Burkholderia thailandensis type strain E264. *Scientific Data* (in review).

## Contact

Ahmed Al-Tohamy, Ahmed.Al-Tohamy@lsu.edu
