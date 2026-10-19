# Paths and settings used by every script in this directory.
# Edit the three paths at the top, then run the scripts in numerical order.

# Raw HiFi reads, one unaligned BAM per SMRT Cell, with kinetic tags intact.
# These are deposited in the SRA; replace the accession below once assigned.
#   SRA accession: PENDING
READS_CELL1=../data/m84378_261006_192612_s3.hifi_reads.bc2020.bam
READS_CELL2=../data/m84378_261013_200442_s4.hifi_reads.bc2020.bam

# Where intermediate and output files are written.
WORK=../work

# Bakta database v6.0, downloaded separately. See docs/databases.md.
BAKTA_DB=/path/to/bakta_db/db

# SMRT Link command line tools, version 10.1.0.119588, from PacBio.
# The conda environment captures only the wrapper package; the tools
# themselves come from the PacBio installer. See docs/databases.md.
SMRTLINK=/path/to/smrtlink/smrtcmds/bin

THREADS=20

mkdir -p "$WORK"

# conda activate is not available in a non-interactive shell until the
# hook is sourced. Adjust the path if conda is installed elsewhere.
source "$(conda info --base)/etc/profile.d/conda.sh"

# Kinetic model used for modification calling. SMRT Link 10.1 ships
# four models, none matching Revio chemistry; SP3-C3 is the closest
# and is specified explicitly because automatic selection fails.
# Find it under the kineticsTools resources directory of the install.
KINETICS_MODEL=/path/to/smrtlink/.../kineticsTools/resources/SP3-C3.npz.gz
