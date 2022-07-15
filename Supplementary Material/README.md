# SARS-CoV-2 build definitions for Africa CDC

Build configuration for Africa CDC SARS-CoV-2 builds hosted at https://nextstrain.org/groups/africa-cdc.

## Instructions for running Africa CDC builds with GISAID data

### Setup Repositories

```
git clone https://github.com/nextstrain/ncov-africa-cdc.git
git clone https://github.com/nextstrain/ncov.git
cd ncov
cp -r ../ncov-africa-cdc/africa_cdc_profile .
```

### Download data

Login to [GISAID (gisaid.org)](https://www.gisaid.org/) and select the “EpiCoV” link from the navigation.

Select the “Downloads” link from the EpiCoV navigation bar. Scroll to the “Genomic epidemiology” section and select the “nextregions” button. Select the “Africa” button. Save the file as `hcov_africa.tar.gz` in the `ncov/data/` workflow directory.

Click "Back" to return to the main "Download" dialog, find the “Download packages” section, and select the “FASTA” button. Save the full GISAID sequences as `data/sequences_fasta.tar.xz`.

Select the “metadata” button from that same “Download packages” section and download the corresponding file as `data/metadata_tsv.tar.xz`.

### Filter data

From an existing `nextstrain` conda environment, install extra tools to extract data from GISAID files.

```
# Install tsvutils and UCSC command to extract sequences.
# You only need to do this once.
conda activate nextstrain
conda install -c conda-forge -c bioconda tsv-utils ucsc-fasomerecords
```

Extract African metadata and sequences from full GISAID downloads.

```
# Get metadata for Africa directly from tarball.
tar xzOf data/metadata_tsv.tar.xz metadata.tsv \
  | tsv-filter -H --str-in-fld Location:Africa \
  | xz -c -2 > data/metadata_africa.tsv.xz

# Get strain names for genomes.
# GISAID uses virus name, collection date, and submission date
# delimited by a pipe character.
xz -c -d data/metadata_africa.tsv.xz \
  | tsv-select -H -f 'Virus\ name','Collection\ date','Submission\ date' \
  | sed 1d \
  | awk -F "\t" '{ print $1"|"$2"|"$3 }' > data/strains_africa.txt

# Get genomes for strain names from tarball.
tar xzOf data/sequences_fasta.tar.xz sequences.fasta \
  | faSomeRecords /dev/stdin data/strains_africa.txt /dev/stdout \
  | xz -c -2 > data/sequences_africa.fasta.xz
```

### Run builds locally

Navigate to the `ncov` workflow directory; these instructions assume this is a sibling directory to this repository.
By default, the following command will run builds for all countries and regions.

```bash
nextstrain build \
  --cpus 4 \
  --memory 8Gib \
  . \
  --configfile ../ncov-africa-cdc/builds_africa.yaml
```

You can specify a comma-delimited subset of builds to run by build name.
The following example only runs builds for Ghana and Kenya.

```bash
nextstrain build \
  --cpus 4 \
  --memory 8Gib \
  . \
  --configfile ../ncov-africa-cdc/builds_africa.yaml \
  --config active_builds=ghana,kenya
```

### Run builds on AWS Batch

To run builds on AWS Batch, the build config and profiles directory must be in the `ncov` workflow directory.
Sync those files into `ncov`.

```bash
# Run from ncov/
rsync -arvz \
  ../ncov-africa-cdc/builds_africa.yaml \
  ../ncov-africa-cdc/africa_cdc_profile \
  .
```

Run all builds on AWS Batch.

```bash
nextstrain build \
  --aws-batch \
  --cpus 36 \
  --memory 70Gib \
  . \
  --configfile builds_africa.yaml
```

### Run builds on Terra

TBD.
