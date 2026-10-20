# Bacterial pan-immunity emerges from defense–growth trade-offs and herd protection

Support scripts for this project.


## Requirements

- Linux
- Conda / Mamba environments per tool: SPAdes, GTDB-Tk, Prodigal, DefenseFinder, QIIME 2, BLAST+, ChronoStrain
- R >= 4.2 with `deSolve`, `FME`, `tidyverse`, `foreach`, `doParallel`
- Python >= 3.8 with `Biopython`

## Repository layout

```
Section_1/   Strain genomes: assembly, taxonomy, gene prediction, defense systems
Section_2/   Full-length 16S rRNA gene profiling of synthetic communities
Section_3/   Dynamical models of defense-growth trade-offs and herd protection
Section_4/   Strain-resolved time-series inference for Thauera aminoaromatica
```

Scripts are numbered in execution order and each carries a header describing its tool, inputs and outputs.

## Sections

### 1) Strain genomes: assembly, taxonomic annotation and defense systems

- **01_Genome_assembly.sh** Assemble the whole-genome sequencing reads of each isolate with [SPAdes](https://github.com/ablab/spades).
- **02_Taxonomic_annotation.sh** Assign species-level taxonomy to each assembly with [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk).
- **03_Gene_prediction.sh** Predict open reading frames with [Prodigal](https://github.com/hyattpd/Prodigal).
- **04_Defense_system_annotation.sh** Annotate antiviral defense systems in the predicted proteins with [DefenseFinder](https://github.com/mdmparis/defense-finder).

### 2) Full-length 16S rRNA gene profiling of synthetic communities

- **01_Full_length_16S.sh** Import PacBio CCS reads, denoise with [DADA2](https://benjjneb.github.io/dada2/) (`denoise-ccs`) in [QIIME 2](https://qiime2.org).

### 3) Dynamical models of defense–growth trade-offs and herd protection

- **01_Monod_growth_fitting.R** Fit a resource-dependent (Monod) growth model to the phage-free monoculture growth curves of each strain.
- **02_Model_1_12h_simulation.R** Simulate Model_1 (defense–growth trade-off only): nutrient-dependent growth plus phage-induced lysis of the susceptible strain.
- **03_Model_2_12h_simulation.R** Simulate Model_2: Model_1 extended with adsorption-mediated phage removal by the resistant strain.
- **04_Chemostat_parameter_scan.R** Scan defense cost against adsorption rate in parallel over a 100 × 100 grid under a chemostat extension of Model_2.


### 4) Strain-resolved time-series inference for *Thauera aminoaromatica*

Database construction:

- **01_Split_metaphlan_markers.py** Split the MetaPhlAn species-level marker FASTA into one file per marker and write the [ChronoStrain](https://github.com/gbaym/chronostrain) seed table.
- **02_Marker_prevalence_QC.sh** Align every marker against each MAG of the species-level cluster with BLASTN.
- **03_Summarize_marker_blast.py** Summarise the alignments into per-MAG coverage, per-marker prevalence and a presence/absence matrix.
- **04_Build_reference_catalog.py** Join the contigs of each MAG into one pseudo-chromosome and write the reference catalog and accession map, using a MAG-to-lineage table supplied by the user.
- **05_Build_marker_subsets.py** Write the seed table, ID list and ID map of the conserved marker set.
- **06_Make_ChronoStrain_database.sh** Write the ChronoStrain configuration, build the database from the conserved-marker subset, and cluster its entries at marker identity 0.998 to define the lineages tracked in the time series.

Time series:

- **07_Prepare_time_series_manifest.py** Locate the paired FASTQs of the monthly samples and write the manifest.
- **08_Count_read_pairs.sh** Count the read pairs of each sample.
- **09_Build_chronostrain_reads.py** Combine the manifest and the read depths into the ChronoStrain read table and a sample metadata table.
- **10_Filter_reads.sh** Run the ChronoStrain read filter and obtain the filtered read table.
- **11_Fix_pair_ids.py** Restore the `/1` and `/2` mate suffixes lost during filtering, rewrite the read table accordingly, and report per-sample pairing statistics.
- **12_Run_advi.sh** Infer the posterior relative abundance of each lineage in every month with ChronoStrain's ADVI estimator.
