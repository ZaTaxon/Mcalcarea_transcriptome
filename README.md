# De novo transcriptome assembly of *Macoma calcarea* (Gmelin, 1791) indicates patterns of adaptive evolution in the intertidal and subtidal bivalves

**Authors:** Evgeny Genelt-Yanovskiy, Sophia Nazarova, Alexey Masharskiy, Polina Drozdova and Olga Bondareva

**Abstract:** Infaunal bivalves are important components of marine soft-bottom ecosystems. In this data article we present the results of transcriptome sequencing of *Macoma calcarea* (Gmelin, 1791), a common clam species in the Subarctic and Arctic seas. Extracts of total RNA from two individuals from the White Sea were pooled to create a de novo assembly of the transcriptome. Over 9 million high-quality paired-end reads were assembled into 82387 transcripts, of which 27082 candidate protein-coding genes were annotated. We found 936 groups of orthologous genes between *M. calcarea* and four typical infaunal bivalve species, of which 141 families were single-copy orthologs. Functional annotation of the *M. calcarea* transcriptome allowed us to find candidate genes potentially involved in contrasting molecular adaptations to subtidal and intertidal habitats. This resource will facilitate adaptive genomic studies of Arctic marine species to test the potential response to climate stresses at various depth ranges.

This repository contains the following files and scripts:


**Files and folders:**


**Data** folder consist of protein fasta files for *Macoma calcarea* and four bivalve reassembled transcriptomes used in the analysis of orthologuous genes.

**calcarea_assembly** folder includes the *de novo* assembly of M.calcarea performed using Trinity v2.8.4 (TSA GIMC00000000) after removal of the remaining adapters and contaminants as suggested by the NCBI TSA piplenine.

**supplementary** - supplementary tables 2 & 3 for the paper in Marine Genomics except for functional annotation


**Scripts:**


**clams_proteinorto.R** - R script used for the gene family expansion and contraction analysis

**Parallel_substitution_looking_terminal_molluscs.py** - Python script for amino acid substitutions detection
