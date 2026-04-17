---
title: 'nowayout: An automated pipeline for taxonomic classification of Eukaryotic mitochondrial reads'
tags:
    - Bioinformatics
    - Food safety
    - Metagenomics
    - Mitochondrial DNA
    - Taxonomic classification
    - Bacterial Genomics
authors:
    - name: Kranti Konganti
      orcid: 0009-0002-3023-1099
      affiliation: "1"
      corresponding: true
    - name: Monica Pava-Ripoll
      orcid: 0000-0001-8401-4044
      affiliation: "1"
    - name: Amanda Windsor
      orcid: 0000-0002-5192-7047
      affiliation: "1"
    - name: Christopher Grim
      orcid: 0000-0002-7839-8309
      affiliation: "1"
    - name: Mark Mammel
      orcid: 0000-0002-8273-6091
      affiliation: "1"
    - name: Padmini Ramachandran
      orcid: 0000-0002-3958-3843
      affiliation: "1"
affiliations:
    - name: Human Foods Program, U.S. Food and Drug Administration, United States
      index: 1
      ror: 034xvzb47
date: 17 April 2026
bibliography: paper.bib
---

# Summary

`nowayout` is an ultra-fast automated software pipeline for taxonomic classification of eukaryotic mitochondrial reads from shotgun metagenomic data. Implemented in Nextflow, the pipeline employs a custom database to identify mitochondrial reads, and then performs taxonomic classification on those reads. The pipeline is currently being developed and evaluated for identifying arthropod taxa in food matrices, providing species-level resolution for both research-based applications and routine metagenomic analysis. Beyond arthropods, `nowayout` can also detect a broad range of eukaryotic DNA in shotgun metagenomic datasets. This allows for the analysis and verification of labeling claims when insects are used as food ingredients, as well as the identification of other eukaryotic taxa that may be present as food ingredients, making it a versatile tool for food safety research, future regulatory monitoring, and other applications where eukaryotic composition is of interest.

# Statement of Need

Food safety and quality assurance programs are increasingly utilizing sequencing and reference databases to improve the detection and identification of organisms in the food chain. This includes the evaluation of metagenomic methods to taxonomically resolve arthropod material from unavoidable crop-associated pests versus avoidable stored-product pests, a critical distinction for food safety assessment. In most conventional foods, DNA extracts are dominated by the food matrix (plant/animal), while arthropod DNA is often rare and fragmented. Consequently, recoverable arthropod DNA is well suited to targeted capture using mitochondrial probe (bait) panels. This strategy is effective because mitochondrial DNA occurs at high copy number, often persists when nuclear DNA is degraded, and is supported by extensive public reference databases for arthropod identification [@Foran; @Ratnasingham]. By contrast, in insect-based foods, the focus shifts from detecting trace arthropod DNA to distinguishing among arthropod species and identifying non-declared components. In both settings, practical analysis requires taxonomic classification of arthropods, whether present as contaminants or ingredients.

# State of the Field

Many marker- or k-mer based metagenomic tools can classify sequences down to the species rank [@mcintyre2017; @monicaetal], but their broader use in eukaryotic metagenomic analysis is limited by the lack of standardized, automated pipelines that minimize manual parameter tuning and support easy visualization of the results. Existing modular platforms such as QIIME 2 [@qiime2] support flexible taxonomic assignment workflows, but they have been primarily developed around amplicon-based microbial profiling, especially, prokaryotic 16S rRNA analysis. Although eukaryotic rRNA databases have recently become available [@eukaryome], they lack Cytochrome c oxidase subunit 1 (COX1) and plant chloroplast loci such as _rbcL_ and _matK_ from voucher specimens, two widely used protein-coding chloroplast genes that together serve as primary DNA barcoding loci for identification of land plant species [@plantdnabarcoding].

`nowayout` addresses this gap by providing a standardized, end-to-end workflow for mitochondrial read identification, taxonomic classification, filtering, and report generation. Rather than requiring users to assemble multiple tools into a custom workflow, `nowayout` emphasizes reproducibility, practical usability, and integrated reporting, including visual summaries that make dense taxonomic results easier to interpret across all input samples. To the best of our knowledge, `nowayout` is one of the first software tools for fully automated analysis of mitochondrial read identification and classification of eukaryotic species from shotgun metagenomic data.

# Methods and Materials

A brief overview of the `nowayout` pipeline is presented in \autoref{fig1}.

![The database preparation step produces the `mitomine` database which is used during KMA alignment stage of the main analysis.\label{fig1}](fig1.pdf){width=100%}

## Database Generation

The `nowayout` pipeline uses a custom database generated from sequences downloaded from NCBI GenBank [@genbank]. The database construction process begins by downloading the NCBI Taxonomy dump [@taxdump] and converting it to lineages using the `ncbitax2lin` tool [@ncbitax2lin]. All sequences catalogued as voucher or cytochrome c oxidase subunit 1 (COX1) are then downloaded from NCBI GenBank [@genbank; @baitspanel] using a keyword search, with separate catalogs maintained for each sequence type.

In the next stage, sequences less than 200 base pairs in length are filtered out of the dataset. CD-HIT [@cd-hit-est] is employed to perform sequence deduplication at 100% identity for each sequence catalog. For each filtered GenBank accession, the corresponding lineage is identified in the NCBI Taxonomy dump and assigned to create a comprehensive taxonomic reference.

In the final stage of database preparation, the voucher and COX1 catalogs are merged and subjected to a final deduplication step at 100% sequence identity using CD-HIT [@cd-hit-est]. This final catalog of GenBank sequences is then indexed using KMA [@kma], creating the main database for all subsequent taxonomic classification tasks. This custom database is named `mitomine`.

## Taxonomic Read Classification

The main analysis steps in `nowayout` (\autoref{fig1}) begin with metagenomic sequencing and read preprocessing using fastp [@fastp] for adapter trimming and quality filtering. Next, mitochondrial reads are identified by aligning to the `mitomine` database using KMA [@kma]. All mitochondrial reads are then extracted, and sketches are created for both the identified mitochondrial reads and accession hits. Reads identified as mtDNA are then classified using the `gather` command from sourmash [@sourmash]. Additionally, Salmon [@salmon] is used to bin the number of reads mapped to each lineage. Finally, a consolidated Krona [@krona] chart and an aggregated MultiQC [@multiqc] report are generated for all samples.

`nowayout` offers three preset threshold filters (strict, mild, relax) for exploring results and optimizing taxonomic classification specificity and stringency trade-offs. The pipeline is also available on the [HFP GalaxyTrakr](https://galaxytrakr.org) platform (version >= 25.x), providing a user-friendly web interface for researchers without command-line expertise.

## Results

To evaluate nowayout, we analyzed sequencing data (FASTQ files) generated from three mock genomic DNA (gDNA) mixtures comprising 23 insect taxa across seven orders (_Coleoptera_, _Blattodea_, _Diptera_, _Hymenoptera_, _Orthoptera_, _Lepidoptera_, and _Hemiptera_) combined with whole wheat flour gDNA as the food matrix. Mixture 1 contained one insect taxon, mixture 2 six, and mixture 3 twenty-two with staggered DNA inputs to mimic an uneven community. In all mixtures, whole wheat flour gDNA was added at four times the total insect gDNA mass (4:1). For each mixture, libraries were prepared in parallel and subjected to hybridization capture using two arthropod bait panel versions (v1 and v2) to enable direct panel-to-panel comparisons on identical mixture inputs.

The sequencing libraries were generated using the KAPA HyperPlus Library Preparation Kit (Roche Diagnostics) following the manufacturer’s instructions. For targeted hybridization capture, amplified libraries with similar concentrations were pooled and enriched for mitochondrial targets using custom arthropod bait panels, applying either panel v1 or panel v2 under the same capture workflow. Enriched libraries were sequenced on an Illumina MiSeq platform.

 Using `nowayout`, all expected insect taxa were detected in mixtures 1 (a) and 2 (b), and most taxa were recovered in mixture 3 (c) (21 of 23), with complete removal of wheat-flour background signal using panel v2 ([Table 1](https://research.foodsafetyrisk.org/nowayout/paper/Table_1.htm)). The nowayout visualizations produced a more interpretable and [informative taxonomic profile](https://research.foodsafetyrisk.org/nowayout/paper/krona_nowayout_paper.html).

# Research Impact Statement

`nowayout` increases the accessibility and reproducibility of food-based metagenomics by automating the identification and classification of eukaryotic mitochondrial reads. The custom `mitomine` database, derived from NCBI voucher specimens and COX1 sequences, provides a curated reference specifically tailored for mitochondrial analysis in complex food matrices. Although optimized for arthropod identification, the pipeline's ability to detect diverse eukaryotic taxa extends its utility to a broader range of food safety and quality control applications. By utilizing the Nextflow DSL2 framework, `nowayout` lowers the technical threshold for laboratories lacking in-house bioinformatics expertise, enabling the consistent application of advanced molecular methods through a standardized workflow with integrated, interpretable reporting.

The `nowayout` pipeline was validated using controlled mock community experiments consisting of three genomic DNA mixtures containing 23 insect taxa across seven orders spiked into a wheat flour matrix at a 4:1 flour-to-insect DNA ratio. The pipeline achieved 100% detection accuracy in mixtures 1 and 2, and identified 21 of 22 expected taxa in mixture 3 (\autoref{fig2}).

![The top taxa identified by `nowayout` in the three DNA mock Mixtures.\label{fig2}](fig2.pdf){width=100%}.

The reduction of wheat flour background signal observed with bait panel v2 highlights the analytical value of integrating high-performance bait panels with automated classification, which streamlines the analysis of complex metagenomic datasets across various samples and conditions. To improve accessibility, `nowayout` has been integrated into the [FDA's GalaxyTrakr platform](https://galaxytrakr.org) (version ≥25.x). This web-based interface removes the requirement for command-line expertise, allowing food safety laboratories without dedicated bioinformatics personnel to perform advanced mitochondrial metagenomic analysis.

# Software Design

The `nowayout` pipeline is implemented in Nextflow [@nextflow] following DSL2 principles and as such can be run on any UNIX based platform. All the individual steps are parallelized and run concurrently for all samples. `nowayout` is released under [MIT license](https://github.com/CFSAN-Biostatistics/nowayout/blob/main/LICENSE.md) and comprehensive documentation is hosted on [GitHub](https://github.com/CFSAN-Biostatistics/nowayout/blob/main/readme/nowayout.md). Current efforts are being undertaken to develop custom algorithms and classification methods to better handle ambiguous read assignments. Additionally, expanding to support Oxford Nanopore long reads in addition to the current Illumina short read capability is planned, enabling analysis of a wider range of sequencing data types.

# AI Usage Disclosure

Generative AI was not used in any aspects of the development of the software, in writing the documentation or during any aspects of the paper authorship process.

# Acknowledgements

We would like to thank the HPC team for providing systems administration support for the HFP Reedling HPC Cluster.

## References
