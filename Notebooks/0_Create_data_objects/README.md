# 0_Create_data_objects

First stage of the analysis pipeline. These notebooks ingest the raw search results and sample annotations and produce the data objects (metadata table and AnnData files) consumed by every downstream folder. Run them in the order below.

## Notebooks

### [00_FASTA_overview.ipynb](00_FASTA_overview.ipynb) (informative)
Summarises the FASTA databases used for the searches. Parses the mixture FASTA (human, E. coli, yeast, contaminants) and the human-only FASTA, then reports per-species and contaminant sequence counts. Contaminant entries are identified by the `Cont_` prefix and species by the `OX=` taxonomy field.

### [01_metadata.ipynb](01_metadata.ipynb)
Builds the single sample-level metadata table used by all downstream notebooks. Reads the per-lab SDRF files and the cellenONE dispensing reports, merges them on plate and well position, applies parsing and imputation, and writes `Data/metadata/metadata.parquet`. Also generates QC overview figures (sample counts per lab and condition, cell morphology distributions, acquisition order) to verify the table before saving.

### [02_Results_to_AnnData_DIANN.ipynb](02_Results_to_AnnData_DIANN.ipynb)
Converts raw DIA-NN report Parquet files into AnnData `.h5ad` objects. Standardises modification notation to UNIMOD, constructs canonical precursor identifiers, joins search results with the metadata table, and pivots to precursor-level and protein-level matrices. Produces four objects (precursor and protein, for both the human samples and the quantitative QC mixtures), each with a `log2` layer and a contaminant flag.

### [02_Results_to_AnnData_Spectronaut.ipynb](02_Results_to_AnnData_Spectronaut.ipynb)
Same conversion as the DIA-NN notebook but for Spectronaut BGS Factory Report output. Produces the equivalent four AnnData objects plus the merged search-result Parquet files.

## Run order
`01_metadata.ipynb` must be run first, since both `02_*` notebooks join their search results against `metadata.parquet`. `00_FASTA_overview.ipynb` is independent and can be run at any time.
