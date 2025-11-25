# UHVDB/proteinfunction
A Nextflow wrapper for predicting the function of proteins predicted from metagenome assemblies.

### Overview
This wrapper performs the following steps:

1. Remove partial proteins from the `--input_faa` file
2. Calculate the sequence hash and length of proteins in the `--input_faa` file
3. Extract unique protein sequences based on the hash + length
4. Cluster dereplicated protein sequences at `--min_id` % amino acid identity and `--min_cov` bidirectional coverage
5. Download UniRef90 FAA file and creates an MMSeqs2 database
6. Align protein cluster representatives to UniRef90
7. Extract unaligned proteins and align to InterPro with HMMER
8. 

### Quick start
In addition to automated downloads and cleanup (limiting disk requirements), this wrapper also makes setup very easy.

First, install Conda/Mamba/Micromamba/Pixi
*Example Micromamba installation*
```
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
```

Then create a Nextflow environment
```
micromamba create -n nextflow -c conda-forge -c bioconda nextflow singularity -y
```

Activate the Nextflow environment
```
micromamba activate nextflow
```

Then just run the pipeline!
```
nextflow run UHVDB/proteinfunction -profile test,<docker/singularity/conda/mamba>
```

### Usage
The only arguments for this tool are:

`--query_fna`: Path (or URL) to a FNA file of query virus genomes.

`--vmr_dmnd`: Path (or URL) to a DMND database created from ICTV genomes (or any set of genomes).

`--vmr_url`: Path (or URL) to an ICTV VMR Excel file.

`--email`: Email address to use for Entrez downloads.

`--chunk_size`: Number (Integer) of host sequences contained in each chunk (default: 10,000)

`--diamond_args`: CLI arguments to use when running DIAMOND query-v-ictv (default: `--masking none -k 1000 -e 1e-3 --very-sensitive`)

`--min_score`: Minimum protein similarity values to output (default: 5.5, decreasing this to 0 will dramatically increase the output file size)

`--output`: Path to output TSV file containing protein similarity values for all queries vs ICTV viruses.

### Output
The only output of this wrapper is a protein similarity TSV file with 3 columns (`query`,`reference`,`protein_similarity`)
```tsv
NC_024375.1	Cinqassovirus_ah1--MG250483.1	0.49
NC_024375.1	Krischvirus_jse--EU863408.1	0.73
NC_024375.1	Krischvirus_georgiaone--EF437941.1	0.25
```

### Credits
This wrapper was made by @CarsonJM. However, primary credit of course goes to the UHGV-classifier developers (https://github.com/snayfach/UHGV-classifier), and their work should be cited if this wrapper is used (https://doi.org/10.1101/2025.11.01.686033).