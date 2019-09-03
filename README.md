# CELLECT

**CELL**-type **E**xpression-specific integration for **C**omplex **T**raits (**CELLECT**) is a computational toolkit for identifying the probability of a cell type being involved in the genetic portion of complex traits. CELLECT uses a variety of models to assign genetic signal from a genome wide association study (GWAS) to single cell transcriptomic data. To integrate the two types of data, we must first calculate how important or specific each gene is to each cell type - the Expression Specificity (ES) of the cell type. Our group has developed a method for doing this called **CELL**-type **EX**pression specificity (**CELLEX**) that can be found here [add link].



## Installation

After cloning this repository there are a few other things that need to be set-up before CELLECT can be run.


One of the models used by CELLECT is LD score regression - it is vital that our forked version of this repository [(found here)](https://github.com/pascaltimshel/ldsc) is also cloned.

CELLECT is built as a Snakemake pipeline and the pipeline utilises conda environments. Therefore the easiest way to get started would be to [install miniconda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) (if conda is not already present on your system) and then either within base conda or a conda environment:
```bash
conda install -c bioconda -c conda-forge snakemake
conda install pandas
```

Within your cloned version of this repository, modify the config file `<src/snakemake/workflow_cts_config.yml>` so that it is specific to the system you are working on.

## Usage

To use CELLECT, specify the GWAS or list of GWAS you wish to provide as input for the genetic signal component in aforementioned config file. 
Also provide a multigeneset file which contains Expression Specificity values for each gene within all the cell types - this multigeneset file can be created by providing CELLEX with single cell count matrix and metadata.



### Example

To demonstrate CELLECT we have a short example that takes roughly 30 minutes to run on ... Our single cell input for this example is the [Mouse Nervous System](https://www.sciencedirect.com/science/article/pii/S009286741830789X) processed with CELLEX to get specificity scores then a subset of cell types taken to speed up the running time. The GWAS input is the UK Biobank's 1.1 million individual study on [Educational Attainment](https://www.nature.com/articles/s41588-018-0147-3)