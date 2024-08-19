# CAS Tools

CAS-Tools is a comprehensive utility package designed to facilitate the effective use and manipulation of the [Cell Annotation Schema](https://github.com/cellannotation/cell-annotation-schema) (CAS) in single-cell transcriptomics data analysis.  CAS supports recording the rationale and evidence for single cell annotation, including gene expression evidence and details of automated annotation transfer. The standard can be saved as a separate JSON file with a resolvable link to a martix (AnnData) file containing annotated data, or embedded in the annotated AnnData file. 

## Installation

You can install CAS-Tools [pypi package](https://pypi.org/project/cas-tools/) using `pip`:

```commandline
pip install cas-tools
```

## Overview

CAS-Tools simplifies the use of CAS by offering a set of programmatically accessible operations including:

- Validate Annotations: Validate annotations against the CAS standard to ensure compliance and consistency.
- Validation of markers against linked AnnData file
- Merge to Anndata - updating annotations in the AnnData file and saving JSON to the AnnData header.
- Reporting: Generation of dataframe reports of CAS content from JSON.
- Export to CAP-Anndata format: Merges into AnnData file following a derived, flattened representation of CAS, used by the Cell Annotation Platform (CAP) ([CAP-Anndata](https://github.com/cellannotation/cell-annotation-schema/blob/main/docs/cap_anndata_schema.md)). 
- Import from [Allen Brain Cell Atlas (ABC) format DataFrames](https://github.com/AllenInstitute/abc_atlas_access): 

## Getting Started

CAS-tools functionality can be accessed via imported object in python or via a command line tool. For CLI tool function details please see:

- [Command line interface](docs/cli.md)

