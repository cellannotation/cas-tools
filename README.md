# CAS Tools

**Note: CAS-Tools is currently under active development. Its behavior and interfaces may change as new features are added and existing ones are refined. Use with caution in production environments.**

CAS-Tools is a comprehensive utility package designed to facilitate the effective use and manipulation of the Cell Annotation Schema (CAS) in single-cell transcriptomics data analysis.

## Overview

The Cell Annotation Schema (CAS) addresses the inherent variability in annotating cell types/classes in single-cell transcriptomics datasets. While annotations lack the reasoning behind the choice of labels, CAS enables the recording of additional metadata about individual cell type annotations. This includes marker genes used as evidence and details of automated annotation transfer.

CAS-Tools simplifies the utilization of CAS by offering a set of programmatically accessible operations to:

- Integration with AnnData: Flatten CAS onto observations in AnnData format for seamless integration and analysis.
- Flatten Onto Dataframes: Decompose the schema into individual tables suitable for use in dataframes/TSVs.
- Validate Annotations: Validate annotations against the CAS standard to ensure compliance and consistency.

## Installation

You can install CAS-Tools [pypi package](https://pypi.org/project/cas-tools/) using `pip`:

```commandline
pip install cas-tools
```

## Getting Started

Please see related guides:

- [Command line interface](docs/cli.md)

