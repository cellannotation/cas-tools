# CLI Operations

Following operations are supported by the CAS commandline interface.

## Validate CAS file

Checks if the provided CAS data files comply with the specified CAS schema. In case of invalid files, the system logs the issues and throws an exception.

```commandline
cas validate --schema bican --data path/to/file
```

**Command-line Arguments:**
- `--schema`    : One of 'base', 'bican' or 'cap'. Identifies the CAS schema to validate data against.
- `--data`   : Path to the data file (or folder) to validate. If given path is a folder, validates all json files inside.

## Flatten onto AnnData

Flattens all content of CAS annotations to `obs` key:value pairs. Flattens all other content to key_value pairs in `uns`. The resulting AnnData object is then saved to a new file.

Key Features:
1. Parses command-line arguments for input JSON file, input AnnData file, and output file.
2. Reads and processes the input JSON file and AnnData file.
3. Updates the AnnData object with information from the JSON annotations and root keys.
4. Writes the modified AnnData object to a specified output file.

Detailed specification about the `flatten` operation can be found in the [related issue](https://github.com/cellannotation/cas-tools/issues/7).

```commandline
cas flatten --json path/to/json_file.json --anndata path/to/anndata_file.h5ad --output path/to/output_file.h5ad
```

**Command-line Arguments:**
- `--json`      : Path to the CAS JSON schema file.
- `--anndata`   : Path to the AnnData file. Ideally, the location will be specified by a resolvable path in the CAS file.
- `--output`    : Optional output AnnData file name. If provided a new flatten anndata file will be created,
    otherwise the inputted anndata file will be updated with the flatten data.

Please check the [related notebook](../notebooks/test_flatten.ipynb) to evaluate the output data format.

## Convert spreadsheet to CAS

Convert a spreadsheet to Cell Annotation Schema (CAS) JSON.

Detailed specification about the `spreadsheet2cas` operation can be found in the [related issue](https://github.com/cellannotation/cell-annotation-schema/issues/25).

```commandline
spreadsheet2cas --spreadsheet  path/to/spreadsheet_file --sheet optional_sheet_name
```

**Command-line Arguments:**
- `--spreadsheet` : Path to the spreadsheet file.
- `--sheet` : Target sheet name in the spreadsheet.
- `--output` : Output CAS file name (default: output.json).
 
## Convert AnnData to CAS

Convert an AnnData file to Cell Annotation Schema (CAS) JSON.

Detailed specification about the `anndata2cas` operation can be found in the 
[related issue](https://github.com/cellannotation/cas-tools/issues/10).

```commandline
cas anndata2cas --anndata path/to/anndata.h5ad --labelsets item1 item2 item3 --output path/to/output_file.json
```

**Command-line Arguments:**
- `--anndata` : Path to the AnnData file.
- `--labelsets` : List of labelsets.
- `--output` : Output CAS file name (default: output.json).
- `--hierarchy`: Flag indicating whether to include hierarchy in the output.

## Convert ABC to CAS

Converts given ABC cluster_annotation files to Cell Annotation Schema (CAS) JSON.

Detailed specification about the `abc2cas` operation can be found in the 
[related issue](https://github.com/cellannotation/cas-tools/issues/22).

```commandline
python -m cas abc2cas --catset path/to/cluster_annotation_term_set.csv --cat path/to/cluster_annotation_term.csv 
--output path/to/output_file.json
```

**Command-line Arguments:**
- `--catset` : Path to the Cluster Annotation Term Set file.
- `--cat` : Path to the Cluster Annotation Term file.
- `--output` : Output CAS file name (default: output.json).

## Convert CAS to ABC

**Status: Incomplete**

Converts given Cell Annotation Schema (CAS) to ABC files: cluster_annotation_term and cluster_annotation_term_set, 
and writes them to files with cat_file_path and cat_set_file_path.

Detailed specification about the `cas2abc` operation can be found in the 
[related issue](https://github.com/cellannotation/cas-tools/issues/22).

```commandline
python -m cas cas2abc --json path/to/json_file.json --catset path/to/cluster_annotation_term_set.csv --cat 
path/to/cluster_annotation_term.csv
```

**Command-line Arguments:**
- `--json` : Path to the CAS JSON schema file.
- `--catset` : Path to the Cluster Annotation Term Set file.
- `--cat` : Path to the Cluster Annotation Term file.

## Merge CAS to AnnData file

Integrates cell annotations from a CAS (Cell Annotation Schema) JSON file into an AnnData object.  It performs validation checks to ensure data consistency between the CAS file and the AnnData file.  The AnnData file location should ideally be specified as a resolvable path in the CAS file.

```commandline
cas merge --json path/to/CAS_schema.json --anndata path/to/input_anndata.h5ad --validate --output path/to/output.h5ad
```

**Command-line Arguments:**
- `--json`      : Path to the CAS JSON schema file.
- `--anndata`   : Path to the AnnData file. Ideally, the location will be specified by a resolvable path in the CAS file.
- `--validate`  : (Optional) Perform validation checks before writing to the output AnnData file.
- `--output`    : Output AnnData file name (default: output.h5ad).

Please check the [related notebook](../notebooks/test_merge.ipynb) to evaluate the output data format.

## Populate Cell IDs

Add/update CellIDs to CAS from matching AnnData file. Checks for alignment between `obs` key:value pairs in AnnData file and labelset:cell_label pairs in CAS for some specified list of `labelsets`. If they are aligned, updates `cell_ids` in CAS.

```commandline
cas populate_cells --json path/to/json_file.json --anndata path/to/anndata_file.h5ad --labelsets Cluster,Supercluster
```

**Command-line Arguments:**
- `--json`      : Path to the CAS JSON schema file.
- `--anndata`   : Path to the AnnData file. Ideally, the location will be specified by a resolvable path in the CAS file.
- `--labelsets` : (Optional) List of labelsets to update with IDs from AnnData. If value is not provided, rank '0' labelset is used.