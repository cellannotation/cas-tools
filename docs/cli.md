# CLI Operations

Following operations are supported by the CAS commandline interface.

## Flatten onto AnnData

Processes and integrates information from a JSON file and an AnnData (Annotated Data) file,  creating a new AnnData object that incorporates the metadata. The resulting AnnData object is then saved to a new file.

Key Features:
1. Parses command-line arguments for input JSON file, input AnnData file, and output file.
2. Reads and processes the input JSON file and AnnData file.
3. Updates the AnnData object with information from the JSON annotations and root keys.
4. Writes the modified AnnData object to a specified output file.

```
cas flatten --json path/to/json_file.json --anndata path/to/anndata_file.h5ad --validate --output path/to/output_file.h5ad
```

**Command-line Arguments:**
- `--json`      : Path to the CAS JSON schema file.
- `--anndata`   : Path to the AnnData file. Ideally, the location will be specified by a resolvable path in the CAS file.
- `--validate`  : (Optional) Perform validation checks before flattening to the output AnnData file.
- `--output`    : Output AnnData file name (default: output.h5ad).

Please check the [related notebook](../notebooks/test_flatten.ipynb) to evaluate the output data format.

## Merge CAS to AnnData file

Integrates cell annotations from a CAS (Cell Annotation Schema) JSON file into an AnnData object.  It performs validation checks to ensure data consistency between the CAS file and the AnnData file.  The AnnData file location should ideally be specified as a resolvable path in the CAS file.

```
cas merge --json path/to/CAS_schema.json --anndata path/to/input_anndata.h5ad --validate --output path/to/output.h5ad
```

**Command-line Arguments:**
- `--json`      : Path to the CAS JSON schema file.
- `--anndata`   : Path to the AnnData file. Ideally, the location will be specified by a resolvable path in the CAS file.
- `--validate`  : (Optional) Perform validation checks before writing to the output AnnData file.
- `--output`    : Output AnnData file name (default: output.h5ad).

Please check the [related notebook](../notebooks/test_merge.ipynb) to evaluate the output data format.