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

Here's your updated **Markdown documentation** to align with the **`export2cap`** naming convention:

---

## **Export CAS to CAP Format in AnnData**

`export2cap` converts CAS annotations into `obs` key-value pairs and stores other CAS content as key-value pairs in `uns`. The resulting AnnData object is then saved to a new file.

### **Key Features:**
1. Parses command-line arguments for an optional input JSON file and/or an AnnData file.
2. **Requires at least one of the following to be supplied:** a CAS JSON file (`--json`) or an AnnData file (`--anndata`).
3. If the CAS JSON file is not provided via `--json`, the tool expects to find the CAS JSON embedded in the AnnData file’s `uns` section.
4. If the AnnData file is not provided via `--anndata`, it will be downloaded using the matrix file ID from the CAS JSON.
5. Updates the AnnData object with information from the CAS JSON annotations and root keys.
6. Writes the modified AnnData object to a specified output file.

A detailed specification about the `export2cap` operation can be found in the [related issue](https://github.com/cellannotation/cas-tools/issues/7).

```commandline
cas export2cap --json path/to/json_file.json --anndata path/to/anndata_file.h5ad --output path/to/output_file.h5ad
```

### **Command-line Arguments:**
- `--json`      : **Optional** path to the CAS JSON schema file.  
  - *If not provided*, the CAS JSON is expected to be embedded in the AnnData file’s `uns` section.
- `--anndata`   : **Optional** path to the AnnData file.  
  - *If not provided*, the AnnData file will be downloaded using the matrix file ID from the CAS JSON.
- `--output`    : Optional output AnnData file name.  
  - If provided, a new AnnData file with CAS exported to CAP format will be created; otherwise, the input AnnData file will be updated in place.
- `--fill-na`   : Optional boolean flag indicating whether to fill missing values in the `obs` field with `pd.NA`.  
  - If provided, missing values will be replaced with `pd.NA`.  
  - If not provided, missing values will remain as empty strings.

**Note:** At least one of `--json` or `--anndata` must be supplied. Additionally, if the CAS JSON is not provided via `--json` and the AnnData file does not contain the CAS JSON in its `uns` section, the operation will fail with an error indicating that the CAS JSON is missing.

Please check the [related notebook](../notebooks/test_export2cap.ipynb) to evaluate the output data format.

---

## Unflatten Operation

Unflattens all content of a flattened AnnData file into a CAS JSON file and creates an unflattened AnnData file.

**Key Features:**
1. Parses command-line arguments for the input AnnData file and optional JSON and output files.
2. Processes the input AnnData file and, optionally, a JSON file.
3. Converts flattened AnnData content back to its unflattened version and creates corresponding CAS JSON files.
4. Annotation Verification and Update:  
   - Uses a lookup dictionary (stored in the `uns` section of the AnnData file and generated in the export2cap step) to verify and update annotations.  
   - **Direct Update:** Annotations are updated when both the labelset-label pair and the generated cell hash (computed using labelset labels and `cell_ids`) match.  
   - **Discarding Mismatches:** If the labelset-label pair matches but the hashes do not, the annotation is discarded.  
   - **Handling Label Changes:** If the cell hash matches without a matching labelset-label pair, it suggests a possible label change, and the annotation is updated accordingly.
5. Saves the unflattened AnnData and CAS JSON files to the specified output locations.

```commandline
cd src
python -m cas unflatten --anndata path/to/anndata_file.h5ad --json path/to/json_file.json --output_anndata path/to/output_file.h5ad --output_json path/to/output_cas.json
python -m cas unflatten --anndata path/to/anndata_file.h5ad
```

**Command-line Arguments:**
- `--anndata`        : Path to the input AnnData file that contains flattened data.
- `--json`           : Optional path to the CAS JSON file. If provided, the 'annotations' within the file will be updated based on lookup dictionary checks; if not provided, a new CAS JSON file will be created.
- `--output_anndata` : Optional output AnnData file name. If not provided, `unflattened.h5ad` will be used as the default name.
- `--output_json`    : Optional output CAS JSON file name. If not provided, `cas.json` will be used as the default name.

### Usage Example

To execute the `unflatten` operation, use the following command:

```commandline
cd src
python -m cas unflatten --anndata path/to/anndata_file.h5ad --json path/to/json_file.json --output_anndata path/to/output_file.h5ad --output_json path/to/output_cas.json
```

## Convert spreadsheet to CAS

Convert a spreadsheet to Cell Annotation Schema (CAS) JSON.

Detailed specification about the `spreadsheet2cas` operation can be found in the [related issue](https://github.com/cellannotation/cell-annotation-schema/issues/25).

```commandline
spreadsheet2cas --spreadsheet  path/to/spreadsheet_file --sheet optional_sheet_name
```

**Command-line Arguments:**
- `--spreadsheet` : Path to the spreadsheet file.
- `--sheet` : Target sheet name in the spreadsheet.
- `--anndata` : Path to the AnnData file. If not provided, AnnData will be downloaded using CxG LINK in spreadsheet.
- `labelsets` : List of names of observation (obs) fields used to record author cell type names,
    which determine the rank of labelsets in a spreadsheet. If not provided, ranks will be determined based on
    the order of the fields specified in the CELL LABELSET NAME column.
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
- `--labelsets` : List of labelsets, which are names of observation (obs) fields used to record author cell type
    names. The labelsets should be provided in order, starting from rank 0 (leaf nodes) and ascending to higher ranks.
- `--output` : Output CAS file name (default: output.json).
- `--hierarchy`: Flag indicating whether to include hierarchy in the output.
- `--accession_columns`: List of columns in the AnnData obs that contain accession ID information.
            This list should match the order and length of the labelsets argument.
            If not provided, accession IDs will be automatically generated using a hash of the cells in each cell set.
            Defaults to None.

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

Integrates cell annotations from a CAS (Cell Annotation Schema) JSON file into an AnnData object. It performs validation checks to ensure data consistency between the CAS file and the AnnData file. The AnnData file location should ideally be specified as a resolvable path in the CAS file.

```commandline
cas merge --json path/to/CAS_schema.json --anndata path/to/input_anndata.h5ad --validate --output path/to/output.h5ad
```

**Command-line Arguments:**
- `--json`      : Path to the CAS JSON schema file.
- `--anndata`   : Path to the AnnData file. If not provided, AnnData will be downloaded using the matrix file ID from the CAS JSON.
- `--validate`  : (Optional) If set, the following validation checks will be performed before writing to the output AnnData file:
    1. Verifies that all cell barcodes (cell IDs) in CAS exist in AnnData and vice versa.
    2. Identifies matching labelset names between CAS and AnnData.
    3. Validates that cell sets associated with each annotation match between CAS and AnnData.
    4. Checks if the cell labels are identical; if not, provides options to update or terminate.
- `--output`    : Output AnnData file name (default: output.h5ad).

Please check the [related notebook](../notebooks/test_merge.ipynb) to evaluate the output data format.

## Populate Cell IDs

Add/update CellIDs to CAS from a matching AnnData file. Checks for alignment between `obs` key-value pairs in the AnnData file and labelset:cell_label pairs in CAS for a specified list of `labelsets`. If they are aligned, updates `cell_ids` in CAS.

```commandline
cas populate_cells --json path/to/json_file.json --anndata path/to/anndata_file.h5ad --labelsets Cluster Supercluster --validate
```

**Command-line Arguments:**
- `--json`      : (Required) Path to the CAS JSON schema file.
- `--anndata`   : (Required) Path to the AnnData file. Ideally, the location will be specified by a resolvable path in the CAS file.
- `--labelsets` : (Optional) A space-separated list of labelsets to update with IDs from AnnData. If not provided, the labelset with rank '0' is used by default. The labelsets should be provided in hierarchical order, starting from rank 0 (leaf nodes) and ascending to higher ranks.
- `--validate`  : (Optional) If set, strict validation is enforced. If validation fails, the program exits immediately with an error code (`sys.exit(1)`). Otherwise, it logs warnings but continues execution.

**Usage Examples:**

Run without validation (default mode):
```commandline
cas populate_cells --json cas.json --anndata data.h5ad --labelsets Cluster Supercluster
```

Run with strict validation (`--validate`):
```commandline
cas populate_cells --json cas.json --anndata data.h5ad --labelsets Cluster Supercluster --validate
```

## Convert CAS data to RDF

Converts the given CAS data to RDF format.

```commandline
cas cas2rdf --schema bican --data path/to/file.json --ontology_ns MTG --ontology_iri https://purl.brain-bican.org/ontology/AIT_MTG/ --out path/to/output.rdf --exclude_cells
```

**Command-line Arguments:**
--schema    : (Optional) Name of the CAS release (such as one of `base`, `cap`, `bican`) or path to the
                CAS schema file or url of the schema file. If not provided, reads the `base` CAS schema from the cas module.
--data   : Path to the json data file
--ontology_ns    : Ontology namespace (e.g. `MTG`)
--ontology_iri    : Ontology IRI (e.g. `https://purl.brain-bican.org/ontology/AIT_MTG/`)
--out    : The output RDF file path.
--skip_validate    : (Optional) Determines if data-schema validation checks will be performed. Validations are performed by default.
--exclude_cells    : (Optional) Determines if cell data will be included in the RDF output. Cell data is exported to RDF by default.

## Add Author Annotations to CAS JSON

This tool processes input CSV and CAS JSON files to add annotation fields to the CAS JSON based on matching columns specified by the user. It can optionally use selected columns from the CSV for annotations and outputs the annotated CAS JSON to a specified file. If no specific columns are provided, all columns from the CSV file will be used.

### Command-line Arguments:
- **--cas_json**: Path to the CAS JSON file that will be updated with annotations. This parameter is required.
- **--csv**: Specifies the path to the CSV file containing the data for annotation. This parameter is required.
- **--join_on**: Specifies the single column name in the CSV used for matching records. Each row must have a unique value in this column.
- **--join_on_cell_set_id**: Use 'cell_set_id' as the column for matching records. This option is triggered with a flag.
- **--join_on_labelset_label**: Use a pair of 'labelset', 'cell_label' columns for matching records. This option is triggered with a flag.
- **--columns**: Optionally specifies which columns in the CSV will be used for annotations. If not provided, all columns are used. Column names containing spaces must be enclosed in quotes (e.g., `"Column Name"`).
- **--output**: Specifies the output file name for the annotated CAS JSON. Defaults to `output.json`.

### Usage Examples:
```commandline
cd src
python -m cas add_author_annotations --cas_json path_to_cas.json --csv path_to_csv --join_on CrossArea_cluster --columns random_annotation_x random_annotation_y --output annotated_output.json
python -m cas add_author_annotations --cas_json path_to_cas.json --csv path_to_csv --join_on_cell_set_id --output annotated_output.json
python -m cas add_author_annotations --cas_json path_to_cas.json --csv path_to_csv --join_on_labelset_label --output annotated_output.json
```

## Split CAS JSON by Cell Label

This tool allows you to split a CAS JSON file based on specified cell accession_id(s). It supports creating multiple output files, each corresponding to one of the specified accession_id(s), or a single output file containing all specified accession_id, depending on the user's choice.

### Command-line Arguments:
- **--cas_json**: Path to the CAS JSON file that will be updated. This parameter is required.
- **--split_on**: Specifies the cell accession_id(s) to split the CAS file. Multiple accession_ids can be provided.
- **--multiple_outputs**: If set, create multiple output files for each term provided in split_on.
                          If not set, a single output file will be created containing all child terms,
                          and it will be named as `split_cas.json`.

### Usage Examples:
```commandline
cd src
python -m cas split_cas --cas_json path/to/cas.json --split_on term1
python -m cas split_cas --cas_json path/to/cas.json --split_on term1 term2
python -m cas split_cas --cas_json path/to/cas.json --split_on term1 term2 --multiple_outputs
```

python -m cas split_cas --cas_json path/to/cas.json --split_on term1 term2 --multiple_outputs
## Split AnnData with CAS JSON

This tool allows you to split an AnnData file based on specified CAS JSON files. It supports creating multiple output files, each corresponding to one of the specified CAS JSON files, or a single output file that contains all cells from the input CAS JSON files, depending on the user's choice.

### Command-line Arguments:
- **--anndata**: Path to the AnnData file. If not provided, AnnData will be downloaded using matrix file id from CAS JSON.
- **--cas_json**: List of CAS JSON file paths that will be used to split the AnnData file. Multiple paths can be provided.
- **--multiple_outputs**: If set, creates multiple output files for each CAS JSON file; if not set, creates a single output file containing all cells from the input CAS JSON files.
- **--compression**: Compression method utilized in anndata write function. It can be `gzip`, `lzf`,  or `None`. Default is "gzip" if flag is provided without a value. If the flag is not provided, defaults to None.

### Usage Examples:
```commandline
cd src
python -m cas split_anndata --anndata path/to/anndata.h5ad --cas_json path/to/cas.json
python -m cas split_anndata --anndata path/to/anndata.h5ad --cas_json path/to/cas.json
python -m cas split_anndata --anndata path/to/anndata.h5ad --cas_json path/to/cas_1.json path/to/cas_2.json --compression lzf
python -m cas split_anndata --anndata path/to/anndata.h5ad --cas_json path/to/cas_1.json path/to/cas_2.json 
--multiple_outputs --compression
```
