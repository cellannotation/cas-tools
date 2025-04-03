import argparse
import os
import pathlib
import warnings

from cas.abc_cas_converter import abc2cas, cas2abc
from cas.add_author_annotations import add_author_annotations_from_file
from cas.anndata_conversion import merge
from cas.anndata_splitter import split_anndata_to_file
from cas.anndata_to_cas import anndata2cas
from cas.cas_splitter import split_cas_to_file
from cas.cas_to_rdf import export_to_rdf
from cas.flatten_data_to_anndata import export2cap, unflatten
from cas.populate_cell_ids import populate_cell_ids
from cas.spreadsheet_to_cas import spreadsheet2cas
from cas.validate import validate as schema_validate

warnings.filterwarnings("always")


def main():
    parser = argparse.ArgumentParser(
        prog="cas", description="Cell Type Annotation Tools cli interface."
    )
    subparsers = parser.add_subparsers(help="Available ctat actions", dest="action")

    create_merge_operation_parser(subparsers)
    create_flatten_operation_parser(subparsers)
    create_unflatten_operation_parser(subparsers)
    create_spreadsheet2cas_operation_parser(subparsers)
    create_anndata2cas_operation_parser(subparsers)
    create_abc2cas_operation_parser(subparsers)
    create_cas2abc_operation_parser(subparsers)
    create_populate_cells_operation_parser(subparsers)
    create_schema_validation_operation_parser(subparsers)
    create_cas2rdf_operation_parser(subparsers)
    create_add_author_annotations_parser(subparsers)
    create_split_cas_parser(subparsers)
    create_split_anndata_parser(subparsers)

    args = parser.parse_args()

    if args.action == "merge":
        json_file_path = args.json
        anndata_file_path = args.anndata
        output_file_path = args.output
        validate = args.validate

        if anndata_file_path == output_file_path:
            raise ValueError("--anndata and --output cannot be the same")

        merge(json_file_path, anndata_file_path, validate, output_file_path)
    elif args.action == "export2cap":
        args = parser.parse_args()
        json_file_path = args.json
        anndata_file_path = args.anndata
        output_file_path = args.output
        fill_na = args.fill_na

        if (
            anndata_file_path
            and output_file_path
            and os.path.abspath(anndata_file_path) == os.path.abspath(output_file_path)
        ):
            raise ValueError("--anndata and --output cannot be the same")
        export2cap(json_file_path, anndata_file_path, output_file_path, fill_na)
    elif args.action == "unflatten":
        args = parser.parse_args()
        json_file_path = args.json
        anndata_file_path = args.anndata
        output_anndata_path = args.output_anndata
        output_json_path = args.output_json

        if output_anndata_path and os.path.abspath(
            anndata_file_path
        ) == os.path.abspath(output_anndata_path):
            raise ValueError("--anndata and --output_anndata cannot be the same")
        if (
            json_file_path
            and output_json_path
            and os.path.abspath(json_file_path) == os.path.abspath(output_json_path)
        ):
            raise ValueError("--json and --output_json cannot be the same")
        unflatten(
            json_file_path, anndata_file_path, output_anndata_path, output_json_path
        )
    elif args.action == "spreadsheet2cas":
        args = parser.parse_args()
        spreadsheet_file_path = args.spreadsheet
        sheet_name = args.sheet
        anndata_file_path = args.anndata
        labelsets = args.labelsets
        schema_name = args.schema
        output_file_path = args.output

        spreadsheet2cas(
            spreadsheet_file_path,
            sheet_name,
            anndata_file_path,
            labelsets,
            schema_name,
            output_file_path,
        )
    elif args.action == "anndata2cas":
        args = parser.parse_args()
        anndata_file_path = args.anndata
        labelsets = args.labelsets
        output_file_path = args.output
        include_hierarchy = args.hierarchy

        anndata2cas(anndata_file_path, labelsets, output_file_path, include_hierarchy)
    elif args.action == "abc2cas":
        args = parser.parse_args()
        cat_set_file_path = args.catset
        cat_file_path = args.cat
        output_file_path = args.output

        abc2cas(cat_set_file_path, cat_file_path, output_file_path)
    # elif args.action == "cas2abc":
    #     args = parser.parse_args()
    #     json_file_path = args.json
    #     cat_set_file_path = args.catset
    #     cat_file_path = args.cat
    #
    #     cas2abc(json_file_path, cat_set_file_path, cat_file_path)
    elif args.action == "populate_cells":
        args = parser.parse_args()
        json_file_path = args.json
        anndata_file_path = args.anndata
        labelsets = args.labelsets
        validate = getattr(args, "validate", False)
        populate_cell_ids(json_file_path, anndata_file_path, labelsets, validate)
    elif args.action == "validate":
        args = parser.parse_args()
        schema = args.schema
        data_file_path = args.data
        schema_validate(schema, data_file_path)
    elif args.action == "cas2rdf":
        args = parser.parse_args()
        export_to_rdf(
            cas_schema=args.schema,
            data=args.data,
            ontology_namespace=args.ontology_ns,
            ontology_iri=args.ontology_iri,
            output_path=args.out,
            validate=not args.skip_validate,
            include_cells=not args.exclude_cells,
        )
    elif args.action == "add_author_annotations":
        args = parser.parse_args()
        # Determine the column(s) to use for joining based on the specified arguments
        if args.join_on:
            join_column = args.join_on
        elif args.join_on_cellset_ids:
            join_column = "cell_set_accession"
        elif args.join_on_labelset_label:
            join_column = ["labelset", "cell_label"]
        else:
            raise ValueError("No valid join column specified.")

        add_author_annotations_from_file(
            args.cas_json, args.csv, join_column, args.columns, args.output
        )
    elif args.action == "split_cas":
        cas_json_path = args.cas_json
        split_terms = args.split_on
        multiple_outputs = args.multiple_outputs

        split_cas_to_file(cas_json_path, split_terms, multiple_outputs)
    elif args.action == "split_anndata":
        anndata_file_path = args.anndata
        cas_json_path_list = args.cas_json
        multiple_outputs = args.multiple_outputs
        compression = args.compression

        split_anndata_to_file(
            anndata_file_path, cas_json_path_list, multiple_outputs, compression
        )


def create_merge_operation_parser(subparsers):
    """
    Command-line Arguments:
    -----------------------
    --json      : Path to the CAS JSON schema file.
    --anndata   : Path to the AnnData file. If not provided, AnnData will be downloaded using matrix file id from CAS JSON.
    --validate  : Perform validation checks before writing to the output AnnData file.
    --output    : Output AnnData file name (default: output.h5ad).

    If `--validate` is used, the following checks will be performed:
        1. Verifies that all cell barcodes (cell IDs) in CAS exist in AnnData and vice versa.
        2. Identifies matching labelset names between CAS and AnnData.
        3. Validates that cell sets associated with each annotation match between CAS and AnnData.
        4. Checks if the cell labels are identical; if not, provides options to update or terminate.

    Usage Example:
    --------------
    cd src
    python -m cas merge --json path/to/CAS_schema.json --anndata path/to/input_anndata.h5ad --validate --output path/to/output.h5ad
    """
    parser_merge = subparsers.add_parser(
        "merge",
        description="The CAS and AnnData merge operation",
        help="Test if CAS can be merged to the AnnData and merges if possible.",
    )
    parser_merge.add_argument(
        "--json", required=True, help="Path to the CAS JSON schema file."
    )
    parser_merge.add_argument("--anndata", help="Path to the AnnData file.")
    # TODO find a better argument name and Help message.
    parser_merge.add_argument(
        "-v",
        "--validate",
        action="store_true",
        help="Perform validation checks before writing to the output AnnData file.",
    )
    parser_merge.add_argument(
        "--output",
        default="output.h5ad",
        help="Output AnnData file name (default: output.h5ad).",
    )
    parser_merge.set_defaults(validate=False)


def create_flatten_operation_parser(subparsers):
    """
    Command-line Arguments:
    -----------------------
    --json      : Optional path to the CAS JSON schema file. If not provided, the CAS JSON will be extracted from the
                  AnnData file's 'uns' section.
    --anndata   : Optional path to the AnnData file. If not provided, the AnnData file will be downloaded using the
                  matrix file id from the CAS JSON.
    --output    : Optional output AnnData file name. If provided, a new flattened AnnData file will be created;
                  otherwise, the input AnnData file will be updated with the flattened data.
    --fill-na   : Optional boolean flag indicating whether to fill missing values in the 'obs' field with pd.NA. If
                  provided, missing values will be replaced with pd.NA; if not provided, they will remain as empty strings.

    Note:
        Either --json or --anndata must be supplied to execute the operation.

    Usage Example:
    --------------
    cd src
    python -m cas export2cap --json path/to/json_file.json --anndata path/to/anndata_file.h5ad --output path/to/output_file.h5ad
    """
    parser_flatten = subparsers.add_parser(
        "export2cap",
        description="Flattens all content of CAS annotations to an AnnData file.",
        help="Flattens all content of CAS annotations to obs key:value pairs.",
    )

    parser_flatten.add_argument(
        "--json",
        required=False,
        help="Optional input JSON file path. If not provided, the CAS JSON will be extracted from the AnnData file's 'uns' section.",
    )
    parser_flatten.add_argument(
        "--anndata",
        required=False,
        help="Optional input AnnData file path. If not provided, the AnnData file will be downloaded using the matrix file id from the CAS JSON.",
    )
    parser_flatten.add_argument(
        "--output",
        required=False,
        help="Output AnnData file name.",
    )
    parser_flatten.add_argument(
        "--fill-na",
        required=False,
        action="store_true",
        help="Boolean flag indicating whether to fill missing values in the 'obs' field with pd.NA. If provided, "
        "missing values will be replaced with pd.NA; if not provided, they will remain as empty strings.",
    )


def create_unflatten_operation_parser(subparsers):
    """
        Command-line Arguments:
    -----------------------
    --anndata           : Path to the AnnData file.
    --json              : Optional path to the CAS JSON file. If provided, the 'annotations'
    within the file will be updated. If not provided, a new CAS JSON file will be created.
    --output_anndata    : Optional output AnnData file name. If not provided, 'unflattened.h5ad' will be used as
    default name.
    --output_json       : Optional output CAS JSON file name. If not provided, 'cas.json' will be used as default name.

    Usage Example:
    --------------
    cd src
    python -m cas unflatten --anndata path/to/anndata_file.h5ad --json path/to/json_file.json --output_anndata
    path/to/output_file.h5ad --output_json path/to/output_cas.json
    """
    parser_unflatten = subparsers.add_parser(
        "unflatten",
        description="Unflattens all content of a flattened AnnData file to a CAS JSON file. Also creates an "
        "unflattened AnnData file.",
        help="Converts flattened AnnData content back to unflattened version and creates CAS JSON files.",
    )

    parser_unflatten.add_argument(
        "--anndata",
        required=True,
        help="Path to the input AnnData file that contains flattened data.",
    )
    parser_unflatten.add_argument(
        "--json",
        required=False,
        help="Optional path to the CAS JSON file. If provided, the 'annotations' within the file "
        "will be updated. If not provided, a new CAS JSON file will be created.",
    )
    parser_unflatten.add_argument(
        "--output_anndata",
        required=False,
        help="Optional output AnnData file name. If not provided, 'unflattened.h5ad' will be used as default name.",
    )
    parser_unflatten.add_argument(
        "--output_json",
        required=False,
        default="cas.json",
        help="Optional output CAS JSON file name. If not provided, 'cas.json' will be used as the default file name.",
    )


def create_spreadsheet2cas_operation_parser(subparsers):
    """
    Command-line Arguments:
    -----------------------
    --spreadsheet   : Path to the spreadsheet file.
    --sheet         : Target sheet name in the spreadsheet.
    --anndata       : Path to the AnnData file. If not provided, AnnData will be downloaded using CxG LINK in
                    spreadsheet.
    --labelsets     : List of names of observation (obs) fields used to record author cell type names,
        which determine the rank of labelsets in a spreadsheet. If not provided, ranks will be determined based on
        the order of the fields specified in the CELL LABELSET NAME column.
    --schema        : Name of cell annotation schema used to in spreadsheet. It can be one of 'base', 'bican' or 'cap'.
    --output        : Output CAS file name (default: output.json).

    Usage Example:
    --------------
    cd src
    python -m cas spreadsheet2cas --spreadsheet path/to/spreadsheet.xlsx --sheet sheet_name --labelsets item1 item2 item3 --output path/to/output_file.json
    """

    parser_spreadsheet2cas = subparsers.add_parser(
        "spreadsheet2cas",
        description="Converts a spreadsheet to CAS JSON.",
        help="Converts a spreadsheet to Cell Annotation Schema (CAS) JSON.",
        usage="cas spreadsheet2cas --spreadsheet path/to/spreadsheet.xlsx --sheet sheet_name --labelsets item1 item2 item3 --output path/to/output_file.json",
    )

    parser_spreadsheet2cas.add_argument(
        "--spreadsheet", required=True, help="Path to the spreadsheet file."
    )
    parser_spreadsheet2cas.add_argument(
        "--sheet", required=False, help="Target sheet name in the spreadsheet."
    )
    parser_spreadsheet2cas.add_argument(
        "--anndata",
        default=None,
        help="Path to the AnnData file. If not provided, AnnData will be downloaded using CxG LINK in spreadsheet.",
    )
    parser_spreadsheet2cas.add_argument(
        "--labelsets",
        default=None,
        nargs="+",
        help="List of names of observation (obs) fields used to record author cell type names, which determine the "
        "rank of labelsets in a spreadsheet. If not provided, ranks will be determined based on the order of the"
        " fields specified in the CELL LABELSET NAME column.",
    )
    parser_spreadsheet2cas.add_argument(
        "--schema",
        default="cap",
        help="Name of cell annotation schema used to in spreadsheet. It can be one of 'base', 'bican' or 'cap'.",
        choices=["base", "bican", "cap"],
    )
    parser_spreadsheet2cas.add_argument(
        "--output",
        default="output.json",
        help="Output CAS file name (default: output.json).",
    )


def create_anndata2cas_operation_parser(subparsers):
    """
    Command-line Arguments:
    -----------------------
    --anndata       : Path to the AnnData file.
    --labelsets     : List of labelsets, which are names of observation (obs) fields used to record author cell type
    names. The labelsets should be provided in order, starting from rank 0 (leaf nodes) and ascending to higher ranks.

    --output        : Output CAS file name (default: output.json).
    --hierarchy     : Flag indicating whether to include hierarchy in the output.


    Usage Example:
    --------------
    cd src
    python -m cas anndata2cas --anndata path/to/anndata.h5ad --labelsets item1 item2 item3 --output
    path/to/output_file.json
    """
    parser_anndata2cas = subparsers.add_parser(
        "anndata2cas",
        description="Converts an anndata to CAS JSON.",
        help="Converts an anndata to Cell Annotation Schema (CAS) JSON.",
        usage="cas anndata2cas --anndata path/to/anndata.h5ad --labelsets item1 item2 item3 --output path/to/output_file.json",
    )

    parser_anndata2cas.add_argument(
        "--anndata",
        required=True,
        help="Path to the AnnData file.",
    )
    parser_anndata2cas.add_argument(
        "--labelsets",
        required=True,
        nargs="+",
        help="List of labelsets, which are names of observation (obs) fields used to record author cell type names. "
        "The labelsets should be provided in order, starting from rank 0 (leaf nodes) and ascending to higher ranks.",
    )
    parser_anndata2cas.add_argument(
        "--output",
        default="output.json",
        help="Output CAS file name (default: output.json).",
    )
    parser_anndata2cas.add_argument(
        "--hierarchy",
        action="store_true",
        help="Include hierarchy in the output.",
    )


def create_abc2cas_operation_parser(subparsers):
    """
    Command-line Arguments:
    -----------------------
    --catset       : Path to the Cluster Annotation Term Set file.
    --cat           : Path to the Cluster Annotation Term file.
    --output        : Output CAS file name (default: output.json).


    Usage Example:
    --------------
    cd src

    python -m cas abc2cas --catset path/to/cluster_annotation_term_set.csv --cat path/to/cluster_annotation_term.csv
    --output path/to/output_file.json
    """
    parser_abc2cas = subparsers.add_parser(
        "abc2cas",
        description="Converts given ABC cluster_annotation files to CAS JSON.",
        help="Converts given ABC cluster_annotation files to Cell Annotation Schema (CAS) JSON.",
    )

    parser_abc2cas.add_argument(
        "--catset",
        required=True,
        help="Path to the Cluster Annotation Term Set file.",
    )
    parser_abc2cas.add_argument(
        "--cat",
        required=True,
        help="Path to the Cluster Annotation Term file.",
    )
    parser_abc2cas.add_argument(
        "--output",
        default="output.json",
        help="Output CAS file name (default: output.json).",
    )


def create_cas2abc_operation_parser(subparsers):
    """
    Command-line Arguments:
    -----------------------
    --json         : Path to the CAS JSON schema file.
    --catset       : Path to the Cluster Annotation Term Set file.
    --cat          : Path to the Cluster Annotation Term file.


    Usage Example:
    --------------
    cd src

    python -m cas cas2abc --json path/to/json_file.json --catset path/to/cluster_annotation_term_set.csv --cat
    path/to/cluster_annotation_term.csv
    """
    parser_cas2abc = subparsers.add_parser(
        "cas2abc",
        description="Converts given CAS JSON to ABC files. NOTE: This command is under development and may change.",
        help="Converts given Cell Annotation Schema (CAS) to ABC files: cluster_annotation_term and "
        "cluster_annotation_term_set, and writes them to files with cat_file_path and cat_set_file_path. (Under Development)",
    )

    parser_cas2abc.add_argument(
        "--json", required=True, help="Input CAS JSON file path"
    )
    parser_cas2abc.add_argument(
        "--catset",
        required=True,
        help="Path to the Cluster Annotation Term Set file.",
    )
    parser_cas2abc.add_argument(
        "--cat",
        required=True,
        help="Path to the Cluster Annotation Term file.",
    )


def create_populate_cells_operation_parser(subparsers):
    """
    Creates a command-line argument parser for the `populate_cells` operation.

    This command populates CellIDs in a CAS JSON file using data from an AnnData file. It ensures
    alignment between labelset values in AnnData's `obs` DataFrame and CAS-defined labelsets before
    updating cell IDs. If the `--validate` flag is enabled, validation errors will cause the
    program to exit.

    Command-line Arguments:
    -----------------------
    --json      (required) : Path to the CAS JSON schema file.
    --anndata   (required) : Path to the AnnData file. Ideally, this should be specified by a
                             resolvable path in the CAS file.
    --labelsets (optional) : A space-separated list of labelsets to update with IDs from AnnData.
                             If not provided, the labelset with rank '0' is used by default. The
                             labelsets should be provided in hierarchical order,
                             starting from rank 0 (leaf nodes) and ascending to higher ranks.
    --validate  (optional) : If set, the program exits with an error if validation fails.
                             Otherwise, it logs warnings and continues execution.

    Usage Example:
    --------------
    cd src
    python -m cas populate_cells --json path/to/json_file.json --anndata path/to/anndata_file.h5ad --labelsets Cluster Supercluster --validate

    This command:
    - Reads the CAS JSON schema from `path/to/json_file.json`.
    - Reads the AnnData file from `path/to/anndata_file.h5ad`.
    - Aligns the labelsets **Cluster** and **Supercluster** with the AnnData data.
    - Runs validation checks.
    - If `--validate` is enabled, the process **stops on validation failure**.

    Returns:
        None
    """
    parser_populate = subparsers.add_parser(
        "populate_cells",
        description="The CAS cell IDs population operation.",
        help="Checks for alignment between obs key:value pairs in AnnData file and labelset:cell_label pairs in CAS for some specified list of labelsets. If they are aligned, updates cell_ids in CAS.",
    )

    parser_populate.add_argument("--json", required=True, help="Input JSON file path")
    parser_populate.add_argument(
        "--anndata", required=True, help="Input AnnData file path"
    )
    parser_populate.add_argument(
        "--labelsets",
        default=None,
        nargs="+",
        help="Space-separated list of labelsets to update with IDs from AnnData",
    )
    parser_populate.add_argument(
        "--validate",
        action="store_true",
        help="Enable strict validation. If validation fails, the process will exit with an error.",
    )



def create_schema_validation_operation_parser(subparsers):
    """
    Command-line Arguments:
    -----------------------
    --schema    : One of 'base', 'bican' or 'cap'. Identifies the CAS schema to validate data against.
    --data   : Path to the data file (or folder) to validate

    Usage Example:
    --------------
    cd src
    python -m cas validate --schema bican --data path/to/file
    """
    parser_validate = subparsers.add_parser(
        "validate",
        description="The CAS file structure validator",
        help="Test if given CAS is compatible with the CAS schema.",
    )

    parser_validate.add_argument(
        "--schema", required=True, help="Schema name: one of 'base', 'bican' or 'cap'"
    )
    parser_validate.add_argument(
        "--data",
        required=True,
        help="Path to the data file (or folder) to validate",
        type=pathlib.Path,
    )


def create_cas2rdf_operation_parser(subparsers):
    """
    Command-line Arguments:
    -----------------------
    --schema    : (Optional) Name of the CAS release (such as one of `base`, `cap`, `bican`) or path to the
                    CAS schema file or url of the schema file. If not provided, reads the `base` CAS schema from the cas module.
    --data   : Path to the json data file
    --ontology_ns    : Ontology namespace (e.g. `MTG`)
    --ontology_iri    : Ontology IRI (e.g. `https://purl.brain-bican.org/ontology/AIT_MTG/`)
    --out    : The output RDF file path.
    --skip_validate    : (Optional) Determines if data-schema validation checks will be performed. Validations are performed by default.
    --exclude_cells    : (Optional) Determines if cell data will be included in the RDF output. Cell data is exported to RDF by default.

    Usage Example:
    --------------
    cd src
    python -m cas cas2rdf --schema bican --data path/to/file.json --ontology_ns MTG --ontology_iri https://purl.brain-bican.org/ontology/AIT_MTG/ --out path/to/output.rdf --exclude_cells
    """
    parser_cas2rdf = subparsers.add_parser(
        "cas2rdf",
        description="CAS to RDF convertor.",
        help="Converts given CAS data into RDF format.",
        usage="cas cas2rdf --schema bican --data path/to/file.json --ontology_ns MTG --ontology_iri https://purl.brain-bican.org/ontology/AIT_MTG/ --labelsets Cluster Subclass Class --out path/to/output.rdf --exclude_cells",
    )

    parser_cas2rdf.add_argument(
        "-s",
        "--schema",
        help="Name of the CAS release (such as one of `base`, `cap`, `bican`) or path to the CAS schema file or url of the schema file. If not provided, reads the `base` CAS schema from the cas module.",
        default="base",
    )
    parser_cas2rdf.add_argument(
        "-d",
        "--data",
        required=True,
        help="Path to the json data file",
        type=pathlib.Path,
    )
    parser_cas2rdf.add_argument(
        "-ns",
        "--ontology_ns",
        required=True,
        help="Ontology namespace (e.g. `MTG`)",
    )
    parser_cas2rdf.add_argument(
        "-iri",
        "--ontology_iri",
        required=True,
        help="Ontology IRI (e.g. `https://purl.brain-bican.org/ontology/AIT_MTG/`)",
    )
    parser_cas2rdf.add_argument(
        "-o",
        "--out",
        required=True,
        help="The output RDF file path.",
        type=pathlib.Path,
    )
    parser_cas2rdf.add_argument(
        "-sv",
        "--skip_validate",
        action="store_true",
        help="Determines if data-schema validation checks will be performed. Validations are performed by default.",
    )
    parser_cas2rdf.add_argument(
        "-ec",
        "--exclude_cells",
        action="store_true",
        help="Determines if cell data will be included in the RDF output. Cell data is exported to RDF by default.",
    )


def create_add_author_annotations_parser(subparsers):
    """
    Adds a command line parser for the operation to add author annotation fields to CAS JSON files using data from
    a specified CSV file. This operation processes the input CSV and CAS JSON files based on the specified matching
    columns, optionally using selected columns for annotations, and outputs the annotated CAS JSON to a specified file.
    If specific columns are not provided, all columns from the CSV file will be used.

    Command-line Arguments:
    -----------------------
    --cas_json    : Path to the CAS JSON file that will be updated with annotations.
    --csv         : Specifies the path to the CSV file.
    --join_on     : Specifies the single column name in the CSV used for matching records. Each row must have a unique value in this column.
    --join_on_cell_set_id     : Use 'cell_set_id' as the column for matching records.
    --join_on_labelset_label : Use a pair of 'labelset', 'cell_label' columns for matching records.
    --columns     : Optionally specifies which columns in the CSV will be used for annotations. If not provided,
    all columns will be used. Column names containing spaces must be enclosed in quotes (e.g., "Columns Name").
    --output      : Specifies the output file name for the annotated CAS JSON (default: output.json).

    Usage Example:
    --------------
    cd src
    python -m cas add_author_annotations --cas_json path_to_cas.json --csv path_to_csv --join_on CrossArea_cluster --columns random_annotation_x random_annotation_y --output annotated_output.json
    python -m cas add_author_annotations --cas_json path_to_cas.json --csv path_to_csv --join_on_cell_set_id --output annotated_output.json
    python -m cas add_author_annotations --cas_json path_to_cas.json --csv path_to_csv --join_on_labelset_label --output annotated_output.json
    """
    parser_add_annotations = subparsers.add_parser(
        "add_author_annotations",
        description="Add author annotation fields to a CAS JSON file from specified CSV columns. If no columns are specified, all columns are used.",
        help="Adds annotation fields from CSV to a specified CAS JSON file, using all columns if none are specified.",
    )
    parser_add_annotations.add_argument(
        "--cas_json", required=True, help="Path to the CAS JSON file to be annotated."
    )
    parser_add_annotations.add_argument(
        "--csv", required=True, help="Path to the CSV file containing annotation data."
    )
    parser_add_annotations.add_argument(
        "--join_on",
        help="Specifies the single column name in the CSV used for matching records. Each row must have a unique value in this column.",
    )
    parser_add_annotations.add_argument(
        "--join_on_cellset_ids",
        action="store_true",
        help="Use 'cell_ids' as the column for matching CAS records.",
    )
    parser_add_annotations.add_argument(
        "--join_on_labelset_label",
        action="store_true",
        help="Use a pair of 'labelset', 'cell_label' columns for matching CAS records.",
    )
    parser_add_annotations.add_argument(
        "--columns",
        nargs="+",
        help="Optional space-separated list of column names in the CSV to be added as annotations. All columns are "
        "used if none are specified. Column names containing spaces must be enclosed in quotes (e.g., 'Columns Name').",
    )
    parser_add_annotations.add_argument(
        "--output",
        default="output.json",
        help="Output CAS file name (default: output.json).",
    )


def create_split_cas_parser(subparsers):
    """
    Command-line Arguments:
    -----------------------
    --cas_json          : Path to the CAS JSON file that will be split.
    --split_on          : Cell accession_id(s) to split the CAS file.
    --multiple_outputs  : If set, create multiple output files for each term provided in split_on.
                          If not set, a single output file will be created containing all child terms,
                          and it will be named as `split_cas.json`.

    Usage Example:
    --------------
    cd src
    python -m cas split_cas --cas_json path/to/cas.json --split_on term1
    python -m cas split_cas --cas_json path/to/cas.json --split_on term1 term2
    python -m cas split_cas --cas_json path/to/cas.json --split_on term1 term2 --multiple_outputs
    """
    parser_split_cas = subparsers.add_parser(
        "split_cas",
        description="Split CAS JSON file based on specified cell label/s.",
        help="Split a CAS JSON file into multiple files based on one or more cell accession_id(s).",
    )
    parser_split_cas.add_argument(
        "--cas_json", required=True, help="Path to the CAS JSON file that will be split"
    )
    parser_split_cas.add_argument(
        "--split_on", nargs="+", help="Cell accession_id(s) to split the CAS file."
    )
    parser_split_cas.add_argument(
        "--multiple_outputs",
        action="store_true",
        help="If set, create multiple output files for each split_on term; if not set, create a single output file named `split_cas.json`.",
    )


def create_split_anndata_parser(subparsers):
    """
    Command-line Arguments:
    -----------------------
    --anndata           : Path to the AnnData file. If not provided, AnnData will be downloaded using matrix file id from CAS JSON.
    --cas_json_list     : List of CAS JSON file paths that will be used to split the AnnData file.
    --multiple_outputs  : If set, creates multiple output files for each term provided in split_on.
                          If not set, creates a single output file containing all cell_ids from the input CAS JSON files.
    --compression       : Compression method utilized in anndata write function. It can be `gzip`, `lzf`,
                            or `None`.Default is "gzip" if flag is provided without a value. If the flag is not provided,
                            defaults to None.
    Usage Example:
    --------------
    cd src
    python -m cas split_anndata --cas_json path/to/cas.json
    python -m cas split_anndata --anndata path/to/anndata.h5ad --cas_json path/to/cas.json
    python -m cas split_anndata --anndata path/to/anndata.h5ad --cas_json path/to/cas_1.json
    path/to/cas_2.json --compression lzf
    python -m cas split_anndata --anndata path/to/anndata.h5ad --cas_json path/to/cas_1.json
    path/to/cas_2.json
    --multiple_outputs --compression
    """
    parser_split_anndata = subparsers.add_parser(
        "split_anndata",
        description="Splits an AnnData file based on specified CAS JSON files.",
        help="Splits an AnnData file into multiple files based on one or more CAS JSON files.",
    )
    parser_split_anndata.add_argument("--anndata", help="Path to the AnnData file.")
    parser_split_anndata.add_argument(
        "--cas_json",
        required=True,
        nargs="+",
        help="List of CAS JSON file paths that will be used to split the AnnData file.",
    )
    parser_split_anndata.add_argument(
        "--multiple_outputs",
        action="store_true",
        help="If set, creates multiple output files for each CAS JSON file; if not set, creates a single output file.",
    )
    parser_split_anndata.add_argument(
        "--compression",
        nargs="?",
        default=None,
        choices=["gzip", "lzf", None],
        action=CompressionAction,
        help="Compression method utilized in anndata write function, can be `gzip`, `lzf`, or `None`.",
    )


class CompressionAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        # If --compression is specified without a value, use 'gzip'
        if values is None:
            setattr(namespace, self.dest, "gzip")
        else:
            setattr(namespace, self.dest, values)


if __name__ == "__main__":
    main()
