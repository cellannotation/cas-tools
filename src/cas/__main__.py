import argparse
import pathlib
import sys
import warnings

from cas.anndata_conversion import merge
from cas.flatten_data_to_anndata import flatten
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
    create_spreadsheet2cas_operation_parser(subparsers)
    create_populate_cells_operation_parser(subparsers)
    create_schema_validation_operation_parser(subparsers)

    args = parser.parse_args()

    if args.action == "merge":
        json_file_path = args.json
        anndata_file_path = args.anndata
        output_file_path = args.output
        validate = args.validate

        if anndata_file_path == output_file_path:
            raise ValueError("--anndata and --output cannot be the same")

        merge(json_file_path, anndata_file_path, validate, output_file_path)
    elif args.action == "flatten":
        args = parser.parse_args()
        json_file_path = args.json
        anndata_file_path = args.anndata
        output_file_path = args.output
        validate = args.validate

        if anndata_file_path == output_file_path:
            raise ValueError("--anndata and --output cannot be the same")

        flatten(json_file_path, anndata_file_path, validate, output_file_path)
    elif args.action == "spreadsheet2cas":
        args = parser.parse_args()
        spreadsheet_file_path = args.spreadsheet
        sheet_name = args.sheet
        anndata_file_path = args.anndata
        labelsets = args.labelsets
        output_file_path = args.output

        spreadsheet2cas(
            spreadsheet_file_path, sheet_name, anndata_file_path, labelsets, output_file_path
        )
    elif args.action == "populate_cells":
        args = parser.parse_args()
        json_file_path = args.json
        anndata_file_path = args.anndata
        labelsets = None
        if "labelsets" in args and args.labelsets:
            labelsets = [item.strip() for item in str(args.labelsets).split(",")]
        populate_cell_ids(json_file_path, anndata_file_path, labelsets)
    elif args.action == "validate":
        args = parser.parse_args()
        schema = args.schema
        data_file_path = args.data
        schema_validate(schema, data_file_path)


def create_merge_operation_parser(subparsers):
    """
    Command-line Arguments:
    -----------------------
    --json      : Path to the CAS JSON schema file.
    --anndata   : Path to the AnnData file. Ideally, the location will be specified by a resolvable path in the CAS file.
    --validate  : Perform validation checks before writing to the output AnnData file.
    --output    : Output AnnData file name (default: output.h5ad).

    Usage Example:
    --------------
    cd src
    python -m cas merge --json path/to/CAS_schema.json --anndata path/to/input_anndata.h5ad --validate --output path/to/output.h5ad
    """
    parser_merge = subparsers.add_parser(
        "merge",
        description="The CAS and AnnData merge parser",
        help="Test if CAS can be merged to the AnnData and merges if possible.",
    )
    parser_merge.add_argument(
        "--json", required=True, help="Path to the CAS JSON schema file."
    )
    parser_merge.add_argument(
        "--anndata", required=True, help="Path to the AnnData file."
    )
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
    --json      : Path to the CAS JSON schema file.
    --anndata   : Path to the AnnData file. Ideally, the location will be specified by a resolvable path in the CAS file.
    --validate  : Perform validation checks before flattening to AnnData file.
    --output    : Output AnnData file name (default: output.h5ad).

    Usage Example:
    --------------
    cd src
    python -m cas flatten --json path/to/json_file.json --anndata path/to/anndata_file.h5ad --output path/to/output_file.h5ad
    """
    parser_flatten = subparsers.add_parser(
        "flatten",
        description="Flattens all content of CAS annotations to an AnnData file.",
        help="Flattens all content of CAS annotations to obs key:value pairs. "
        "Flattens all other content to key_value pairs in uns.",
    )

    parser_flatten.add_argument("--json", required=True, help="Input JSON file path")
    parser_flatten.add_argument(
        "--anndata", required=True, help="Input AnnData file path"
    )
    parser_flatten.add_argument(
        "-v",
        "--validate",
        action="store_true",
        help="Perform validation checks before writing to the output AnnData file.",
    )
    parser_flatten.add_argument(
        "--output",
        help="Output AnnData file name (default: output.h5ad)",
        default="output.h5ad",
    )
    parser_flatten.set_defaults(validate=False)


def create_spreadsheet2cas_operation_parser(subparsers):
    """
    Command-line Arguments:
    -----------------------
    --spreadsheet   : Path to the spreadsheet file.
    --sheet         : Target sheet name in the spreadsheet.
    --anndata       : Path to the AnnData file. If not provided anndata will be downloaded using CxG LINK in
                    spreadsheet.
    --labelsets     : List to determine the rank of labelsets in spreadsheet. If not provided ranks will be
                    determined using order of CELL LABELSET NAME.
    --output        : Output CAS file name (default: output.json).

    Usage Example:
    --------------
    cd src
    python -m cas spreadsheet2cas --spreadsheet path/to/spreadsheet.xlsx --sheet
    sheet_name --labelsets item1 item2 item3 --output path/to/output_file.json
    """
    parser_spreadsheet2cas = subparsers.add_parser(
        "spreadsheet2cas",
        description="Converts a spreadsheet to CAS JSON.",
        help="Converts a spreadsheet to Cell Annotation Schema (CAS) JSON.",
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
        nargs='+',
        help="List to determine the rank of labelsets in spreadsheet. If not provided ranks will be determined using "
             "order of CELL LABELSET NAME.",
    )
    parser_spreadsheet2cas.add_argument(
        "--output",
        default="output.json",
        help="Output CAS file name (default: output.json).",
    )
    parser_spreadsheet2cas.set_defaults(validate=False)


def create_populate_cells_operation_parser(subparsers):
    """
    Command-line Arguments:
    -----------------------
    --json      : Path to the CAS JSON schema file.
    --anndata   : Path to the AnnData file. Ideally, the location will be specified by a resolvable path in the CAS file.
    --labelsets : List of labelsets to update with IDs from AnnData. If value is not provided, rank '0' labelset is used.

    Usage Example:
    --------------
    cd src
    python -m cas populate_cells --json path/to/json_file.json --anndata path/to/anndata_file.h5ad --labelsets Cluster,Supercluster
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
        help="List of labelsets to update with IDs from AnnData",
        default="",
    )
    parser_populate.set_defaults(validate=False)


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


if __name__ == "__main__":
    main()
