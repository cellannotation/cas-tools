import sys
import argparse
import pathlib

from cas.anndata_conversion import merge
from cas.flatten_data_to_anndata import flatten
from cas.populate_cell_ids import populate_cell_ids


def main():
    parser = argparse.ArgumentParser(prog="cas", description='Cell Type Annotation Tools cli interface.')
    subparsers = parser.add_subparsers(help='Available ctat actions', dest='action')

    create_merge_operation_parser(subparsers)
    create_flatten_operation_parser(subparsers)
    create_populate_cells_operation_parser(subparsers)

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
    elif args.action == "populate_cells":
        args = parser.parse_args()
        json_file_path = args.json
        anndata_file_path = args.anndata
        labelsets = None
        if "labelsets" in args and args.labelsets:
            labelsets = [item.strip() for item in str(args.labelsets).split(",")]
        populate_cell_ids(json_file_path, anndata_file_path, labelsets)


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
    parser_merge = subparsers.add_parser("merge",
                                         description="The CAS and AnnData merge parser",
                                         help="Test if CAS can be merged to the AnnData and merges if possible.")
    parser_merge.add_argument(
        "--json", required=True, help="Path to the CAS JSON schema file."
    )
    parser_merge.add_argument("--anndata", required=True, help="Path to the AnnData file.")
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
    parser_flatten = subparsers.add_parser("flatten",
                                         description="The CAS and AnnData merge parser",
                                         help="Test if CAS can be merged to the AnnData and merges if possible.")

    parser_flatten.add_argument("--json", required=True, help="Input JSON file path")
    parser_flatten.add_argument("--anndata", required=True, help="Input AnnData file path")
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
    parser_populate = subparsers.add_parser("populate_cells",
                                         description="The CAS and AnnData merge parser",
                                         help="Test if CAS can be merged to the AnnData and merges if possible.")

    parser_populate.add_argument("--json", required=True, help="Input JSON file path")
    parser_populate.add_argument("--anndata", required=True, help="Input AnnData file path")
    parser_populate.add_argument(
        "--labelsets",
        help="List of labelsets to update with IDs from AnnData",
        default="",
    )
    parser_populate.set_defaults(validate=False)


if __name__ == "__main__":
    main()
