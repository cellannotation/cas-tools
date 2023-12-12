# Script to generate Siletti CAS json representation.

import os

from cas.model import (CellTypeAnnotation, Annotation, Labelset, AnnotationTransfer,
                       UserAnnotation, AutomatedAnnotation)
from cas.file_utils import read_csv_to_dict, write_json_file
from cas.accession.incremental_accession_manager import IncrementalAccessionManager

CLUSTERS_TSV = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                            "../spreadsheet/Siletti_for review - all_clusters.tsv")
N_SUPER_CLUSTERS_TSV = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                    "../spreadsheet/Siletti_for review - neuron_supercluster.tsv")
NN_SUPER_CLUSTERS_TSV = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                     "../spreadsheet/Siletti_for review - non-neuronal_supercluster.tsv")
NOMENCLATURE_TABLE = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                     "../nomenclature_table_CS202210140.csv")


def main():
    cas_neuronal = CellTypeAnnotation("", list())
    cas_non_neuronal = CellTypeAnnotation("", list())

    labelsets = list()
    cluster_ls = Labelset("Cluster_name")
    cluster_ls.rank = "0"
    supercluster_ls = Labelset("supercluster_term")
    supercluster_ls.rank = "1"
    labelsets.append(cluster_ls)
    labelsets.append(supercluster_ls)
    cas_neuronal.labelsets = labelsets
    cas_non_neuronal.labelsets = labelsets

    nm_headers, nm_records = read_csv_to_dict(NOMENCLATURE_TABLE, id_column_name="cell_set_preferred_alias")
    last_accession_id = str(nm_records[list(nm_records.keys())[-1]]["cell_set_accession"]).split("_")
    accession_prefix = last_accession_id[0] + "_"
    last_accession_number = int(last_accession_id[1])

    accession_manager = IncrementalAccessionManager(accession_prefix=accession_prefix, last_accession_id=last_accession_number)

    annotations_neuronal = list()
    annotations_non_neuronal = list()

    nsc_headers, nsc_records = read_csv_to_dict(N_SUPER_CLUSTERS_TSV, delimiter="\t")
    for record_key in nsc_records:
        record = nsc_records[record_key]
        annotation = Annotation("", "")
        annotation.labelset = "supercluster_term"
        extract_annotation_data(record, record_key, annotation, nm_records, accession_manager)
        annotations_neuronal.append(annotation)

    nnsc_headers, nnsc_records = read_csv_to_dict(NN_SUPER_CLUSTERS_TSV, delimiter="\t")
    for record_key in nnsc_records:
        record = nnsc_records[record_key]
        annotation = Annotation("", "")
        annotation.labelset = "supercluster_term"
        extract_annotation_data(record, record_key, annotation, nm_records, accession_manager)
        annotations_non_neuronal.append(annotation)

    headers, records = read_csv_to_dict(CLUSTERS_TSV, id_column_name="Cluster ID", delimiter="\t")
    for record_key in records:
        record = records[record_key]
        annotation = Annotation("", "")
        annotation.labelset = "Cluster_name"

        extract_annotation_data(record, record_key, annotation, nm_records, accession_manager)

        if record.get("Supercluster", "") in nsc_records:
            annotation.parent_cell_set_name = record.get("Supercluster", "")
            annotation.parent_cell_set_accession = get_accession(annotations_neuronal, annotation.parent_cell_set_name)
            annotations_neuronal.append(annotation)
        elif record.get("Supercluster", "") in nnsc_records:
            annotation.parent_cell_set_name = record.get("Supercluster", "")
            annotation.parent_cell_set_accession = get_accession(annotations_non_neuronal, annotation.parent_cell_set_name)
            annotations_non_neuronal.append(annotation)
        else:
            raise Exception("Supercluster not found: {} at {}".format(record.get("Supercluster", ""), annotation.cell_label))

    cas_neuronal.annotations = annotations_neuronal
    cas_non_neuronal.annotations = annotations_non_neuronal

    write_json_file(cas_neuronal, "../Siletti_all_neurons.json")
    write_json_file(cas_non_neuronal, "../Siletti_all_non_neuronal_cells.json")


def extract_annotation_data(record, record_key, annotation, nm_records, accession_manager):
    if "Cluster name" in record:
        annotation.cell_label = record.get("Cluster name", "")
    elif "SUPERCLUSTER" in record:
        annotation.cell_label = record.get("SUPERCLUSTER", "")
    else:
        raise Exception("Cell Label couldn't be identified for: {}".format(record_key))
    annotation.cell_ontology_term_id = record.get("CELL_ONTOLOGY_TERM_ID", "")
    annotation.cell_ontology_term = record.get("CELL_ONTOLOGY_TERM", "")
    annotation.rationale = record.get("RATIONALE", "")
    annotation.rationale_dois = [x.strip() for x in record.get("RATIONALE_DOI", "").split(",")]
    annotation.marker_gene_evidence = [x.strip() for x in record.get("POSITIVE_GENE_EVIDENCE", "").split(",")]
    if "Cluster ID" in record:
        annotation.add_user_annotation("Cluster ID", record["Cluster ID"])

    if "Cluster name" in record:
        if record["Cluster name"] and record["Cluster name"] in nm_records:
            annotation.cell_set_accession = nm_records[record["Cluster name"]]["cell_set_accession"]
        else:
            raise Exception("AccessionId couldn't be identified for: {}".format(record_key))
    elif "SUPERCLUSTER" in record:
        if record["SUPERCLUSTER"] and record["SUPERCLUSTER"] in nm_records:
            annotation.cell_set_accession = nm_records[record["SUPERCLUSTER"]]["cell_set_accession"]
        else:
            annotation.cell_set_accession = accession_manager.generate_accession_id()


def get_accession(annotations, cell_label):
    for annotation in annotations:
        if cell_label == annotation.cell_label:
            return annotation.cell_set_accession
    return ""


if __name__ == "__main__":
    main()
