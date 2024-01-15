import pandas as pd


def get_all_annotations(
    cas: dict, show_cell_ids: bool = False, labels: list = None
) -> pd.DataFrame:
    """
    Lists all annotations.

    Args:
        cas: Cell Annotation Schema json object.
        show_cell_ids: identifies if result have 'cell_ids' column. Default value is false
        labels: list of key(labelset), value(cell_label) pairs to filter annotations
    Returns:
        Annotations data frame
    """
    df = pd.json_normalize(cas["annotations"])
    if not show_cell_ids:
        df = df.drop("cell_ids", axis=1, errors="ignore")

    if labels:
        mask = df.apply(lambda row: filter_by_label(row, labels), axis=1)
        return df[mask]
    else:
        return df


def filter_by_label(row, labels):
    """
    Filters and shows row if it
    """
    filter_row = False
    for lbl in labels:
        if row["labelset"] == lbl[0] and row["cell_label"] == lbl[1]:
            filter_row = True
    return filter_row
