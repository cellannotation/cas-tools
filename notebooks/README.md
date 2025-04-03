# CAS-tools Tutorials and Snippets

This guide explains how to set up a virtual environment, install the necessary dependencies, and launch Jupyter Notebook to work with CAS-tools.

## 1. Setting Up Your Virtual Environment

1. **Navigate to the `notebooks` directory:**

    ```bash
    cd notebooks
    ```

2. **Create a virtual environment:**

    ```bash
    python3 -m venv venv
    ```

3. **Activate the virtual environment:**

    ```bash
    source venv/bin/activate
    ```

## 2. Installing Dependencies

### a. CAS-tools

Install the CAS-tools package:

```bash
pip install cas-tools
```

### b. IPython Kernel

Install `ipykernel` and set up a kernel for your virtual environment:

```bash
pip install ipykernel
python3 -m ipykernel install --user --name=venv
```

### c. Jupyter Notebook

Install Jupyter Notebook:

```bash
pip install notebook
```

## 3. Running Jupyter Notebook

Launch Jupyter Notebook with the following command:

```bash
jupyter notebook
```

Once started, Jupyter Notebook will open in your default web browser. You can now navigate through your notebooks and begin working with CAS-tools.

## Notebooks

* [CAS-tutorial](CAS-tutorial.ipynb):  A basic guide to the Cell Annotation chema, how it works with h5ad files and how to use CAS-tools to access it.
* [CAS_support_for_hierarchy](CAS_support_for_hierarchy.ipynb):  How CAS supports heirarchical annotation of cell type.
* [allen_spreadsheet_to_cas_to_cap]([allen_spreadsheet_to_cas_to_cap.ipynb): Annotation metadata is typically mainteined in spreadsheets, this notebook shows how to covert a spreadsheet to CAS, embed CAS in an h5ad file of annotated data and how to convert the resulting h5ad file for loading to CAP.
* [CAS-CAP_Roundtrip](CAS-CAP_Roundtrip.ipynb): How CAS-tools supports round-tripping between CAS and CAP, and how it deals
* TBA: 
* TBA: 