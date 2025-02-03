### Notes on Tutorial:

- **WMB-demo.ipynb**: Ugur's interpretation of this ticket: [https://github.com/cellannotation/cas-tools/issues/143](https://github.com/cellannotation/cas-tools/issues/143)
- **CAS_Tutorial.ipynb**: New tutorial doc with extensive markdown - follows starting points structure detailed below.

---

### Work So Far:
1. **File Issues**
   1. **WMB CB Glut** - [http://cellular-semantics.cog.sanger.ac.uk/public/merged_CS20230722_CLAS_29.h5ad](http://cellular-semantics.cog.sanger.ac.uk/public/merged_CS20230722_CLAS_29.h5ad) issues:
      1. Doesn't have Neurotransmitter labelset.
      2. Limited metadata (no ontology IDs for anything except tissue).
      3. Only ENSEMBL Gene IDs in `var` (non-standard field name) --> better to modify to CxG standard?
      4. Should be updated from the latest CAS (which needs cluster metadata).
      5. UMAP embedding name doesn't work with Scanpy default standard (`X-UMAP` --> `umap`).

2. **CAS Functionality Needed**
   1. `cas flatten` should use CAS in the header by default.
   2. `cas flatten` should store CAS in the header.
   3. `flatten` should be renamed to `export2CAP` and should only flatten CAP fields.
   4. It should be possible to update with `obs`<->`cas` without going through a flattening step.

3. **Starting Points**
   1. **Introduction**
      - What is CAS / CAS-tools? Show pre-rolled `.h5ad` file in use.
      - Reporting.
      - Add a column (author annotation): [Add Author Annotations to CAS JSON](https://github.com/cellannotation/cas-tools/blob/main/docs/cli.md#add-author-annotations-to-cas-json).
      - Export to CAP: Needs introductory text to explain - `!cas flatten` (currently needs a separate file).
   2. **CxG Dataset**
      - Generate CAS from AnnData (e.g., from CxG): [Convert AnnData to CAS](https://github.com/cellannotation/cas-tools/blob/main/docs/cli.md#convert-anndata-to-cas) --> embedded `.h5ad` file.
   3. **Allen Style Taxonomy Spreadsheet + h5ad**
      - Point users to the TDT repository.
      - Steps (using WMB example):
         1. Generate CAS.
         2. Validate CAS against `.h5ad`.
         3. Populate IDs.
         4. Merge.
   4. **CAS Repository**
      - Point to TDT documentation for this.
   5. **Round-Tripping**
      - See existing notebook: `CAS-CAP round-tripping demo.ipynb` - Needs better annotation.

4. **To Demo**
   1. Export to CAP.
