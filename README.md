# Single-Cell Viewer (Streamlit + Scanpy)

An interactive web app to explore single-cell RNA-seq data using [Streamlit](https://streamlit.io/) and [Scanpy](https://scanpy.readthedocs.io/).  

This app is built around a demo dataset of mouse brown adipose tissue (BAT), but also allows users to upload their own `.h5ad` AnnData files for exploration.

---

## Features
- **Try it online**: No setup needed — test the app with the included minimal dataset here:  
  [shamsi-2021-snrnaseq-visualizer.streamlit.app](https://shamsi-2021-snrnaseq-visualizer.streamlit.app/)  
- **Subset by cell type**: Choose which cell types to visualize via a sidebar multiselect.  
- **Summary metrics**: Displays mean reads/cell and mean genes/cell.  
- **Interactive visualizations**:  
  - UMAP colored by cell type  
  - Violin plots of selected gene expression per cell type  
  - Dot plots of selected gene expression by condition (if `cond` metadata available)  
- **Metadata preview**: View the first rows of cell-level metadata (`.obs`).

---

## Repository structure
```
.
├── streamlit_app.py        # Main Streamlit application
├── requirements.txt        # Python dependencies
├── adata_shrinker.ipynb    # Notebook to generate a small demo dataset
└── README.md               # This file
```

---

## Quickstart

### 1. Clone this repository
```bash
git clone https://github.com/camara-h/single-nuclei-visualizer.git
cd single-nuclei-visualizer
```

### 2. Install dependencies
We recommend using a virtual environment (e.g. `venv` or `conda`).  
Install requirements with:
```bash
pip install -r requirements.txt
```

### 3. Run the app
```bash
streamlit run streamlit_app.py
```

Open the link that appears in your terminal (usually `http://localhost:8501`).

---

## Demo and full datasets

- **Minimal dataset**:  
  A small `.h5ad` file (`shamsi_adata_demo_mini.h5ad`) is included for testing.  
  - Subsampled to ~100 cells per cell type  
  - Top 200 highly variable genes retained  
  - Lightweight enough for GitHub hosting  

- **Full dataset**:  
  The full dataset (`shamsi_adata.h5ad`) will be deposited publicly soon..  
  To use it, download the file and replace this line in `streamlit_app.py`:
  ```python
  adata = sc.read("shamsi_adata_demo_mini.h5ad")
  ```
  with:
  ```python
  adata = sc.read("shamsi_adata.h5ad")
  ```

- **Custom datasets**:  
  You can also upload your own `.h5ad` file through the app’s sidebar.  
  - Must contain `obs['cell_type_name']` for cell-type annotations  
  - Optional: `obs['cond']` for experimental conditions  
  - If `.raw` is available, it will be used for gene expression visualization

---

## Dependencies
See [`requirements.txt`](requirements.txt) for full list. Main libraries include:
- Streamlit ≥ 1.34  
- Scanpy ≥ 1.10  
- Matplotlib ≥ 3.8  
- Pandas, NumPy, AnnData  

---

## Citation
The demo dataset is derived from:  
[**Comprehensive analysis of intercellular communication in the thermogenic adipose niche. 
Shamsi F, Zheng R, Ho LL, Chen K, Tseng YH. 
Commun Biol. 2023 Jul 21;6(1):761. doi: 10.1038/s42003-023-05140-2. PMID: 37479789; PMCID: PMC10361964.**](https://doi.org/10.1038/s42003-023-05140-2)

---

## Acknowledgments
- Built with [Streamlit](https://streamlit.io/) and [Scanpy](https://scanpy.readthedocs.io/).  
- Demo dataset courtesy of the authors above.  
