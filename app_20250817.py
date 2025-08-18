# -------------------------------------------------------------------------------------------------
# Single-Cell Viewer (Streamlit + Scanpy) — fully annotated
#
# What this app does:
# - Loads an AnnData (.h5ad) single-cell RNA-seq dataset (mouse BAT)
# - Lets you subset cells by cell type via a sidebar multiselect
# - Provides quick summary metrics (mean reads/genes per cell)
# - Visualizes UMAP, a violin plot for a selected gene, and a dot plot by condition
# - Supports optional Google Drive download + MD5 integrity checking (via gdown)
#
# Notes for readers:
# - This annotation explains why each piece of code exists, Streamlit behaviors (e.g., session_state),
#   and Scanpy plotting choices.
# - No functional changes were made — only comments were added.
# -------------------------------------------------------------------------------------------------

import os                   # Standard library: filesystem checks (exists, remove), paths, etc.
import hashlib              # Standard library: computing MD5 checksums to verify file integrity
import streamlit as st      # Streamlit: UI framework for building the web app
import scanpy as sc         # Scanpy: single-cell analysis; provides AnnData and plotting utilities
import gdown                # gdown: convenient Google Drive downloads by URL
import matplotlib.pyplot as plt  # Matplotlib: used by Scanpy for plotting; Streamlit renders via st.pyplot

# Configure Streamlit page: title in browser tab + wide layout for more room
st.set_page_config(page_title="Single-Cell Viewer", layout="wide")

# === Settings ===
# Google Drive file link for the demo .h5ad. Replace with your own file if needed.
# (Comment "seu ID" from original: your Google Drive file ID is embedded in this link.)
DRIVE_URL = "https://drive.google.com/file/d/1HYiv3_81VGgyBHpkOGKp6QvXu9HydvQI/view?usp=sharing"  # seu ID

# Local cache path where the dataset will live. If the file exists, we reuse it (unless MD5 check fails).
LOCAL_PATH = "./data/shamsi_adata_demo_mini.h5ad"

# Optional MD5 checksum. If provided, we verify the downloaded file to guard against corruption.
# If None, we skip MD5 validation (faster, but less safe).
EXPECTED_MD5 = None  # opcional: coloque o md5 do arquivo para checagem

# === Helper functions ===
def md5sum(path, chunk=1024*1024):
    """
    Compute the MD5 checksum of a file in chunks (default: 1 MB),
    which keeps memory usage low even for large files.
    """
    h = hashlib.md5()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:         # EOF
                break
            h.update(b)
    return h.hexdigest()

def ensure_local_file(url, path, expected_md5=None):
    """
    Ensure `path` exists and (optionally) matches `expected_md5`.
    - If the file exists and checksum matches (when provided), we return immediately.
    - Otherwise, we use gdown to download it from Drive (supports both cached_download with md5 or plain download).
    - If an expected_md5 is provided and post-download check fails, we delete the file and raise.
    """
    if os.path.exists(path):
        if expected_md5:
            try:
                if md5sum(path) == expected_md5:
                    return path  # File is present and correct
            except Exception:
                # If checksum computation fails for any reason, fall through to re-download
                pass
        else:
            return path  # File exists and no MD5 requested → accept it

    # File missing or MD5 mismatch → download
    if expected_md5:
        # cached_download will avoid re-downloading if it already has a matching cache entry and MD5
        gdown.cached_download(url=url, path=path, md5=expected_md5, quiet=False)
    else:
        # Plain download without MD5 verification
        gdown.download(url=url, output=path, quiet=False)

    # If we asked for MD5 checking, re-verify. If it still doesn't match, remove and raise.
    if expected_md5 and md5sum(path) != expected_md5:
        os.remove(path)
        raise RuntimeError("MD5 mismatch on downloaded file.")
    return path

@st.cache_resource(show_spinner=True)
def load_adata():
    """
    Lazily download (if needed) and load the AnnData file.
    - st.cache_resource caches the returned object across reruns and widget interactions.
    - show_spinner=True displays a spinner the first time this runs, which is nice UX for large files.
    """
    path = ensure_local_file(DRIVE_URL, LOCAL_PATH, EXPECTED_MD5)
    return sc.read(path)

# === App ===
# Title + brief context with citation
st.title("Single-Cell simple visualizer")
st.markdown(
    "Explore the single-cell RNAseq from mouse BAT, published in https://doi.org/10.1038/s42003-023-05140-2. "
    "Use filters on the left to refine your analysis."
)

# Two equivalent ways to load data:
#  - Using the cached function (recommended for production): `adata = load_adata()`
#  - Direct local read (used here): relies on the file already being present at LOCAL_PATH.
# The line below is commented to show the alternative approach without changing behavior.
# adata = load_adata()
adata = sc.read("shamsi_adata_demo_mini.h5ad")

# --- Sidebar filters ---
st.sidebar.header("Filters")

# Build the list of available cell types from AnnData.obs['cell_type_name'].
# - Drop NaNs to avoid UI noise and convert to a sorted Python list for deterministic ordering.
cell_types = sorted(adata.obs['cell_type_name'].dropna().unique().tolist())
# Example of a simpler multiselect (kept commented to preserve original code paths):
# selected_cell_types = st.sidebar.multiselect("Cell Type", cell_types, default=cell_types)

# --- TEST OF SELECT ALL ---
# We use st.session_state to persist selections across reruns.
# This is important because every interaction (e.g., button click) causes a rerun of the script.
if "selected_cell_types" not in st.session_state:
    st.session_state.selected_cell_types = cell_types  # default to "everything selected"

# Optional "Select all" checkbox as a convenience for users
select_all = st.sidebar.checkbox("Select all cell types", value=True)
if select_all:
    preselected = cell_types          # When checked, preselect all options
else:
    preselected = st.session_state.selected_cell_types  # Otherwise, use the last applied selection

# Multiselect control (does not auto-apply).
# Users can change the selection, but the filter only applies when they click the "Apply filter" button.
selected_cell_types = st.sidebar.multiselect(
    "Cell Type",
    cell_types,
    default=preselected
)

# Apply filter button:
# - Subsets the AnnData object and stores the filtered copy in session_state
# - We keep the selection as the "committed" state
if st.sidebar.button("Apply filter"):
    st.session_state.selected_cell_types = selected_cell_types
    # Important: .copy() creates a new AnnData object; avoids view/setwithcopy surprises later
    st.session_state.filtered = adata[adata.obs['cell_type_name'].isin(selected_cell_types)].copy()

# Clear filter button:
# - Removes the filtered object from session_state and resets selection to all cell types
if st.sidebar.button("Clear filter"):
    st.session_state.pop("filtered", None)
    st.session_state.selected_cell_types = cell_types

st.subheader("Summary data")
st.markdown("<hr style='margin-top:0; margin-bottom:10px;'>", unsafe_allow_html=True)

# Decide which AnnData to use downstream:
# - If a filtered dataset exists in session_state, use it and show a "Filtered dataset" summary
# - Otherwise, operate on the full dataset and show its dimensions
if "filtered" in st.session_state:
    adata_current = st.session_state.filtered
    # Note: this line assumes adata_current.raw is not None; otherwise .raw.n_vars would raise.
    # The dataset used here includes .raw, but keep this caveat in mind if you swap datasets.
    st.write(f"Filtered dataset: {adata_current.n_obs:,} cells × {adata_current.raw.n_vars:,} genes")
else:
    adata_current = adata
    # Same caveat about .raw applies here.
    st.write(f"Full dataset: {adata_current.n_obs:,} cells × {adata_current.raw.n_vars:,} genes")
# --- END OF TEST SELECT ALL ---

# Gene selector:
# - Use .raw.var_names when available (common in Scanpy pipelines: .raw holds log-normalized counts, etc.)
# - Fallback to .var_names if .raw is None
gene_list = sorted(adata.raw.var_names.tolist()) if adata.raw is not None else sorted(adata.var_names.tolist())
selected_gene = st.sidebar.selectbox("Gene of interest", gene_list)

# The code below shows an earlier alternative flow (kept commented to retain original structure):
# - It applied the filter immediately and then chose "current" dataset based on session_state
# # Apply filter button
# if st.sidebar.button("Apply filter"):
#     st.session_state.filtered = adata[adata.obs['cell_type_name'].isin(selected_cell_types)].copy()
#
# # Choose which dataset is "current"
# if "filtered" in st.session_state:
#     adata_current = st.session_state.filtered
#     st.write(f"**Filtered dataset** - showing {adata_current.n_obs} cells × {adata_current.raw.n_vars if adata_current.raw else adata_current.n_vars} genes")
# else:
#     adata_current = adata
#     st.write(f"**No filter applied** — showing all data: {adata_current.n_obs} cells × {adata_current.raw.n_vars if adata_current.raw else adata_current.n_vars} genes")



# --- Summary metrics ---
# Quick QC-like stats from .obs:
# - nCount_RNA: often total reads/UMIs per cell (naming varies by pipeline; guarded with dict check)
# - nFeature_RNA: number of detected genes per cell
mean_reads = adata_current.obs['nCount_RNA'].mean() if 'nCount_RNA' in adata_current.obs else float('nan')
umi_per_cell = adata_current.obs['nFeature_RNA'].mean() if 'nFeature_RNA' in adata_current.obs else float('nan')

# Streamlit metrics nicely display key figures at-a-glance in two columns
col1, col2 = st.columns(2)
col1.metric("Mean reads/cell", f"{mean_reads:,.0f}")
col2.metric("Mean genes/cell", f"{umi_per_cell:,.0f}")

st.subheader("Main plots")
st.markdown("<hr style='margin-top:0; margin-bottom:10px;'>", unsafe_allow_html=True)

# --- Main visualizations ---
# 1) UMAP:
#    - When a filter is applied, we plot the full adata and pass 'groups=selected_cell_types'
#      so Scanpy colors the selected groups and greys out the rest. This keeps global context.
#    - Otherwise, we plot the current (unfiltered) dataset, colored by 'cell_type_name'.
with plt.rc_context({"figure.figsize": (6, 6)}):
    if "filtered" in st.session_state:
        sc.pl.umap(adata, color="cell_type_name", groups=selected_cell_types, show=False, frameon=False)
    else:
        sc.pl.umap(adata_current, color="cell_type_name", show=False, frameon=False)
    plt.title("UMAP Plot", fontsize=14)   # <-- added title
    st.pyplot(plt.gcf())

# 2) Violin plot for the selected gene across cell types (using the current dataset)
#    - Common way to visualize expression distribution per cluster/category
if selected_gene:
    with plt.rc_context({"figure.figsize": (6, 4)}):
        sc.pl.violin(adata_current, keys=selected_gene, groupby="cell_type_name", show=False)
        plt.xticks(rotation=45, ha="right")
        plt.title(f"Gene Expression per Cell Type", fontsize=12)
        st.pyplot(plt.gcf())

# 3) Dot plot for the selected gene grouped by condition (if 'cond' metadata exists)
#    - Useful to compare expression patterns across experimental conditions
if selected_gene and "cond" in adata_current.obs:
    plot_width = max(4, len([selected_gene]) * 1.5)
    with plt.rc_context({"figure.figsize": (plot_width, 2)}):
        sc.pl.dotplot(adata_current, var_names=[selected_gene], groupby="cond", show=False)
        st.pyplot(plt.gcf(), use_container_width=False)  # prevents full-width stretching

