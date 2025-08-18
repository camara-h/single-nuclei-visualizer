# -------------------------------------------------------------------------------------------------
# Single-Cell Viewer (Streamlit + Scanpy) — simplified version
#
# What this app does:
# - Loads a demo AnnData (.h5ad) single-cell RNA-seq dataset (mouse BAT)
# - Lets you subset cells by cell type via a sidebar multiselect
# - Provides quick summary metrics (mean reads/genes per cell)
# - Visualizes UMAP, a violin plot for a selected gene, and a dot plot by condition
# - Allows upload of a full .h5ad dataset (optional)
# -------------------------------------------------------------------------------------------------

import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt

# Configure Streamlit page
st.set_page_config(page_title="Single-Cell Viewer", layout="wide")

# === App title and description ===
st.title("Single-Cell simple visualizer")
st.markdown(
    "Explore the single-cell RNAseq from mouse BAT, published in https://doi.org/10.1038/s42003-023-05140-2. "
    "Upload your own dataset or use the demo below."
)

# === Data loading ===
uploaded_file = st.file_uploader("Upload your full .h5ad dataset")
if uploaded_file is not None:
    adata = sc.read(uploaded_file)
else:
    adata = sc.read("shamsi_adata_demo_mini.h5ad")

# === Sidebar filters ===
st.sidebar.header("Filters")

# Cell type options
cell_types = sorted(adata.obs['cell_type_name'].dropna().unique().tolist())

# Keep last selection in session_state
if "selected_cell_types" not in st.session_state:
    st.session_state.selected_cell_types = cell_types

# Multiselect (user must hit Apply)
selected_cell_types = st.sidebar.multiselect(
    "Cell Type",
    cell_types,
    default=st.session_state.selected_cell_types
)

# Apply filter button
if st.sidebar.button("Apply filter"):
    st.session_state.selected_cell_types = selected_cell_types
    st.session_state.filtered = adata[adata.obs['cell_type_name'].isin(selected_cell_types)].copy()

# Decide which AnnData to use
if "filtered" in st.session_state:
    adata_current = st.session_state.filtered
    st.write(f"Filtered dataset: {adata_current.n_obs:,} cells × {adata_current.raw.n_vars:,} genes")
else:
    adata_current = adata
    st.write(f"Full dataset: {adata_current.n_obs:,} cells × {adata_current.raw.n_vars:,} genes")

# Gene selector
gene_list = sorted(adata.raw.var_names.tolist()) if adata.raw is not None else sorted(adata.var_names.tolist())
selected_gene = st.sidebar.selectbox("Gene of interest", gene_list)

# === Summary data ===
st.subheader("Summary data")
st.markdown("<hr style='margin-top:0; margin-bottom:10px;'>", unsafe_allow_html=True)

mean_reads = adata_current.obs['nCount_RNA'].mean() if 'nCount_RNA' in adata_current.obs else float('nan')
umi_per_cell = adata_current.obs['nFeature_RNA'].mean() if 'nFeature_RNA' in adata_current.obs else float('nan')

col1, col2 = st.columns(2)
col1.metric("Mean reads/cell", f"{mean_reads:,.0f}")
col2.metric("Mean genes/cell", f"{umi_per_cell:,.0f}")

# === Main plots ===
st.subheader("Main plots")
st.markdown("<hr style='margin-top:0; margin-bottom:10px;'>", unsafe_allow_html=True)

# 1) UMAP
with plt.rc_context({"figure.figsize": (6, 6)}):
    if "filtered" in st.session_state:
        sc.pl.umap(adata, color="cell_type_name", groups=selected_cell_types, show=False, frameon=False)
    else:
        sc.pl.umap(adata_current, color="cell_type_name", show=False, frameon=False)
    plt.title("UMAP Plot", fontsize=14)
    st.pyplot(plt.gcf())

# 2) Violin plot
if selected_gene:
    with plt.rc_context({"figure.figsize": (6, 4)}):
        sc.pl.violin(adata_current, keys=selected_gene, groupby="cell_type_name", show=False)
        plt.xticks(rotation=45, ha="right")
        plt.title(f"Gene Expression per Cell Type", fontsize=12)
        st.pyplot(plt.gcf())

# 3) Dot plot
if selected_gene and "cond" in adata_current.obs:
    plot_width = max(4, len([selected_gene]) * 1.5)
    with plt.rc_context({"figure.figsize": (plot_width, 2)}):
        sc.pl.dotplot(adata_current, var_names=[selected_gene], groupby="cond", show=False)
        st.pyplot(plt.gcf(), use_container_width=False)
