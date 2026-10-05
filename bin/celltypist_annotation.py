
import argparse
import glob
import os
import re

import numpy as np
import pandas as pd
import scanpy as sc
import celltypist
from celltypist import models


def calc_percent_mt(adata):
    """
    Calculate mitochondrial percentage.
    Equivalent to Seurat:
    PercentageFeatureSet(object, pattern = "^MT-")
    """

    mt_mask = adata.var_names.str.upper().str.startswith("MT-")

    total_counts = np.asarray(
        adata.X.sum(axis=1)
    ).ravel()

    if mt_mask.sum() > 0:

        mt_counts = np.asarray(
            adata[:, mt_mask].X.sum(axis=1)
        ).ravel()

        percent_mt = np.divide(
            mt_counts,
            total_counts,
            out=np.zeros_like(total_counts, dtype=float),
            where=total_counts > 0
        ) * 100

    else:

        percent_mt = np.zeros(
            adata.n_obs,
            dtype=float
        )

    return percent_mt


def make_major_type(label):

    """
    Convert CellTypist fine-grained labels
    into major cell types.

    Examples:
        Oligo MOG OPALIN -> Oligodendrocyte
        Micro P2RY12 CCL3 -> Microglia
        Astro GFAP FABP7 -> Astrocyte
    """

    label = str(label)

    if label.startswith("Oligo"):
        return "Oligodendrocyte"

    elif label.startswith("OPC"):
        return "OPC"

    elif label.startswith("COP"):
        return "COP"

    elif label.startswith("Astro"):
        return "Astrocyte"

    elif label.startswith("Micro"):
        return "Microglia"

    elif label.startswith("Macro"):
        return "Macrophage"

    elif label.startswith("Endo"):
        return "Endothelial"

    elif label.startswith("VLMC"):
        return "VLMC"

    elif label.startswith("SMC"):
        return "SMC"

    elif label.startswith("Per"):
        return "Pericyte"

    elif label.startswith("PC"):
        return "Pericyte"

    elif label.startswith("L2"):
        return "Excitatory neuron"

    elif label.startswith("L3"):
        return "Excitatory neuron"

    elif label.startswith("L4"):
        return "Excitatory neuron"

    elif label.startswith("L5"):
        return "Excitatory neuron"

    elif label.startswith("L6"):
        return "Excitatory neuron"

    elif label.startswith("InN"):
        return "Inhibitory neuron"

    elif label.startswith("B"):
        return "B cell"

    elif label.startswith("T"):
        return "T cell"

    elif label.startswith("Myeloid"):
        return "Myeloid"

    elif label.startswith("RB"):
        return "Red blood cell"

    else:
        return "Other"


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--h5_dir",
        required=True
    )

    parser.add_argument(
        "--samplesheet",
        required=True
    )

    parser.add_argument(
        "--model",
        required=True
    )

    parser.add_argument(
        "--output",
        required=True
    )

    parser.add_argument(
        "--predictions",
        required=True
    )

    parser.add_argument(
        "--counts",
        required=True
    )

    args = parser.parse_args()

    print("=" * 70)
    print("CellTypist automatic annotation")
    print("=" * 70)

    # ============================================================
    # 1. Find H5 files
    # ============================================================

    h5_files = sorted(
        glob.glob(
            os.path.join(
                args.h5_dir,
                "*_filtered_feature_bc_matrix.h5"
            )
        )
    )

    if len(h5_files) == 0:
        raise RuntimeError(
            "No *_filtered_feature_bc_matrix.h5 files found."
        )

    print("\nH5 files:")

    for f in h5_files:
        print("  ", os.path.basename(f))

    print(
        f"\nTotal H5 files: {len(h5_files)}"
    )

    # ============================================================
    # 2. Read samplesheet
    # ============================================================

    sample_df = pd.read_csv(
        args.samplesheet
    )

    print("\nSamplesheet:")
    print(sample_df)

    condition_map = dict(
        zip(
            sample_df["library"].astype(str),
            sample_df["condition"].astype(str)
        )
    )

    # ============================================================
    # 3. Read H5 files
    # ============================================================

    adata_list = []

    for h5 in h5_files:

        basename = os.path.basename(h5)

        library = re.sub(
            r"_filtered_feature_bc_matrix\.h5$",
            "",
            basename
        )

        print("\n" + "-" * 70)
        print(f"Reading: {basename}")
        print(f"Library: {library}")

        adata = sc.read_10x_h5(
            h5,
            gex_only=True
        )

        print(
            f"Original: "
            f"{adata.n_obs} cells × "
            f"{adata.n_vars} genes"
        )

        # --------------------------------------------------------
        # Make gene names unique
        # --------------------------------------------------------

        adata.var_names_make_unique()

        # --------------------------------------------------------
        # Prefix cell barcode with library
        # --------------------------------------------------------

        adata.obs_names = [
            f"{library}_{barcode}"
            for barcode in adata.obs_names
        ]

        adata.obs["library"] = library

        adata.obs["condition"] = condition_map.get(
            library,
            "Unknown"
        )

        # --------------------------------------------------------
        # QC metrics
        # --------------------------------------------------------

        adata.obs["nCount_RNA"] = np.asarray(
            adata.X.sum(axis=1)
        ).ravel()

        adata.obs["nFeature_RNA"] = np.asarray(
            (adata.X > 0).sum(axis=1)
        ).ravel()

        adata.obs["percent.mt"] = calc_percent_mt(
            adata
        )

        print(
            f"Before QC: {adata.n_obs} cells"
        )

        # --------------------------------------------------------
        # Same QC as Seurat
        #
        # nFeature_RNA >= 200
        # nFeature_RNA <= 3000
        # percent.mt <= 5
        # --------------------------------------------------------

        keep = (
            (adata.obs["nFeature_RNA"] >= 200)
            &
            (adata.obs["nFeature_RNA"] <= 3000)
            &
            (adata.obs["percent.mt"] <= 5)
        )

        adata = adata[keep].copy()

        print(
            f"After QC: {adata.n_obs} cells"
        )

        adata_list.append(adata)

    # ============================================================
    # 4. Merge all libraries
    # ============================================================

    print("\n" + "=" * 70)
    print("Merging libraries")
    print("=" * 70)

    adata = sc.concat(
        adata_list,
        join="outer",
        merge="same",
        index_unique=None
    )

    print(
        f"Merged object: "
        f"{adata.n_obs} cells × "
        f"{adata.n_vars} genes"
    )

    # ============================================================
    # 5. Preserve raw counts
    # ============================================================

    adata.layers["counts"] = adata.X.copy()

    # ============================================================
    # 6. Normalize
    # ============================================================

    print("\nNormalizing expression...")

    sc.pp.normalize_total(
        adata,
        target_sum=10000
    )

    sc.pp.log1p(adata)

    # ============================================================
    # 7. CellTypist model
    # ============================================================

    print("\n" + "=" * 70)
    print("Loading CellTypist model")
    print("=" * 70)

    model = models.Model.load(
        model=args.model
    )

    print(model)

    # ============================================================
    # 8. Gene overlap
    # ============================================================

    model_genes = set(
        model.features
    )

    data_genes = set(
        adata.var_names
    )

    overlap = len(
        model_genes.intersection(data_genes)
    )

    print(
        f"\nModel genes: {len(model_genes)}"
    )

    print(
        f"Dataset genes: {len(data_genes)}"
    )

    print(
        f"Gene overlap: {overlap}"
    )

    if overlap < 1000:

        raise RuntimeError(
            f"Gene overlap is too low: "
            f"{overlap} genes."
        )

    # ============================================================
    # 9. CellTypist annotation
    # ============================================================

    print("\n" + "=" * 70)
    print("Running CellTypist")
    print("=" * 70)

    predictions = celltypist.annotate(
        adata,
        model=model,
        majority_voting=False
    )

    print(
        "\nCellTypist annotation finished."
    )

    # ============================================================
    # 10. Extract prediction
    # ============================================================

    pred_df = (
        predictions.predicted_labels.copy()
    )

    print("\nPrediction table:")
    print(pred_df.head())

    if "predicted_labels" in pred_df.columns:

        label_col = "predicted_labels"

    else:

        label_col = pred_df.columns[0]

    if "conf_score" in pred_df.columns:

        conf_col = "conf_score"

    else:

        conf_col = None

    adata.obs["celltypist_label"] = (
        pred_df[label_col].values
    )

    if conf_col is not None:

        adata.obs["celltypist_conf_score"] = (
            pred_df[conf_col].values
        )

    else:

        adata.obs["celltypist_conf_score"] = np.nan

    # ============================================================
    # 11. Major cell type
    # ============================================================

    adata.obs["celltypist_major"] = (
        adata.obs["celltypist_label"]
        .map(make_major_type)
    )

    # ============================================================
    # 12. Prediction table
    # ============================================================

    result = adata.obs[
        [
            "library",
            "condition",
            "nCount_RNA",
            "nFeature_RNA",
            "percent.mt",
            "celltypist_label",
            "celltypist_major",
            "celltypist_conf_score"
        ]
    ].copy()

    result.insert(
        0,
        "cell",
        result.index
    )

    result.to_csv(
        args.predictions,
        index=False
    )

    # ============================================================
    # 13. Cell type counts
    # ============================================================

    celltype_counts = (
        result["celltypist_label"]
        .value_counts()
        .rename_axis("cell_type")
        .reset_index(name="n_cells")
    )

    celltype_counts.to_csv(
        args.counts,
        index=False
    )

    # ============================================================
    # 14. PCA
    # ============================================================

    print("\n" + "=" * 70)
    print("Computing PCA")
    print("=" * 70)

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=2000,
        flavor="seurat"
    )

    sc.pp.scale(
        adata,
        max_value=10
    )

    sc.tl.pca(
        adata,
        n_comps=30,
        use_highly_variable=True,
        svd_solver="arpack"
    )

    # ============================================================
    # 15. Neighbors
    # ============================================================

    print("\nComputing neighbors...")

    sc.pp.neighbors(
        adata,
        n_neighbors=15,
        n_pcs=30
    )

    # ============================================================
    # 16. UMAP
    # ============================================================

    print("\nComputing UMAP...")

    sc.tl.umap(
        adata,
        random_state=42
    )

    # ============================================================
    # 17. UMAP: fine-grained CellTypist labels
    # ============================================================

    print(
        "\nGenerating CellTypist UMAP..."
    )

    sc.pl.umap(
        adata,
        color="celltypist_label",
        legend_loc="right margin",
        frameon=False,
        size=8,
        title="CellTypist annotation",
        show=False,
        save="_celltypist.pdf"
    )

    # ============================================================
    # 18. UMAP: major cell types
    # ============================================================

    print(
        "Generating major cell type UMAP..."
    )

    sc.pl.umap(
        adata,
        color="celltypist_major",
        legend_loc="right margin",
        frameon=False,
        size=8,
        title="CellTypist major cell types",
        show=False,
        save="_celltypist_major.pdf"
    )

    # ============================================================
    # 19. UMAP: confidence
    # ============================================================

    print(
        "Generating confidence UMAP..."
    )

    sc.pl.umap(
        adata,
        color="celltypist_conf_score",
        color_map="viridis",
        frameon=False,
        size=8,
        title="CellTypist confidence score",
        show=False,
        save="_celltypist_confidence.pdf"
    )

    # ============================================================
    # 20. UMAP: library
    # ============================================================

    print(
        "Generating library UMAP..."
    )

    sc.pl.umap(
        adata,
        color="library",
        legend_loc="right margin",
        frameon=False,
        size=8,
        title="Library",
        show=False,
        save="_library.pdf"
    )

    # ============================================================
    # 21. Rename Scanpy generated figures
    # ============================================================

    scanpy_fig_dir = os.path.join(
        "figures"
    )

    if os.path.exists(scanpy_fig_dir):

        rename_map = {

            "umap_celltypist.pdf":
                "celltypist_umap.pdf",

            "umap_celltypist_major.pdf":
                "celltypist_major_umap.pdf",

            "umap_celltypist_confidence.pdf":
                "celltypist_confidence_umap.pdf",

            "umap_library.pdf":
                "celltypist_library_umap.pdf"
        }

        for old_name, new_name in rename_map.items():

            old_path = os.path.join(
                scanpy_fig_dir,
                old_name
            )

            new_path = os.path.join(
                ".",
                new_name
            )

            if os.path.exists(old_path):

                os.rename(
                    old_path,
                    new_path
                )

    # ============================================================
    # 22. Save annotated AnnData
    # ============================================================

    print(
        "\nSaving annotated AnnData..."
    )

    adata.write_h5ad(
        args.output
    )

    # ============================================================
    # 23. Summary
    # ============================================================

    print("\n" + "=" * 70)
    print("CellTypist summary")
    print("=" * 70)

    print(
        f"Final cells: {adata.n_obs}"
    )

    print(
        f"Genes: {adata.n_vars}"
    )

    print("\nMajor cell types:")

    print(
        adata.obs["celltypist_major"]
        .value_counts()
        .to_string()
    )

    print("\nOutput files:")
    print(
        f"  {args.predictions}"
    )

    print(
        f"  {args.counts}"
    )

    print(
        f"  {args.output}"
    )

    print(
        "  celltypist_umap.pdf"
    )

    print(
        "  celltypist_major_umap.pdf"
    )

    print(
        "  celltypist_confidence_umap.pdf"
    )

    print(
        "  celltypist_library_umap.pdf"
    )

    print("\nFinished.")


if __name__ == "__main__":
    main()
