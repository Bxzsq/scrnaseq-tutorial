
process ANNOTATION_AUTO {

    tag "CellTypist_auto"

    container 'celltypist-auto:1.0'

    cpus 4
    memory '32 GB'

    publishDir "${params.outdir}/annotation_auto", mode: 'copy'

    input:
    path matrix_h5
    path samplesheet
    path annotation_script

    output:
    path "celltypist_predictions.csv"
    path "celltypist_celltype_counts.csv"
    path "celltypist_annotated.h5ad"

    path "celltypist_umap.pdf"
    path "celltypist_major_umap.pdf"
    path "celltypist_confidence_umap.pdf"
    path "celltypist_library_umap.pdf"

    script:
    """
    echo "============================================================"
    echo "CellTypist automatic annotation"
    echo "============================================================"

    echo "Input H5 files:"
    ls -lh *.h5

    echo
    echo "Running CellTypist..."

    python ${annotation_script} \
        --h5_dir . \
        --samplesheet ${samplesheet} \
        --model Adult_Human_PrefrontalCortex.pkl \
        --output celltypist_annotated.h5ad \
        --predictions celltypist_predictions.csv \
        --counts celltypist_celltype_counts.csv

    echo
    echo "CellTypist annotation finished."

    echo
    echo "Output files:"
    ls -lh *.csv *.h5ad *.pdf
    """
}
