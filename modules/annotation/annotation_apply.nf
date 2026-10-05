process ANNOTATION_APPLY {

    tag "Apply_manual_annotation"

    container 'satijalab/seurat:5.5.1'

    cpus 4
    memory '32 GB'

    publishDir "${params.outdir}/annotation", mode: 'copy'

    input:
    path seurat_rds
    path annotation_file
    path annotation_apply_script

    output:
    path "annotated_manual.rds", emit: annotated_rds
    path "celltype_annotation.csv"
    path "celltype_umap.pdf"

    script:
    """
    echo "============================================================"
    echo "Apply manual cell type annotation"
    echo "============================================================"

    echo "Input Seurat object:"
    ls -lh ${seurat_rds}

    echo
    echo "Annotation file:"
    cat ${annotation_file}

    echo
    echo "Annotation script:"
    ls -lh ${annotation_apply_script}

    Rscript ${annotation_apply_script} \
        --input ${seurat_rds} \
        --annotation ${annotation_file} \
        --output annotated_manual.rds \
        --annotation_output celltype_annotation.csv \
        --umap_output celltype_umap.pdf

    echo
    echo "Manual annotation applied successfully."
    """
}