process ANNOTATION {

    tag "Cell_annotation"

    container 'satijalab/seurat:5.5.1'

    cpus 4
    memory '32 GB'

    publishDir "${params.outdir}/annotation", mode: 'copy'

    input:
    path seurat_rds
    path annotation_script

    output:
    path "annotated_manual.rds", emit: annotated_rds
    path "cluster_markers.csv"
    path "top20_markers.csv"
    path "marker_dotplot.pdf"
    path "marker_heatmap.pdf"
    path "cluster_annotation.csv"

    script:
    """
    echo "============================================================"
    echo "Cell type annotation"
    echo "============================================================"

    Rscript ${annotation_script} \
        --input ${seurat_rds} \
        --output annotated_manual.rds \
        --markers cluster_markers.csv \
        --top_markers top20_markers.csv \
        --dotplot marker_dotplot.pdf \
        --heatmap marker_heatmap.pdf \
        --annotation cluster_annotation.csv

    echo
    echo "Cell annotation finished."
    """
}