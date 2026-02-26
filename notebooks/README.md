## Tutorial 
```

Run DEG_comparison.ipynb (R)
Run Phenotype_evaluation.ipynb (R)
Run scRNA_scanpy.ipynb (python)
Run gene_expression_notebook.md (R)
Run single_cell_integration.ipynb (R)

```


## Code to reproduce analyses and figures 

### Data processing 
- [RNA-seq and WGS processing](rna_seq_data_processing.md)
- [Genome assessment](genome_analysis.md)

### Bulk expression analysis
- [Armadillo - human comparisons](human_notebook.md)
- [Gene-expression analysis](gene_expression_notebook.md)

### Validation 
- [DEG comparison with Hess 2022](DEG_comparison.ipynb)
- [Phenotype visualization](Phenotype_evaluation.ipynb)
- [Identity prediction](Identity_prediction.ipynb)

## Single-cell analysis

- scRNA-seq

- scATAC-seq

- Seurat (joint integration)

## Correspondence table 

| Figure 1 | File |
|----------|------|
| B-D      | [Phenotype visualization](Phenotype_evaluation.ipynb) |
| G-H      | [marker_heatmap_integrated.R](../scripts/multi_omics/marker_heatmap_integrated.R) |
| I        | [HVGs cross-quad](gene_expression_notebook.md) | 
    
| Figure 2 | File |
|----------|------|
| A-B      | [clustering_consistency.R](../scripts/sc_clustering/clustering_consistency.R) |
| D-F      | [gene_expression_notebook.md](gene_expression_notebook.md) |
| G-H      | [marker_heatmap_integrated.R](../scripts/multi_omics/marker_heatmap_integrated.R) |

| Figure 3 | File |
|----------|------|
| A-C      | [gene_expression_notebook.md](gene_expression_notebook.md) |
| E-F      | [gene_expression_notebook.md](gene_expression_notebook.md) |
| G        | [pred_identity.R](../scripts/multi_omics/pred_identity.R)  |
| H-I      | [analyze_gene_set_and_predictors.R](../scripts/multi_omics/analyze_gene_set_and_predictors.R)  |

| Figure 4 | File |
|----------|------|
| A-F        | [seurat_obj_analysis.R](../scripts/sc_clustering/seurat_obj_analysis.R)|
| G        | [pred_identity.R](../scripts/multi_omics/pred_identity.R)            |
| H-I      | [celltype_inference.R](../scripts/sc_clustering/celltype_inference.R)|
| J        | [snap_obj_analysis.R](../scripts/sc_clustering/snap_obj_analysis.R)  |


| Figure 5 | File |
|----------|------|
| A        | [seurat_obj_analysis.R](../scripts/sc_clustering/seurat_obj_analysis.R)|
| B-C      | [analyze_gene_set_and_predictors.R](../scripts/multi_omics/analyze_gene_set_and_predictors.R)|
| D-E      | [marker_heatmap_integrated.R](../scripts/multi_omics/marker_heatmap_integrated.R) |

## Supplementary Figure

| Figure   | File |
|----------|------|
| Figure S1 | [human_notebook.md](humannotebook.md)                      |
| Figure S2 | [gene_expression_notebook.md](gene_expression_notebook.md) |
| Figure S3 | [DEG comparison with Hess 2022](DEG_comparison.ipynb)      |
| Figure S4 | [gene_expression_notebook.md](gene_expression_notebook.md) |
| Figure S5 | [gene_expression_notebook.md](gene_expression_notebook.md) |
| Figure S6 | [gene_expression_notebook.md](gene_expression_notebook.md) |
| Figure S7 | [analyze_gene_set_and_predictors.R](../scripts/multi_omics/analyze_gene_set_and_predictors.R) |
| Figure S8 | [analyze_gene_set_and_predictors.R](../scripts/multi_omics/analyze_gene_set_and_predictors.R) |
| Figure S9 | [pred_identity.R](../scripts/multi_omics/pred_identity.R) |
| Figure S10 | [scRNA_scanpy.ipynb](scRNA_scanpy.ipynb)
| Figure S11 | [snap_basic_stats.R](../scripts/sc_clustering/snap_basic_stats.R) and bamPEFragmentSize (deeptools)|
| Figure S12 | [snap_basic_stats.R](../scripts/sc_clustering/snap_basic_stats.R) |
| Figure S13 | [Phenotype_evaluation.ipynb](Phenotype_evaluation.ipynb) |
| Figure S14 | [snap_obj_analysis.R](../scripts/sc_clustering/snap_obj_analysis.R) |
| Figure S15 | [seurat_obj_analysis.R](../scripts/sc_clustering/seurat_obj_analysis.R) |




