### TODO

- final cell type labels: 'results/celltypist_labels_from_protein_Immune_All_High.csv'
- geni più rilevanti per ogni cell type in 'results/celltypist_markers_from_protein_Immune_All_High.json' (50 top genes)

- [X] 11_BBKNN_GSE239626_GSE138266.ipynb
    - [X] Memory error in sc.tl.ingest: FIX
    - [X] or find an other way to integrate the new dataset with the reference dataset
    - [X] tempo ingest+BBKNN su 1/8 dataset: 25 min

- [X] dopo aver controllato labels in 08_cell_types_from_proteins_clustering.ipynb:

    - [X] 03_GSE239626_labels_from_proteins.ipynb
        - [X] Infer cell types for dataset GSE23962 from clustering protein and save in csv

    - [X] BBKNN (NON FARE, non integra correttamente i dataset!)
        - [X] run 10_BBKNN_GSE239626_GSE194078.ipynb and 11_BBKNN_GSE239626_GSE138266.ipynb
        - [X] SAVE labels for datasets GSE194078 and GSE138266
    
    - [X] BBKNN non integra correttamente, usare altra lib (in R)

    - [X] 12_transcriptomic_dataset.ipynb
        - [X] add new genes to GSE194078
        - [X] concatenate all transcriptomic dataset
        - [X] add cell types labels
        - [X] split CSF and PBMC
        - [X] split HC and patiens
        - [X] split train and test with stratification on labels and HC/P (both for CSF and PBMC) (NOT TO DO)

    - [X] 13_EDA_transcriptomic_dataset.ipynb
        - [X] EDA tutto dataset
        
- [ ] cambiare parametri clustering (?) in notebook 14_clustering_and_finding_markers_gene.ipynb

- [X] differential expression analysis: 16_differential_expression_analysis.ipynb
- [ ] ridurre rumore cellule e ri runnare nb 16_differential_expression_analysis.ipynb (o usare solo cell type più frequenti)

- [X] differential expression analysis (all cell types): NOTEBOOK 17, 18, 19, 20 + results 21_markers.ipynb

- [ ] gene integration: 25_microRNA_GSE159033.ipynb

- [X] MRI? (NOT TO DO)
    - [ ] MRI_data.ipynb
    - [ ] brain lesion also in GSE173787

- [X] Metabolomic? (NOT TO DO)

- [ ] REPORT, files:
    - clustering: 01_proteins_clustering.ipynb + 02_assess_clustering_robustness.ipynb
    - labels: 03_GSE239626_labels_from_proteins.ipynb
    - labels: 10_BBKNN_GSE239626_GSE194078.ipynb + 11_BBKNN_GSE239626_GSE138266.ipynb
    - EDA: 13_EDA_transcriptomic_dataset.ipynb
    - differential expression analysis: 22_differential_expression_per_cell_types.ipynb.ipynb + 23_most_expressive_genes.ipynb
        thresholds: p-value adjusted < 1e-50, abs(logFC) > 0.5
    - gene integration (T cells): 24_T_cells_genes_integration_GSE173787.ipynb