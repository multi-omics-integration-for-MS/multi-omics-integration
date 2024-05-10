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

    - [ ] 12_transcriptomic_dataset.ipynb
        - [ ] add new genes to GSE194078
        - [ ] concatenate all transcriptomic dataset
        - [X] add cell types labels
        - [X] split CSF and PBMC
        - [X] split HC and patiens
        - [X] split train and test with stratification on labels and HC/P (both for CSF and PBMC) (NOT TO DO)

    - [ ] 13_EDA_transcriptomic_dataset.ipynb
        - [ ] EDA tutto dataset

- [ ] differential expression analysis

- [ ] gene integration

- [X] MRI? (NOT TO DO)
    - [ ] 14_MRI_data.ipynb
    - [ ]brain lesion also in GSE173787

- [X] Metabolomic? (NOT TO DO)
