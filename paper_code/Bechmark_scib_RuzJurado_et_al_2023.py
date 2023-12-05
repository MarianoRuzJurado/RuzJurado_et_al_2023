import warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)
import scanpy as sc
import pandas as pd
import numpy as np
import scib
import os
import anndata as ad


'For the metrics need an unintegrated AnnData Object from my samples'
folder_path = '/media/Storage/R/Revision_scib/Conversion/unintegrated'
ann_data_dict  = {}
sce_obj_list = os.listdir(folder_path)
for file_name in sce_obj_list:
    if file_name.endswith(".h5ad"):
        name=os.path.splitext(file_name)[0]

        file_path = os.path.join(folder_path, file_name)
        adata = sc.read_h5ad(file_path)

        ann_data_dict[name] = adata

adata_unint = ad.concat(ann_data_dict.values(), join="outer", label='orig.ident', keys=ann_data_dict.keys())
sc.pp.pca(adata_unint, n_comps=30)



'CALCULATE METRICS FOR ORTHOINTEGRATE'

# Read in sce object, convert to anndata obj
adata = sc.read_h5ad('/media/Storage/R/Revision_scib/Conversion/SeuratObject.combined.integrated.annotated_10_08_22.h5ad')

# Read in UMAP embeddings
UMAPEmbeddings = pd.read_csv('/media/Storage/R/Integrated/UMAPEmbedings.csv', skiprows=0, header=0 ,index_col=0)
'Remake the Seurat Clustering since its not saved in the SCE object'
PCAEmbeddings = pd.read_csv('/media/Storage/R/Integrated/PCAEmbedings.csv', skiprows=0, header=0, index_col=0)

# Check order
tmpcheck = UMAPEmbeddings.index == adata.obs_names
np.unique(tmpcheck, return_counts=True)
tmpcheck = PCAEmbeddings.index == adata.obs_names
np.unique(tmpcheck, return_counts=True)

#Assign UMAP parameters in anndata object
adata.obsm['X_umap'] = UMAPEmbeddings.values
sc.pl.umap(adata, color=['species'], legend_loc='on data')
adata.obsm['X_pca'] = PCAEmbeddings.values

#scib package metric calculation for principal component regression and silhouette score per batch
pcr_comparison_ortho = scib.metrics.pcr_comparison(adata_unint, adata, covariate='sample', n_comps=30)
silhouetteBatch_ortho = scib.metrics.silhouette_batch(adata, batch_key="sample", label_key="cell_type", embed="X_umap", return_all=True)

#calculate graph connectivity score
sc.pp.neighbors(adata, n_pcs=10)
graph_connectivity_ortho = scib.metrics.graph_connectivity(adata, label_key="cell_type")

#calculate average of k-nearest neighbor batch effect on downsampled score, like proposed by Büttner et al. (2019) for big datasets
subset_size = 0.1
num_cells_to_subsample = int(subset_size * adata.n_obs)
subset_indices = np.random.choice(adata.n_obs, num_cells_to_subsample, replace=False)
subsampled_adata = adata[subset_indices, :]
kBET_ortho = scib.metrics.kBET(subsampled_adata, batch_key="sample", label_key="cell_type", type_="embed", embed="X_umap")
LISI_score_ortho = scib.metrics.ilisi_graph(adata, batch_key="sample", type_="embed", use_rep="X_umap", n_cores=30)

#Biological Conservation Metrics
cell_cycle_ortho = scib.metrics.cell_cycle(adata_unint,adata, batch_key="sample", organism='human')
clisi_graph_ortho = scib.metrics.clisi_graph(adata, label_key='cell_type', type_='embed', use_rep='X_umap',n_cores=40)
isolated_labels_f1_ortho = scib.metrics.isolated_labels_f1(subsampled_adata, label_key="cell_type", batch_key="sample", embed='X_umap', iso_threshold=50)



#Build Dataframe with results for OrthoIntegrate
output_ortho = pd.DataFrame()
output_ortho.loc['graph_connectivity', 'OrthoIntegrate'] = graph_connectivity_ortho
output_ortho.loc['pcr_comparison', 'OrthoIntegrate'] = pcr_comparison_ortho # lower scores mean less impact of each sample to our pca analysis, so it become more robust to individual samples
output_ortho.loc['silhouette_batch', 'OrthoIntegrate'] = silhouetteBatch_ortho[0]
output_ortho.loc['kBET', 'OrthoIntegrate'] = kBET_ortho
output_ortho.loc['LISI', 'OrthoIntegrate'] = LISI_score_ortho
output_ortho.loc['Cell cycle conservation', 'OrthoIntegrate'] = cell_cycle_ortho
output_ortho.loc['Cell-type LISI', 'OrthoIntegrate'] = clisi_graph_ortho
output_ortho.loc['isolated_labels_f1', 'OrthoIntegrate'] = isolated_labels_f1_ortho
output_ortho.loc['Species mixing score', 'OrthoIntegrate'] = (graph_connectivity_ortho + pcr_comparison_ortho + silhouetteBatch_ortho[0] + kBET_ortho + LISI_score_ortho)/5 # Species Score
output_ortho.T.to_csv('/media/Storage/R/Integrated/metrics_OrthoIntegrate_cell_type.csv')

'CALCULATE METRICS FOR Inparanoid'

# Read in sce object, convert to anndata obj
adata_para = sc.read_h5ad('/media/Storage/R/Revision_scib/Conversion/SeuratObject.inparanoid_10_08_22.h5ad')
# Read in UMAP embeddings
UMAPEmbeddings_para = pd.read_csv('/media/Storage/R/Inparanoid/UMAPEmbedings.csv', skiprows=0, header=0 ,index_col=0)
'Remake the Seurat Clustering since its not saved in the SCE object'
PCAEmbeddings_para = pd.read_csv('/media/Storage/R/Inparanoid/PCAEmbedings.csv', skiprows=0, header=0, index_col=0)

# Check order
tmpcheck = UMAPEmbeddings_para.index == adata_para.obs_names
np.unique(tmpcheck, return_counts=True)
tmpcheck = PCAEmbeddings_para.index == adata_para.obs_names
np.unique(tmpcheck, return_counts=True)

#Assign UMAP parameters in anndata object
adata_para.obsm['X_umap'] = UMAPEmbeddings_para.values
sc.pl.umap(adata, color=['species'], legend_loc='on data')
adata_para.obsm['X_pca'] = PCAEmbeddings_para.values

#scib package metric calculation for principal component regression and silhouette score per batch
pcr_comparison_para = scib.metrics.pcr_comparison(adata_unint, adata_para, covariate='sample', n_comps=30)
silhouetteBatch_para = scib.metrics.silhouette_batch(adata_para, batch_key="sample", label_key="cell_type", embed="X_umap", return_all=True)

#calculate graph connectivity score
sc.pp.neighbors(adata_para, n_pcs=10)
graph_connectivity_para = scib.metrics.graph_connectivity(adata_para, label_key="cell_type")

#calculate average of k-nearest neighbor batch effect on downsampled score, like proposed by Büttner et al. (2019) for big datasets
subset_size = 0.1
num_cells_to_subsample = int(subset_size * adata_para.n_obs)
subset_indices = np.random.choice(adata_para.n_obs, num_cells_to_subsample, replace=False)
subsampled_adata = adata_para[subset_indices, :]
kBET_para = scib.metrics.kBET(subsampled_adata, batch_key="sample", label_key="cell_type", type_="embed", embed="X_umap")

LISI_score_para = scib.metrics.ilisi_graph(adata_para, batch_key="sample", type_="embed", use_rep="X_umap", n_cores=30)
cell_cycle_para = scib.metrics.cell_cycle(adata_unint,adata_para, batch_key="sample", organism='human')
clisi_graph_para = scib.metrics.clisi_graph(adata_para, label_key='cell_type', type_='embed', use_rep='X_umap',n_cores=40)
isolated_labels_f1_para = scib.metrics.isolated_labels_f1(subsampled_adata, label_key="cell_type", batch_key="sample", embed='X_umap', iso_threshold=50)

#Build Dataframe with results for inParanoid
output_para = pd.DataFrame()
output_para.loc['graph_connectivity', 'inParanoid'] = graph_connectivity_para
output_para.loc['pcr_comparison', 'inParanoid'] = pcr_comparison_para # lower scores mean less impact of each sample to our pca analysis, so it become more robust to individual samples
output_para.loc['silhouette_batch', 'inParanoid'] = silhouetteBatch_para[0]
output_para.loc['kBET', 'inParanoid'] = kBET_para
output_para.loc['LISI', 'inParanoid'] = LISI_score_para
output_para.loc['Cell cycle conservation', 'inParanoid'] = cell_cycle_para
output_para.loc['Cell-type LISI', 'inParanoid'] = clisi_graph_para
output_para.loc['isolated_labels_f1', 'inParanoid'] = isolated_labels_f1_para
output_para.loc['Species mixing score', 'inParanoid'] = (graph_connectivity_para + pcr_comparison_para + silhouetteBatch_para[0] + kBET_para + LISI_score_para)/5 # Species Score
output_para.T.to_csv('/media/Storage/R/Inparanoid/metrics_Inparanoid_cell_type.csv')

'CALCULATE METRICS FOR just_biomart'

# Read in sce object, convert to anndata obj
adata_biomart = sc.read_h5ad('/media/Storage/R/Revision_scib/Conversion/SeuratObject.biomart_10_08_22.h5ad')
# Read in UMAP embeddings
UMAPEmbeddings_biomart = pd.read_csv('/media/Storage/R/just_biomart/UMAPEmbedings.csv', skiprows=0, header=0 ,index_col=0)
'Remake the Seurat Clustering since its not saved in the SCE object'
PCAEmbeddings_biomart = pd.read_csv('/media/Storage/R/just_biomart/PCAEmbedings.csv', skiprows=0, header=0, index_col=0)

# Check order
tmpcheck = UMAPEmbeddings_biomart.index == adata_biomart.obs_names
np.unique(tmpcheck, return_counts=True)
tmpcheck = PCAEmbeddings_biomart.index == adata_biomart.obs_names
np.unique(tmpcheck, return_counts=True)

#Assign UMAP parameters in anndata object
adata_biomart.obsm['X_umap'] = UMAPEmbeddings_biomart.values
sc.pl.umap(adata_biomart, color=['species'], legend_loc='on data')
adata_biomart.obsm['X_pca'] = PCAEmbeddings_biomart.values

#scib package metric calculation for principal component regression and silhouette score per batch
pcr_comparison_just_biomart = scib.metrics.pcr_comparison(adata_unint, adata_biomart, covariate='sample', n_comps=30)
silhouetteBatch_just_biomart = scib.metrics.silhouette_batch(adata_biomart, batch_key="sample", label_key="cell_type", embed="X_umap", return_all=True)

#calculate graph connectivity score
sc.pp.neighbors(adata_biomart, n_pcs=10)
graph_connectivity_just_biomart = scib.metrics.graph_connectivity(adata_biomart, label_key="cell_type")

#calculate average of k-nearest neighbor batch effect on downsampled score, like proposed by Büttner et al. (2019) for big datasets
subset_size = 0.1
num_cells_to_subsample = int(subset_size * adata_biomart.n_obs)
subset_indices = np.random.choice(adata_biomart.n_obs, num_cells_to_subsample, replace=False)
subsampled_adata = adata_biomart[subset_indices, :]
kBET_just_biomart = scib.metrics.kBET(subsampled_adata, batch_key="sample", label_key="cell_type", type_="embed", embed="X_umap")

LISI_score_just_biomart = scib.metrics.ilisi_graph(adata_biomart, batch_key="sample", type_="embed", use_rep="X_umap", n_cores=30)
cell_cycle_biomart = scib.metrics.cell_cycle(adata_unint,adata_biomart, batch_key="sample", organism='human')
clisi_graph_biomart = scib.metrics.clisi_graph(adata_biomart, label_key='cell_type', type_='embed', use_rep='X_umap',n_cores=40)
#subsample, to make calculation possible with big datasets
isolated_labels_f1_biomart = scib.metrics.isolated_labels_f1(subsampled_adata, label_key="cell_type", batch_key="sample", embed='X_umap', iso_threshold=50)

#Build Dataframe with results for biomart
output_biomart = pd.DataFrame()
output_biomart.loc['graph_connectivity', 'just_biomart'] = graph_connectivity_just_biomart
output_biomart.loc['pcr_comparison', 'just_biomart'] = pcr_comparison_just_biomart
output_biomart.loc['silhouette_batch', 'just_biomart'] = silhouetteBatch_just_biomart[0]
output_biomart.loc['kBET', 'just_biomart'] = kBET_just_biomart
output_biomart.loc['LISI', 'just_biomart'] = LISI_score_just_biomart
output_biomart.loc['Cell cycle conservation', 'just_biomart'] = cell_cycle_biomart
output_biomart.loc['Cell-type LISI', 'just_biomart'] = clisi_graph_biomart
output_biomart.loc['isolated_labels_f1', 'just_biomart'] = isolated_labels_f1_biomart
output_biomart.loc['Species mixing score', 'just_biomart'] = (graph_connectivity_just_biomart + pcr_comparison_just_biomart + silhouetteBatch_just_biomart[0] + kBET_just_biomart + LISI_score_just_biomart)/5 # Species Score
output_biomart.T.to_csv('/media/Storage/R/just_biomart/metrics_just_biomart_cell_type.csv')

'CALCULATE METRICS FOR OMA'

# Read in sce object, convert to anndata obj
adata_OMA = sc.read_h5ad('/media/Storage/R/Revision_scib/Conversion/SeuratObject.oma_10_08_22.h5ad')
# Read in UMAP embeddings
UMAPEmbeddings_OMA = pd.read_csv('/media/Storage/R/OMA/UMAPEmbedings.csv', skiprows=0, header=0 ,index_col=0)
'Remake the Seurat Clustering since its not saved in the SCE object'
PCAEmbeddings_OMA = pd.read_csv('/media/Storage/R/OMA/PCAEmbedings.csv', skiprows=0, header=0, index_col=0)

# Check order
tmpcheck = UMAPEmbeddings_OMA.index == adata_OMA.obs_names
np.unique(tmpcheck, return_counts=True)
tmpcheck = PCAEmbeddings_OMA.index == adata_OMA.obs_names
np.unique(tmpcheck, return_counts=True)

#Assign UMAP parameters in anndata object
adata_OMA.obsm['X_umap'] = UMAPEmbeddings_OMA.values
sc.pl.umap(adata_OMA, color=['species'], legend_loc='on data')
adata_OMA.obsm['X_pca'] = PCAEmbeddings_OMA.values

#scib package metric calculation for principal component regression and silhouette score per batch
pcr_comparison_OMA = scib.metrics.pcr_comparison(adata_unint, adata_OMA, covariate='sample', n_comps=30)
silhouetteBatch_OMA = scib.metrics.silhouette_batch(adata_OMA, batch_key="sample", label_key="cell_type", embed="X_umap", return_all=True)

#calculate graph connectivity score
sc.pp.neighbors(adata_OMA, n_pcs=10)
graph_connectivity_OMA = scib.metrics.graph_connectivity(adata_OMA, label_key="cell_type")

#calculate average of k-nearest neighbor batch effect on downsampled score, like proposed by Büttner et al. (2019) for big datasets
subset_size = 0.1
num_cells_to_subsample = int(subset_size * adata_OMA.n_obs)
subset_indices = np.random.choice(adata_OMA.n_obs, num_cells_to_subsample, replace=False)
subsampled_adata = adata_OMA[subset_indices, :]
kBET_OMA = scib.metrics.kBET(subsampled_adata, batch_key="sample", label_key="cell_type", type_="embed", embed="X_umap")

LISI_score_OMA = scib.metrics.ilisi_graph(adata_OMA, batch_key="sample", type_="embed", use_rep="X_umap", n_cores=30)
cell_cycle_OMA= scib.metrics.cell_cycle(adata_unint,adata_OMA, batch_key="sample", organism='human')
clisi_graph_OMA = scib.metrics.clisi_graph(adata_OMA, label_key='cell_type', type_='embed', use_rep='X_umap',n_cores=40)
isolated_labels_f1_OMA = scib.metrics.isolated_labels_f1(subsampled_adata, label_key="cell_type", batch_key="sample", embed='X_umap', iso_threshold=50)

#Build Dataframe with results for biomart
output_OMA = pd.DataFrame()
output_OMA.loc['graph_connectivity', 'OMA'] = graph_connectivity_OMA
output_OMA.loc['pcr_comparison', 'OMA'] = pcr_comparison_OMA
output_OMA.loc['silhouette_batch', 'OMA'] = silhouetteBatch_OMA[0]
output_OMA.loc['kBET', 'OMA'] = kBET_OMA
output_OMA.loc['LISI', 'OMA'] = LISI_score_OMA
output_OMA.loc['Cell cycle conservation', 'OMA'] = cell_cycle_OMA
output_OMA.loc['Cell-type LISI', 'OMA'] = clisi_graph_OMA
output_OMA.loc['isolated_labels_f1', 'OMA'] = isolated_labels_f1_OMA
output_OMA.loc['Species mixing score', 'OMA'] = (graph_connectivity_OMA + pcr_comparison_OMA + silhouetteBatch_OMA[0] + kBET_OMA + LISI_score_OMA)/5 # Species Score
output_OMA.T.to_csv('/media/Storage/R/OMA/metrics_OMA_cell_type.csv')

output_ALL = pd.DataFrame()
output_ALL = output_ALL.append(output_ortho.T, ignore_index=True)
output_ALL = output_ALL.append(output_OMA.T, ignore_index=True)
output_ALL = output_ALL.append(output_biomart.T, ignore_index=True)
output_ALL = output_ALL.append(output_para.T, ignore_index=True)
output_ALL = output_ALL.set_index([pd.Index(['OrthoIntegrate', 'OMA', 'biomart','InParanoid'])])
output_ALL.to_csv('/media/Storage/R/Revision_scib/metric_comb_cell_type.csv')

metric_comb =  pd.read_csv('/media/Storage/R/Revision_scib/metric_comb_cell_type.csv', index_col=0)
#Calculate Bio Conservation Score Ortho
OrthoIntegrate_Bio_Score = metric_comb.loc['OrthoIntegrate']
OrthoIntegrate_Bio_Score = (OrthoIntegrate_Bio_Score['Cell cycle conservation'] + OrthoIntegrate_Bio_Score['Species-type LISI'] + OrthoIntegrate_Bio_Score['isolated_labels_f1'] + OrthoIntegrate_Bio_Score['Silhouette Coefficient'])/4
OrthoIntegrate_Total_Score = (OrthoIntegrate_Bio_Score * 0.6 + output_ortho.loc['Species mixing score'] * 0.4)
#Calculate Bio Conservation Score OMA
OMA_Bio_Score = metric_comb.loc['OMA']
OMA_Bio_Score = (OMA_Bio_Score['Cell cycle conservation'] + OMA_Bio_Score['Species-type LISI'] + OMA_Bio_Score['isolated_labels_f1'] + OMA_Bio_Score['Silhouette Coefficient'])/4
OMA_Total_Score = (OMA_Bio_Score * 0.6 + output_OMA.loc['Species mixing score'] * 0.4)
#Calculate Bio Conservation Score biomart
biomart_Bio_Score = metric_comb.loc['biomart']
biomart_Bio_Score = (biomart_Bio_Score['Cell cycle conservation'] + biomart_Bio_Score['Species-type LISI'] + biomart_Bio_Score['isolated_labels_f1'] + biomart_Bio_Score['Silhouette Coefficient'])/4
biomart_Total_Score = (biomart_Bio_Score * 0.6 + output_biomart.loc['Species mixing score'] * 0.4)
#Calculate Bio Conservation Score InParanoid
InParanoid_Bio_Score = metric_comb.loc['InParanoid']
InParanoid_Bio_Score = (InParanoid_Bio_Score['Cell cycle conservation'] + InParanoid_Bio_Score['Species-type LISI'] + InParanoid_Bio_Score['isolated_labels_f1'] + InParanoid_Bio_Score['Silhouette Coefficient'])/4
InParanoid_Total_Score = (InParanoid_Bio_Score * 0.6 + output_para.loc['Species mixing score'] * 0.4)

#Add column
metric_comb['Bio Conservation Score'] = [OrthoIntegrate_Bio_Score, OMA_Bio_Score, biomart_Bio_Score, InParanoid_Bio_Score]
metric_comb['Total Score'] = [OrthoIntegrate_Total_Score.values[0], OMA_Total_Score.values[0], biomart_Total_Score.values[0], InParanoid_Total_Score.values[0]]
metric_comb.to_csv('/media/Storage/R/Revision_scib/metric_comb_cell_type.csv')

