from pathlib import Path
from collections import Counter
import snapatac2 as snap

#This combines the QC'd samples (9 with intMo excluded) and then does feature selection and spectral embedding

#1 Initializing
QC_dir = Path("/home/projects/icell/sumseq_pilot/ATAC/snapatac2_out")
Out_dir = QC_dir / "clustering"
Plot_dir = Out_dir / "plots"
Out_dir.mkdir(parents=True, exist_ok=True)
Plot_dir.mkdir(parents=True, exist_ok=True)

#Loading the samples
Samples = ["cMo", "Ctrl3", "iMono", "MV1", "MV2", "MV3", "MV5", "MV6", "ncMo"]

#snapatac2 uses 250,000 number of tiles kept as features
N_features = 250000

#Batch Correction NOT included, each sample is a different cell type/from a different donor so it does not make sense

#Leiden Resolution set to 0.5 initially, can be adjusted
Resolution = 0.5

#2 Combining Samples into AnnData Object 
adatas = [(name, snap.read(QC_dir / f"{name}_qc.h5ad")) for name in Samples]
data = snap.AnnDataSet(adatas=adatas, filename=str(Out_dir / "combined.h5ads"))
print(data)

data.obs_names = [f"{s}:{bc}" for s, bc in zip(data.obs["sample"], data.obs_names)]

#Bringing QC metrics
for col in ["tsse", "n_fragment"]:
    try:
        data.obs[col] = data.adatas.obs[col]
    except Exception as e:
        print(f"Could not copy {col} to combined obs: {e}")
 
print("Cells per sample:", Counter(list(data.obs["sample"])))

#3 Dimensionality Reduction
snap.pp.select_features(data, n_features=N_features)
snap.tl.spectral(data)        
use_rep = "X_spectral"

#4 Clustering and UMAP

snap.pl.umap(data, color="sample", interactive=False,
             out_file=str(Plot_dir / "umap_sample.png"))
#Colour by cluster
snap.pl.umap(data, color="leiden", interactive=False,
             out_file=str(Plot_dir / "umap_leiden.png"))

#I am coloring the map by QC (so lower quality cells are distinguished)
for col in ["tsse", "n_fragment"]:
    try:
        snap.pl.umap(data, color=col, interactive=False,
                     out_file=str(Plot_dir / f"umap_{col}.png"))
    except Exception as e:
        print(f"Skipped UMAP coloured by {col}: {e}")

#5 Saving

obs = data.obs
try:
    obs.to_csv(Out_dir / "cells_clusters.csv")
except AttributeError:
    obs.write_csv(str(Out_dir / "cells_clusters.csv"))   #in case obs comes back as polars
 
data.close()
