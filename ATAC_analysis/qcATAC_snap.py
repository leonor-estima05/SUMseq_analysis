from pathlib import Path
import pandas as pd 
import snapatac2 as snap

#1 Initializing 

Frag_dir = "/home/projects/icell/sumseq_pilot/ATAC/ATAC_outdir/ATAC/3.fragmentDir"

Samples = {
    "cMo":   f"{Frag_dir}/cMo_S8_fragments.bed.gz",
    "Ctrl3": f"{Frag_dir}/Ctrl3_S6_fragments.bed.gz",
    "iMono": f"{Frag_dir}/iMono_S7_fragments.bed.gz",
    "intMo": f"{Frag_dir}/intMo_S10_fragments.bed.gz",
    "MV1":   f"{Frag_dir}/MV1_S1_fragments.bed.gz",
    "MV2":   f"{Frag_dir}/MV2_S2_fragments.bed.gz",
    "MV3":   f"{Frag_dir}/MV3_S3_fragments.bed.gz",
    "MV5":   f"{Frag_dir}/MV5_S4_fragments.bed.gz",
    "MV6":   f"{Frag_dir}/MV6_S5_fragments.bed.gz",
    "ncMo":  f"{Frag_dir}/ncMo_S9_fragments.bed.gz",
}

Genome = snap.genome.hg38

#First run of this script will be to just see the basic QC graphs and counts with NO filtering, second run will be with filtering according to graphs from first run
Run_Filtering = False

Whitelists = {}

#Thresholds for after run 1
Thresholds = {
  "cMo": {"min_counts": 1000, "max_counts": 100000, "min_tsse": 3}, #lenient mintsse for initial filtering due to low library depth for primary monos
  "Ctrl3": {"min_counts": 1000, "max_counts": 100000, "min_tsse": 3},
  "iMono": {"min_counts": 1000, "max_counts": 100000, "min_tsse": 3},
  "intMo": {"min_counts": 1000, "max_counts": 100000, "min_tsse": 3},
  "MV1": {"min_counts": 1000, "max_counts": 100000, "min_tsse": 3},
  "MV2": {"min_counts": 1000, "max_counts": 100000, "min_tsse": 3},
  "MV3": {"min_counts": 1000, "max_counts": 100000, "min_tsse": 3},
  "MV5": {"min_counts": 1000, "max_counts": 100000, "min_tsse": 3},
  "MV6": {"min_counts": 1000, "max_counts": 100000, "min_tsse": 3},
  "ncMo":{"min_counts": 1000, "max_counts": 100000, "min_tsse": 3},
}

Out_dir = Path("/home/projects/icell/sumseq_pilot/ATAC/snapatac2_out")
Plot_dir = Out_dir / "qc_plots"
Out_dir.mkdir(parents=True, exist_ok=True)
Plot_dir.mkdir(parents=True, exist_ok=True)

#2 Extra

def load_whitelist(path):
    """Read a one-barcode-per-line file into a list, or return None."""
    if path is None:
        return None
    with open(path) as f:
        return [line.strip() for line in f if line.strip()]
 
 
def qc_sample(name, fragment_file):
    print(f"\n=== {name} ===")
 
    #Importing
    adata = snap.pp.import_fragments(
        fragment_file,
        chrom_sizes=Genome,
        sorted_by_barcode=False,
        min_num_fragments=200,
        whitelist=load_whitelist(Whitelists.get(name)),
    )
    print(f"Imported {adata.n_obs} barcodes")
 
    #TSS enrichment (different calculation than ArchR)
    snap.metrics.tsse(adata, Genome)
 
    #Plotting for Run 1
    snap.pl.frag_size_distr(adata, interactive=False,
                            out_file=str(Plot_dir / f"{name}_frag_size.png"))
    #Tsse vs number of fragments (good quality cells sit on the top right)
    snap.pl.tsse(adata, interactive=False,
                 out_file=str(Plot_dir / f"{name}_tsse.png"))
    #Per-cell QC table (n_fragment, frac_dup, frac_mito, tsse)
    adata.obs.to_csv(Out_dir / f"{name}_qc_metrics_prefilter.csv")
    print(adata.obs[["n_fragment", "tsse"]].describe().round(2))
 
    if not Run_Filtering:
        print("Run_Filtering = False -> stopping after plots.")
        #Pre-filter summary so all 10 samples can be compared side by side
        return {"sample": name, "barcodes": adata.n_obs,
                "median_fragments": adata.obs["n_fragment"].median(),
                "median_tsse": adata.obs["tsse"].median(),
                "median_frac_mito": adata.obs["frac_mito"].median()}
 
    #3 Filtering
    t = Thresholds[name]
    n_before = adata.n_obs
    snap.pp.filter_cells(adata,
                         min_counts=t["min_counts"],
                         max_counts=t["max_counts"],
                         min_tsse=t["min_tsse"])
    print(f"Cell filter: {n_before} -> {adata.n_obs} cells")
 
    #4 Doublet Removal Using Scrublet (may fail on samples with very few cells, e.g. low-depth primary monos)
    snap.pp.add_tile_matrix(adata)
    snap.pp.select_features(adata)
    try:
        snap.pp.scrublet(adata)
        n_before = adata.n_obs
        snap.pp.filter_doublets(adata)
        print(f"Doublet filter: {n_before} -> {adata.n_obs} cells")
    except Exception as e:
        print(f"Scrublet skipped for {name}: {e}")
 
    #5 Save QC Run
    adata.obs["sample"] = name
    adata.obs.to_csv(Out_dir / f"{name}_qc_metrics_postfilter.csv")
    adata.write(Out_dir / f"{name}_qc.h5ad")
    return {"sample": name, "cells_kept": adata.n_obs,
            "median_fragments": adata.obs["n_fragment"].median(),
            "median_tsse": adata.obs["tsse"].median()}

#6 Run
if __name__ == "__main__":
    summary = []
    for name, frag in Samples.items():
        result = qc_sample(name, frag)
        if result:
            summary.append(result)
 
    if summary:
        summary = pd.DataFrame(summary)
        stage = "postfilter" if Run_Filtering else "prefilter"
        summary.to_csv(Out_dir / f"qc_summary_{stage}.csv", index=False)
        print("\n", summary.to_string(index=False))


  
