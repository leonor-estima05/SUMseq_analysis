from pathlib import Path
from collections import Counter
import logging
import pandas as pd
import snapatac2 as snap

for name in ("kaleido", "choreographer"):
  logging.getLogger(name).setLevel(logging.WARNING)

#1 Initializing 
Cluster_dir = Path("/home/projects/icell/sumseq_pilot/ATAC/snapatac2_out/clustering")
Out_dir = Cluster_dir.parent / "peaks_motifs"
Plot_dir = Out_dir / "plots"
Out_dir.mkdir(parents=True, exist_ok=True)
Plot_dir.mkdir(parents=True, exist_ok=True)

Genome = snap.genome.hg38

#The hiPSC village samples will be pooled together for peak calling since they are all villages and this will allow for better peaks

Groups = {
  "iMono": "iMono",
  "cMo": "cMo",
  "ncMo": "ncMo",
  "MV1": "hiPSC", "MV2": "hiPSC", "MV3": "hiPSC", "MV5": "hiPSC", "MV6": "hiPSC", "Ctrl3": "hiPSC",
}

Marker_pvalue = 0.01 #same value as used in single cell best practices

Motif_max_fdr = 1e-5
#Only show motifs below this FDR, in the book they use 1e-20 but I need to use a more lenient cutoff since this data is more low depth than theirs

#2 Loading Clustered Dataset and Assigning Groups
data = snap.read_dataset(str(Cluster_dir / "combined.h5ads"))
samples = list(data.obs["sample"])


labels = [Groups[s] for s in samples]
print("Cells per group:", Counter(labels))

groups_to_call = sorted(set(labels))

#3 Peak Calling with MACS3
snap.tl.macs3(data, groupby=labels, selections=set(groups_to_call))
peaks = snap.tl.merge_peaks(data.uns["macs3"], Genome.chrom_sizes)
print(f"Merged peak set: {peaks.shape[0]} peaks")

#Cell and peak matrix
peak_matrix = snap.pp.make_peak_matrix(data, use_rep=peaks["Peaks"])
data.close()

peak_matrix.obs["group"] = labels
peak_matrix = peak_matrix[peak_matrix.obs["group"]].copy()
peak_matrix.write(Out_dir / "peak_matrix.h5ad")
print(peak_matrix)

#4 Group Specific Marker Peaks (comparing accessibility between groups, this is not a significance test)
markers = snap.tl.marker_regions(peak_matrix, groupby="group", pvalue=Marker_pvalue)
for g, regions in markers.items():
    print(f"{g}: {len(regions)} marker peaks")
    with open(Out_dir / f"marker_peaks_{g}.txt", "w") as f:
        f.write("\n".join(regions) + "\n")
 
snap.pl.regions(peak_matrix, groupby="group", peaks=markers, interactive=False,
                out_file=str(Plot_dir / "marker_peaks_heatmap.png"))

#5 Motif Enrichment 
motifs = snap.datasets.cis_bp(unique=True)
enrichment = snap.tl.motif_enrichment(
    motifs=motifs,
    regions=markers,
    genome_fasta=Genome,
)
 
for g, df in enrichment.items():
    out = Out_dir / f"motif_enrichment_{g}.csv"
    try:
        df.write_csv(str(out))      
    except AttributeError:
        df.to_csv(out)              
 
snap.pl.motif_enrichment(enrichment, max_fdr=Motif_max_fdr, height=1200, interactive=False,
                         out_file=str(Plot_dir / "motif_enrichment.png"))
 
print("Finished with no errors")
