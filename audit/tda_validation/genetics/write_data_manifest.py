"""Write data_manifest.json: source URL, sha256, size of every raw/derived data file used (data live outside git).
Command: python write_data_manifest.py"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tda_common import DATA_ROOT, sha256

GEO = "https://ftp.ncbi.nlm.nih.gov/geo"
SRC = {
    "buettner/G1_singlecells_counts.txt": "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/805/E-MTAB-2805/Files/G1_singlecells_counts.txt",
    "buettner/S_singlecells_counts.txt": "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/805/E-MTAB-2805/Files/S_singlecells_counts.txt",
    "buettner/G2M_singlecells_counts.txt": "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/805/E-MTAB-2805/Files/G2M_singlecells_counts.txt",
    "buettner/E-MTAB-2805.sdrf.txt": "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/805/E-MTAB-2805/Files/E-MTAB-2805.sdrf.txt",
    "genelists/Mus_musculus.csv": "https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Mus_musculus.csv",
    "genelists/Homo_sapiens.csv": "https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Homo_sapiens.csv",
    "u2os/GSE146773_Counts.csv.gz": GEO + "/series/GSE146nnn/GSE146773/suppl/GSE146773_Counts.csv.gz",
    "u2os/GSE146773_fucci_coords.csv.gz": GEO + "/series/GSE146nnn/GSE146773/suppl/GSE146773_fucci_coords.csv.gz",
    "caulobacter/GSM1120445_Laublab_BglII_HiC_NA1000_swarmer_cell_untreated_replicate1_overlap_after_normalization.txt.gz":
        GEO + "/samples/GSM1120nnn/GSM1120445/suppl/GSM1120445_Laublab_BglII_HiC_NA1000_swarmer_cell_untreated_replicate1_overlap_after_normalization.txt.gz",
    "ecoli/GSM2870407_mat_BC70_TACT_wt_MM_30C.txt.gz": GEO + "/samples/GSM2870nnn/GSM2870407/suppl/GSM2870407_mat_BC70_TACT_wt_MM_30C.txt.gz",
    "gm12878/gm12878_insitu_combined_chr1_145000000_249250621_250kb_KR_observed.npy":
        "DERIVED by fetch_gm12878_chr1q.py: hicstraw.straw('observed','KR','https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic','1:145000000:249250621','1:145000000:249250621','BP',250000) (Rao et al. 2014, GSE63525); the .hic file was read remotely, not downloaded",
    "influenza/genomeset.dat.gz": "https://ftp.ncbi.nih.gov/genomes/INFLUENZA/genomeset.dat.gz",
    "influenza/influenza_na.dat.gz": "https://ftp.ncbi.nih.gov/genomes/INFLUENZA/influenza_na.dat.gz",
    "influenza/influenza.fna.gz": "https://ftp.ncbi.nih.gov/genomes/INFLUENZA/influenza.fna.gz",
    "influenza/README": "https://ftp.ncbi.nih.gov/genomes/INFLUENZA/README",
    "rcsb_reference/1CRN.pdb": "https://files.rcsb.org/download/1CRN.pdb",
    "rcsb_reference/1UBQ.pdb": "https://files.rcsb.org/download/1UBQ.pdb",
    "rcsb_reference/4OBE.pdb": "https://files.rcsb.org/download/4OBE.pdb",
}
for s in range(1, 9):
    SRC[f"influenza/work/aln_seg{s}.npy"] = "DERIVED by test3_influenza.py prep (FAMSA via pyfamsa 0.7.0, default parameters, 300 genomes)"
SRC["influenza/work/prep_info.json"] = "DERIVED by test3_influenza.py prep (strain names, accessions, subtypes of the 300 genomes)"
out = {"data_root": DATA_ROOT, "files": {}, "failed_or_unused_urls": {
    "https://raw.githubusercontent.com/scverse/scanpy_usage/master/180209_cell_cycle/data/regev_lab_cell_cycle_genes.txt": "HTTP 404 (1 try); not needed: the tinyatlas lists were used instead",
    "https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/README.md": "HTTP 404 (1 try); documentation only"}}
for rel, src in SRC.items():
    p = os.path.join(DATA_ROOT, rel)
    out["files"][rel] = {"source": src, "sha256": sha256(p), "size_bytes": os.path.getsize(p)} if os.path.exists(p) else {"source": src, "status": "MISSING"}
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "data_manifest.json"), "w"), indent=1)
print(sum(1 for v in out["files"].values() if "sha256" in v), "files hashed;", [k for k, v in out["files"].items() if "status" in v])
