import pandas as pd
import scanpy as sc

from mcDETECT.utils import *
from mcDETECT.model import *

import warnings
warnings.filterwarnings("ignore")
sc.settings.verbosity = 0

for dataset in ["MERSCOPE_9m_AD_R1", "MERSCOPE_9m_AD_R2", "MERSCOPE_9m_WT_R1", "MERSCOPE_9m_WT_R2"]:
    
    # File paths
    data_path = f"../data/{dataset}/"
    output_path = f"../output/{dataset}/"
    
    # ==================== Read data ==================== #
    
    # Transcripts
    transcripts = pd.read_parquet(data_path + "processed_data/transcripts.parquet")

    # Genes
    genes = pd.read_csv(data_path + "processed_data/genes.csv")
    genes = list(genes.iloc[:, 0])

    # Negative control markers
    nc_genes = pd.read_csv(data_path + "processed_data/negative_controls.csv")
    nc_genes = list(nc_genes["Gene"])

    # Spots
    spots = sc.read_h5ad(data_path + "intermediate_data/spots_raw.h5ad")

    # Markers
    syn_genes = ["Camk2a", "Cplx2", "Slc17a7", "Ddn", "Syp", "Map1a", "Shank1", "Syn1", "Gria1", "Gria2", "Cyfip2", "Vamp2", "Bsn", "Slc32a1", "Nfasc", "Syt1", "Tubb3", "Nav1", "Shank3", "Mapt"]
    print("Number of syn_genes: ", len(syn_genes))

    # ==================== Fine detection ==================== #

    print(f"Running fine detection for {dataset}...")

    mc = mcDETECT(type = "discrete", transcripts = transcripts, gnl_genes = syn_genes, nc_genes = nc_genes, eps = 1.5,
                minspl = 3, grid_len = 1, cutoff_prob = 0.95, alpha = 10, low_bound = 3, size_thr = 4.0,
                in_soma_thr = 0.1, l = 1, rho = 0.2, s = 1, nc_top = 20, nc_thr = 0.1)
    granules = mc.detect()
    granules.to_parquet(output_path + "granules.parquet")
    print("Granules shape: ", granules.shape)