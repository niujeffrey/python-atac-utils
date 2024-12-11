from os.path import join, isfile, isdir
from typing import Literal
import logging

import anndata as ad
import pandas as pd
import scanpy as sc
import muon.atac as ac


def read_10x_filtered_peak_matrix(
    dir: str,
) -> ad.AnnData:
    """Produce a matrix given by 10X cellranger

    Parameters
    ----------
    dir
        The directory containing filtered cellranger output (ends with outs/)
    """
    if not isdir(join(dir, 'filtered_peak_bc_matrix')):
        raise ValueError('Directory `filtered_peak_bc_matrix` not found in the supplied dir')

    adata = sc.read_mtx(join(dir, 'filtered_peak_bc_matrix/matrix.mtx')).T

    obs = pd.read_csv(
        join(dir, 'filtered_peak_bc_matrix/barcodes.tsv'),
        header=None,
        index_col=0,
        names=['barcode']
    )

    if isfile(join(dir, 'singlecell.csv')):
        metadata = pd.read_csv(join(dir, 'singlecell.csv'), index_col=0)
        metadata = metadata.iloc[1:] # Remove the extraneous row
    else:
        logging.warn('singlecell.csv not found in the supplied dir')

    obs = obs.join(metadata)

    var = pd.read_csv(
        join(dir, 'filtered_peak_bc_matrix/peaks.bed'),
        header=None,
        sep='\t',
        names=['chrom', 'start', 'end']
    )
    var.index = var.chrom + ':' + var.start.astype(str) + '-' + var.end.astype(str)

    adata.obs = obs
    adata.var = var

    annot_path = join(dir, "peak_annotation.tsv")
    if isfile(annot_path):
        # Copied from muon
        try:
            ac.tl.add_peak_annotation(adata, annot_path)
            print(
                f"Added peak annotation from {annot_path} to .uns['atac']['peak_annotation']"
            )
        except AttributeError:
            logging.warn(
                f"Peak annotation from {annot_path} could not be added. Please check the annotation file is formatted correctly."
            )

    frag_path = join(dir, "fragments.tsv.gz")
    if isfile(frag_path):
        print(f"Located fragments file: {frag_path}")
        try:
            ac.tl.locate_fragments(adata, frag_path)
        except ImportError:
            logging.warn(
                "Pysam is not installed. To work with the fragments file please install pysam (pip install pysam)."
            )
            if "files" not in adata.uns:
                adata.uns["files"] = dict()
            adata.uns["files"]["fragments"] = frag_path

    return adata
