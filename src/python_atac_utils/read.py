from os.path import join, isfile, isdir
from typing import Literal
import logging

import anndata as ad
import pandas as pd
import scanpy as sc


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

    if isfile(join(dir, 'peak_annotation.tsv')):
        annot = pd.read_csv(join(dir, 'peak_annotation.tsv'), sep='\t')
        annot.index = annot.chrom + ':' + annot.start.astype(str) + '-' + annot.end.astype(str)
        annot.index.name = 'peak'

        adata.uns['peak_annotation'] = annot
    else:
        logging.warn('peak_annotation.tsv not found in the supplied dir')

    if isfile(join(dir, 'fragments.tsv.gz')):
        adata.uns['fragments'] = join(dir, 'fragments.tsv.gz')
    else:
        logging.warn('fragments.tsv.gz not found in the supplied dir')

    adata.obs = obs
    adata.var = var

    return adata
