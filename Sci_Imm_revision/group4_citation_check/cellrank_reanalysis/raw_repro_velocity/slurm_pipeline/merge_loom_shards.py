#!/usr/bin/env python3
"""Merge chromosome-sharded velocyto loom outputs into one loom matching
standard velocyto run10x output conventions (CellID = '<capture>:<barcode>x').
"""
import sys
import gzip
import numpy as np
import loompy

def load_barcodes(barcodes_gz):
    with gzip.open(barcodes_gz, "rt") as f:
        bcs = [line.strip().split("-")[0] for line in f]
    return bcs

def merge(capture_name, barcodes_gz, shard_loom_paths, out_path):
    master_bcs = load_barcodes(barcodes_gz)
    master_index = {bc: i for i, bc in enumerate(master_bcs)}
    n_cells = len(master_bcs)

    all_gene_ids, all_gene_rows, layer_names = [], [], None
    layer_mats = {}

    for path in shard_loom_paths:
        with loompy.connect(path, validate=False) as ds:
            if layer_names is None:
                layer_names = [l for l in ds.layers.keys() if l != ""]
                for l in layer_names:
                    layer_mats[l] = []

            shard_bcs = [cid.split(":")[1].rstrip("x") for cid in ds.ca["CellID"]]
            col_map = np.full(ds.shape[1], -1, dtype=np.int64)
            for j, bc in enumerate(shard_bcs):
                if bc in master_index:
                    col_map[j] = master_index[bc]

            n_genes = ds.shape[0]
            all_gene_ids.append(ds.ra["Accession"].copy())
            all_gene_rows.append({k: ds.ra[k].copy() for k in ds.ra.keys()})

            for l in layer_names:
                dense = ds.layers[l][:, :]
                out = np.zeros((n_genes, n_cells), dtype=dense.dtype)
                valid = col_map >= 0
                out[:, col_map[valid]] = dense[:, valid]
                layer_mats[l].append(out)

    merged_layers = {l: np.concatenate(layer_mats[l], axis=0) for l in layer_names}
    ra = {k: np.concatenate([g[k] for g in all_gene_rows]) for k in all_gene_rows[0].keys()}
    ca = {"CellID": np.array([f"{capture_name}:{bc}x" for bc in master_bcs])}

    main_layer = merged_layers.pop("spliced")
    loompy.create(out_path, main_layer, ra, ca)
    with loompy.connect(out_path) as ds:
        ds.layers["spliced"] = main_layer
        for l, mat in merged_layers.items():
            ds.layers[l] = mat

    print(f"Merged {len(shard_loom_paths)} shards -> {out_path}")
    print(f"  shape: {n_cells} cells x {sum(len(g) for g in all_gene_ids)} genes")

if __name__ == "__main__":
    capture_name = sys.argv[1]
    barcodes_gz = sys.argv[2]
    out_path = sys.argv[3]
    shard_paths = sys.argv[4:]
    merge(capture_name, barcodes_gz, shard_paths, out_path)
