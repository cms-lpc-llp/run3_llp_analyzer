#!/usr/bin/env python3
"""
ROOT to HDF5 Converter for DGCNN Inference
==========================================

Converts ROOT files from CMS ntuples to HDF5 format compatible with DGCNN inference.

This script processes CSC rechit cluster data from ROOT files and outputs HDF5 files
that can be directly used with the inference pipeline.

Usage:
    # Convert a single ROOT file
    python convert_root_to_h5.py --input data.root --output data.h5
    
    # Convert a directory of ROOT files
    python convert_root_to_h5.py --input /path/to/root/files --output /path/to/output
    
    # Process as signal (LLP) data
    python convert_root_to_h5.py --input data.root --output data.h5 --label signal
    
    # Process as background data
    python convert_root_to_h5.py --input data.root --output data.h5 --label background

For help:
    python convert_root_to_h5.py --help
"""

import argparse
import numpy as np
import h5py
import os
import re
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Dict, Tuple, Optional, Any
from tqdm import tqdm

try:
    import uproot
    import awkward as ak
except ImportError:
    print("Error: This script requires 'uproot' and 'awkward' packages.")
    print("Install with: pip install uproot awkward")
    exit(1)


# Label mapping
LABEL_MAP = {
    'signal': 0,      # LLP signal (ggH, etc.)
    'background': 1,  # Background (QCD, etc.)
    'ggH': 0,         # Alias for signal
    'QCD': 1,         # Alias for background
    'Gun': 2,         # Gun samples (if used)
}


def get_root_files(directory: str) -> List[str]:
    """Recursively search for .root files in a directory."""
    root_files = []
    for dirpath, _, filenames in os.walk(directory):
        for filename in filenames:
            if filename.endswith('.root'):
                root_files.append(os.path.join(dirpath, filename))
    return sorted(root_files)


def flatten_rechit_cluster_data(cluster_ids, data_column, label: str, cluster_is_llp, size_cutoff: int):
    """
    Flatten rechit data by cluster, optionally filtering by LLP status.
    
    Args:
        cluster_ids: Cluster ID assignments for each rechit
        data_column: Data column to flatten
        label: Label category (determines if LLP filtering is applied)
        cluster_is_llp: Boolean mask for LLP clusters (only used for signal)
        size_cutoff: Minimum cluster size to keep
    """
    filtered_data = []
    max_cluster = ak.max(cluster_ids)
    if max_cluster is None or max_cluster < 0:
        return ak.Array([])
    
    for i in range(int(max_cluster) + 1):
        # Pull out all rechits belonging to this cluster
        altered_rechits = ak.unflatten(data_column[cluster_ids == i], counts=1)
        
        # For signal (LLP) data, filter out clusters not matched to LLPs
        if LABEL_MAP.get(label, 1) == 0 and cluster_is_llp is not None:
            altered_rechits = altered_rechits[cluster_is_llp[:, i:i+1]]
        
        filtered_data.append(altered_rechits)
    
    if not filtered_data:
        return ak.Array([])
    
    filtered_data = ak.concatenate(filtered_data, axis=1)
    filtered_data = ak.flatten(filtered_data, axis=1)
    filtered_data = filtered_data[ak.num(filtered_data) > size_cutoff]
    return filtered_data


def flatten_cluster_data(cluster_data, size_column, label: str, cluster_is_llp, size_cutoff: int):
    """
    Flatten cluster-level data with size and LLP filtering.
    """
    if LABEL_MAP.get(label, 1) == 0 and cluster_is_llp is not None:
        truth = cluster_is_llp & (size_column > size_cutoff)
    else:
        truth = size_column > size_cutoff
    
    filtered_data = cluster_data[truth]
    flattened_data = ak.flatten(filtered_data, axis=1)
    return flattened_data


def process_csc_data(input_file: str, label: str, size_cutoff: int = 50) -> Tuple:
    """
    Process CSC rechit and cluster data from a ROOT file.
    
    Args:
        input_file: Path to ROOT file
        label: Data label ('signal' or 'background')
        size_cutoff: Minimum cluster size
        
    Returns:
        Tuple of (rechit_columns, cluster_columns, rechit_record, cluster_record)
    """
    with uproot.open(input_file) as f:
        if 'ntuples/llp' not in f:
            raise ValueError(f"Expected 'ntuples/llp' tree not found in {input_file}")
        tree = f['ntuples/llp']
    
    # Determine LLP matching for signal samples
    cluster_is_llp = None
    if LABEL_MAP.get(label, 1) == 0:
        try:
            gLLP_in_csc = ak.any(tree['gLLP_csc'].array(library='ak'), axis=1)
            cscRechitCluster_match_gParticle_minDeltaR = tree['cscRechitCluster_match_gParticle_minDeltaR'].array(library='ak')
            cscRechitCluster_match_gParticle_id = tree['cscRechitCluster_match_gParticle_id'].array(library='ak')
            
            # LLP cluster matching criteria:
            # - Event has gLLP in CSC
            # - Cluster is matched to particle with dR < 0.4
            # - Matched particle has LLP ID (9000006)
            cluster_is_llp = (
                gLLP_in_csc & 
                (cscRechitCluster_match_gParticle_minDeltaR < 0.4) & 
                (abs(cscRechitCluster_match_gParticle_id) == 9000006)
            )
        except Exception as e:
            print(f"Warning: Could not determine LLP matching for {input_file}: {e}")
            print("Processing as background instead.")
            label = 'background'
    
    # --- Process Rechits ---
    rechit_columns = [name for name in tree.keys() 
                      if name.startswith('cscRechit') and 'Cluster' not in name]
    useless_rechit_columns = ['cscRechitsE', 'cscRechitsChannels']
    rechit_columns = [col for col in rechit_columns if col not in useless_rechit_columns]
    
    data_rechit_columns = {col: tree[col].array(library='ak') for col in rechit_columns}
    cscRechitsClusterId = tree['cscRechitsClusterId'].array(library='ak')
    
    flattened_rechit_columns = {
        col: flatten_rechit_cluster_data(
            cscRechitsClusterId, data_rechit_columns[col], 
            label, cluster_is_llp, size_cutoff
        )
        for col in rechit_columns
    }
    
    rechit_record = ak.zip(flattened_rechit_columns)
    
    # --- Process Clusters ---
    cluster_columns = [name for name in tree.keys() 
                       if name.startswith('cscRechit') and 'Cluster' in name]
    useless_cluster_columns = [
        'cscRechitClusterGenMuonDeltaR', 'cscRechitClusterCaloJetVeto',
        'cscRechitClusterVertexR', 'cscRechitClusterVertexZ',
        'cscRechitClusterVertexDis', 'cscRechitClusterVertexChi2',
        'cscRechitClusterVertexN1', 'cscRechitClusterVertexN5',
        'cscRechitClusterVertexN10', 'cscRechitClusterVertexN15',
        'cscRechitClusterVertexN20', 'cscRechitClusterVertexN',
        'cscRechitCluster_match_gParticle_index',
        'cscRechitCluster_match_gParticle_minDeltaR',
        'cscRechitCluster_match_gParticle_id', 'cscRechitsClusterId'
    ]
    cluster_columns = [col for col in cluster_columns if col not in useless_cluster_columns]
    
    data_cluster_columns = {col: tree[col].array(library='ak') for col in cluster_columns}
    
    flattened_cluster_columns = {
        col: flatten_cluster_data(
            data_cluster_columns[col], 
            data_cluster_columns['cscRechitClusterSize'],
            label, cluster_is_llp, size_cutoff
        )
        for col in cluster_columns
    }
    
    # Apply timing cuts: keep clusters with time in [-5, 12.5] ns
    if len(flattened_cluster_columns.get('cscRechitClusterTimeWeighted', [])) > 0:
        mask = [
            -5 <= weighted <= 12.5 and -5 <= time <= 12.5
            for weighted, time in zip(
                flattened_cluster_columns['cscRechitClusterTimeWeighted'],
                flattened_cluster_columns['cscRechitClusterTime']
            )
        ]
        
        flattened_cluster_columns = {
            key: [val for val, include in zip(values, mask) if include]
            for key, values in flattened_cluster_columns.items()
        }
        rechit_record = [val for val, include in zip(rechit_record, mask) if include]
    
    cluster_record = ak.zip(flattened_cluster_columns)
    
    return rechit_columns, cluster_columns, rechit_record, cluster_record


def convert_to_inference_h5(
    input_file: str,
    output_file: str,
    label: str,
    size_cutoff: int = 50,
    include_labels: bool = True
) -> int:
    """
    Convert a ROOT file to HDF5 format for inference.
    
    This creates an HDF5 file compatible with the inference data loader,
    where each graph represents one cluster.
    
    Args:
        input_file: Path to input ROOT file
        output_file: Path to output HDF5 file
        label: Data label ('signal' or 'background')
        size_cutoff: Minimum cluster size
        include_labels: Whether to include ground truth labels
        
    Returns:
        Number of graphs (clusters) written
    """
    # Process ROOT file
    rechit_cols, cluster_cols, rechit_record, cluster_record = process_csc_data(
        input_file, label, size_cutoff
    )
    
    if len(rechit_record) == 0:
        print(f"Warning: No valid clusters found in {input_file}")
        return 0
    
    if len(rechit_record) != len(cluster_record):
        raise ValueError("Rechit and cluster records have different lengths")
    
    # Define which features to use for positions and node features
    # Position features: eta, phi, energy (for k-NN graph construction)
    pos_features = ['cscRechitsEta', 'cscRechitsPhi', 'cscRechitsE']
    pos_features = [f for f in pos_features if f in rechit_cols]
    if len(pos_features) < 3:
        # Fallback to first 3 rechit features
        pos_features = rechit_cols[:3]
    
    # All rechit features for node features
    x_features = rechit_cols
    
    # Extract file identifier for graph naming
    match = re.search(r'(\d+)(?:\.root)?$', os.path.basename(input_file))
    file_id = match.group(1) if match else '0'
    
    # Create output directory if needed
    os.makedirs(os.path.dirname(output_file) or '.', exist_ok=True)
    
    num_graphs = 0
    
    with h5py.File(output_file, 'w') as h5file:
        # Write metadata
        metadata = h5file.create_group('metadata')
        metadata.attrs['rechit feature names'] = str(x_features)
        metadata.attrs['position feature names'] = str(pos_features)
        metadata.attrs['cluster feature names'] = str(cluster_cols)
        metadata.attrs['source_file'] = input_file
        metadata.attrs['label'] = label
        
        # Process each cluster
        for idx, (rechit_cluster, cluster_info) in enumerate(zip(rechit_record, cluster_record)):
            # Extract position features
            pos_data = []
            for feat in pos_features:
                if feat in rechit_cols:
                    pos_data.append(np.asarray(rechit_cluster[feat], dtype=np.float32))
            
            if len(pos_data) == 0:
                continue
            
            pos = np.stack(pos_data, axis=1)  # [N, 3]
            
            # Extract all features for x
            x_data = []
            for feat in x_features:
                x_data.append(np.asarray(rechit_cluster[feat], dtype=np.float32))
            
            x = np.stack(x_data, axis=1)  # [N, F]
            
            # Skip small clusters
            if x.shape[0] < size_cutoff:
                continue
            
            # Create graph group
            graph_name = f'graph_{file_id}_{num_graphs}'
            graph_group = h5file.create_group(graph_name)
            
            # Store data
            graph_group.create_dataset('pos', data=pos)
            graph_group.create_dataset('x', data=x)
            
            if include_labels:
                y = LABEL_MAP.get(label, 1)
                graph_group.create_dataset('y', data=y)
            
            num_graphs += 1
        
        # Update metadata with graph count
        metadata.attrs['num_graphs'] = num_graphs
    
    return num_graphs


def process_single_file(args: Tuple) -> Tuple[str, int]:
    """Worker function for parallel processing."""
    input_file, output_file, label, size_cutoff = args
    try:
        num_graphs = convert_to_inference_h5(
            input_file, output_file, label, size_cutoff
        )
        return (input_file, num_graphs)
    except Exception as e:
        print(f"Error processing {input_file}: {e}")
        return (input_file, 0)


def convert_directory(
    input_dir: str,
    output_dir: str,
    label: str,
    size_cutoff: int = 50,
    num_workers: int = 4
) -> None:
    """
    Convert all ROOT files in a directory to HDF5 format.
    
    Args:
        input_dir: Directory containing ROOT files
        output_dir: Output directory for HDF5 files
        label: Data label for all files
        size_cutoff: Minimum cluster size
        num_workers: Number of parallel workers
    """
    root_files = get_root_files(input_dir)
    
    if not root_files:
        print(f"No ROOT files found in {input_dir}")
        return
    
    print(f"Found {len(root_files)} ROOT files")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Prepare arguments for parallel processing
    args_list = []
    for root_file in root_files:
        # Preserve relative directory structure
        rel_path = os.path.relpath(root_file, input_dir)
        output_file = os.path.join(output_dir, rel_path.replace('.root', '.h5'))
        os.makedirs(os.path.dirname(output_file) or '.', exist_ok=True)
        args_list.append((root_file, output_file, label, size_cutoff))
    
    # Process files
    total_graphs = 0
    if num_workers > 1:
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            futures = [executor.submit(process_single_file, args) for args in args_list]
            for future in tqdm(as_completed(futures), total=len(futures), desc="Converting files"):
                input_file, num_graphs = future.result()
                total_graphs += num_graphs
    else:
        for args in tqdm(args_list, desc="Converting files"):
            input_file, num_graphs = process_single_file(args)
            total_graphs += num_graphs
    
    print(f"\nConversion complete!")
    print(f"  Files processed: {len(root_files)}")
    print(f"  Total graphs (clusters): {total_graphs}")
    print(f"  Output directory: {output_dir}")


def main():
    parser = argparse.ArgumentParser(
        description='Convert ROOT files to HDF5 format for DGCNN inference',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Convert a single ROOT file (as background)
    python convert_root_to_h5.py --input data.root --output data.h5
    
    # Convert signal (LLP) data
    python convert_root_to_h5.py --input signal.root --output signal.h5 --label signal
    
    # Convert a directory of ROOT files
    python convert_root_to_h5.py --input /path/to/root/ --output /path/to/h5/
    
    # Use multiple workers for faster processing
    python convert_root_to_h5.py --input /path/to/root/ --output /path/to/h5/ --workers 8
        """
    )
    
    parser.add_argument(
        '--input', '-i',
        type=str,
        required=True,
        help='Input ROOT file or directory containing ROOT files'
    )
    parser.add_argument(
        '--output', '-o',
        type=str,
        required=True,
        help='Output HDF5 file or directory'
    )
    parser.add_argument(
        '--label', '-l',
        type=str,
        default='background',
        choices=['signal', 'background', 'ggH', 'QCD', 'Gun'],
        help='Data label (default: background). Use "signal" for LLP data.'
    )
    parser.add_argument(
        '--size-cutoff', '-s',
        type=int,
        default=50,
        help='Minimum cluster size to include (default: 50)'
    )
    parser.add_argument(
        '--workers', '-w',
        type=int,
        default=4,
        help='Number of parallel workers (default: 4)'
    )
    parser.add_argument(
        '--no-labels',
        action='store_true',
        help='Do not include ground truth labels in output'
    )
    
    args = parser.parse_args()
    
    input_path = Path(args.input)
    output_path = Path(args.output)
    
    if input_path.is_file():
        # Single file conversion
        if not input_path.suffix == '.root':
            print(f"Warning: Input file does not have .root extension")
        
        output_file = str(output_path)
        if output_path.is_dir():
            output_file = str(output_path / input_path.with_suffix('.h5').name)
        
        print(f"Converting {input_path} -> {output_file}")
        print(f"  Label: {args.label}")
        print(f"  Size cutoff: {args.size_cutoff}")
        
        num_graphs = convert_to_inference_h5(
            str(input_path),
            output_file,
            args.label,
            args.size_cutoff,
            include_labels=not args.no_labels
        )
        
        print(f"\nConversion complete!")
        print(f"  Graphs (clusters) written: {num_graphs}")
        print(f"  Output: {output_file}")
    
    elif input_path.is_dir():
        # Directory conversion
        convert_directory(
            str(input_path),
            str(output_path),
            args.label,
            args.size_cutoff,
            args.workers
        )
    
    else:
        print(f"Error: Input path not found: {input_path}")
        exit(1)


if __name__ == '__main__':
    main()

