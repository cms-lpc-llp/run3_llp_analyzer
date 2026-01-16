"""
Utility functions for DGCNN inference.
"""

import numpy as np
import h5py
from pathlib import Path
from typing import List, Tuple, Optional, Union
import torch


def create_h5_from_arrays(
    output_path: str,
    features_list: List[np.ndarray],
    positions_list: List[np.ndarray],
    labels_list: Optional[List[int]] = None,
    feature_names: Optional[List[str]] = None,
    file_id: str = "0"
) -> None:
    """
    Create an HDF5 file from numpy arrays in the expected format.
    
    Args:
        output_path: Path to output HDF5 file
        features_list: List of feature arrays, each [num_nodes, num_features]
        positions_list: List of position arrays, each [num_nodes, 3]
        labels_list: Optional list of integer labels
        feature_names: Optional list of feature names
        file_id: Identifier for the file (used in graph naming)
    
    Example:
        >>> features = [np.random.randn(100, 14) for _ in range(10)]
        >>> positions = [np.random.randn(100, 3) for _ in range(10)]
        >>> create_h5_from_arrays('output.h5', features, positions)
    """
    if feature_names is None:
        feature_names = [
            'rechit_eta', 'rechit_phi', 'rechit_energy', 'rechit_time',
            'rechit_layer', 'rechit_x', 'rechit_y', 'rechit_z',
            'rechit_rho', 'rechit_pt', 'rechit_fraction', 'rechit_chi2',
            'rechit_flags', 'rechit_size'
        ]
    
    num_graphs = len(features_list)
    
    with h5py.File(output_path, 'w') as f:
        # Create metadata group
        metadata = f.create_group('metadata')
        metadata.attrs['num_graphs'] = num_graphs
        metadata.attrs['rechit feature names'] = str(feature_names)
        
        # Create graph groups
        for idx, (features, positions) in enumerate(zip(features_list, positions_list)):
            group_name = f'graph_{file_id}_{idx}'
            graph_group = f.create_group(group_name)
            
            graph_group.create_dataset('x', data=features.astype(np.float32))
            graph_group.create_dataset('pos', data=positions.astype(np.float32))
            
            if labels_list is not None:
                graph_group.create_dataset('y', data=labels_list[idx])
    
    print(f"Created HDF5 file: {output_path}")
    print(f"  - {num_graphs} graphs")
    print(f"  - {features_list[0].shape[1]} features per node")


def load_checkpoint_info(checkpoint_path: str) -> dict:
    """
    Load and display information about a checkpoint file.
    
    Args:
        checkpoint_path: Path to .pth checkpoint file
        
    Returns:
        Dictionary with checkpoint information
    """
    checkpoint = torch.load(checkpoint_path, map_location='cpu')
    
    info = {
        'keys': list(checkpoint.keys()),
    }
    
    if 'epoch' in checkpoint:
        info['epoch'] = checkpoint['epoch']
    if 'batch_step' in checkpoint:
        info['batch_step'] = checkpoint['batch_step']
    if 'config' in checkpoint:
        info['config'] = checkpoint['config']
    if 'model_state_dict' in checkpoint:
        state_dict = checkpoint['model_state_dict']
        info['num_parameters'] = sum(p.numel() for p in state_dict.values())
        info['model_layers'] = list(state_dict.keys())
    
    return info


def verify_data_format(h5_path: str) -> Tuple[bool, List[str]]:
    """
    Verify that an HDF5 file has the expected format.
    
    Args:
        h5_path: Path to HDF5 file
        
    Returns:
        Tuple of (is_valid, list_of_issues)
    """
    issues = []
    
    try:
        with h5py.File(h5_path, 'r') as f:
            # Check metadata
            if 'metadata' not in f:
                issues.append("Missing 'metadata' group")
            else:
                if 'num_graphs' not in f['metadata'].attrs:
                    issues.append("Missing 'num_graphs' attribute in metadata")
                else:
                    num_graphs = f['metadata'].attrs['num_graphs']
                    
                    # Check each graph
                    for idx in range(num_graphs):
                        # Try different naming patterns
                        found = False
                        for key in f.keys():
                            if key.startswith('graph_') and key != 'metadata':
                                found = True
                                graph_group = f[key]
                                
                                if 'x' not in graph_group:
                                    issues.append(f"Missing 'x' dataset in {key}")
                                if 'pos' not in graph_group:
                                    issues.append(f"Missing 'pos' dataset in {key}")
                                
                                # Check shapes
                                if 'x' in graph_group and 'pos' in graph_group:
                                    x_shape = graph_group['x'].shape
                                    pos_shape = graph_group['pos'].shape
                                    
                                    if len(x_shape) != 2:
                                        issues.append(f"{key}/x should be 2D, got {len(x_shape)}D")
                                    if len(pos_shape) != 2:
                                        issues.append(f"{key}/pos should be 2D, got {len(pos_shape)}D")
                                    if pos_shape[1] != 3:
                                        issues.append(f"{key}/pos should have 3 columns, got {pos_shape[1]}")
                                break
                        
                        if not found and idx == 0:
                            issues.append("No graph groups found")
                            break
    
    except Exception as e:
        issues.append(f"Error reading file: {e}")
    
    return len(issues) == 0, issues


def print_h5_structure(h5_path: str, max_graphs: int = 3) -> None:
    """
    Print the structure of an HDF5 file.
    
    Args:
        h5_path: Path to HDF5 file
        max_graphs: Maximum number of graphs to display details for
    """
    print(f"\nHDF5 Structure: {h5_path}")
    print("=" * 60)
    
    with h5py.File(h5_path, 'r') as f:
        def print_attrs(name, obj):
            indent = "  " * name.count('/')
            if isinstance(obj, h5py.Group):
                print(f"{indent}📁 {name}/")
                for key, val in obj.attrs.items():
                    val_str = str(val)[:50] + "..." if len(str(val)) > 50 else str(val)
                    print(f"{indent}    @{key}: {val_str}")
            elif isinstance(obj, h5py.Dataset):
                print(f"{indent}📊 {name}: {obj.dtype} {obj.shape}")
        
        # Print metadata
        if 'metadata' in f:
            print_attrs('metadata', f['metadata'])
        
        # Print graphs (limited)
        graph_count = 0
        for key in sorted(f.keys()):
            if key.startswith('graph_'):
                if graph_count < max_graphs:
                    print_attrs(key, f[key])
                    for subkey in f[key].keys():
                        print_attrs(f"{key}/{subkey}", f[key][subkey])
                graph_count += 1
        
        if graph_count > max_graphs:
            print(f"\n  ... and {graph_count - max_graphs} more graphs")
    
    print("=" * 60)


# Example usage and testing
if __name__ == '__main__':
    import sys
    
    if len(sys.argv) > 1:
        h5_path = sys.argv[1]
        print_h5_structure(h5_path)
        
        is_valid, issues = verify_data_format(h5_path)
        if is_valid:
            print("\n✅ Data format is valid!")
        else:
            print("\n❌ Data format issues found:")
            for issue in issues:
                print(f"  - {issue}")
    else:
        print("Usage: python utils.py <path_to_h5_file>")
        print("\nOr import and use programmatically:")
        print("  from utils import create_h5_from_arrays, verify_data_format")

