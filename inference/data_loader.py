"""
Data loading utilities for DGCNN inference.
Supports HDF5 files in the MDS-ML format.
"""

import h5py
import torch
import os
import re
from torch.utils.data import Dataset
from torch_geometric.data import Data
from pathlib import Path
import atexit


# Global cache for HDF5 handles (per-worker process)
_worker_h5_handles = {}
_worker_id_for_cleanup_debug = -1


def _close_all_h5_handles():
    """Closes all HDF5 handles cached by a worker process."""
    for path_str, handle in list(_worker_h5_handles.items()):
        try:
            if handle and hasattr(handle, 'close') and callable(handle.close):
                handle.close()
        except Exception as e:
            print(f"Error closing {path_str}: {e}", flush=True)
    _worker_h5_handles.clear()


def worker_init_fn(worker_id):
    """PyTorch DataLoader worker_init_fn to register HDF5 handle cleanup."""
    global _worker_id_for_cleanup_debug
    _worker_id_for_cleanup_debug = worker_id
    atexit.register(_close_all_h5_handles)


class H5GraphDataset(Dataset):
    """
    Dataset for loading graph data from HDF5 files.
    
    Expected HDF5 structure:
        /metadata
            attrs: num_graphs, rechit feature names, etc.
        /graph_{file_id}_{graph_idx}
            /pos: Position features [num_nodes, 3]
            /x: Node features [num_nodes, num_features]
            /y: Label (scalar)
    
    Args:
        file_paths: List of paths to HDF5 files
        return_labels: Whether to include labels (set False if labels unavailable)
    """
    
    def __init__(self, file_paths, return_labels=True):
        self.file_paths = [str(fp) for fp in file_paths]
        self.return_labels = return_labels
        self.graph_map = []
        self.metadata = {}
        
        # Build graph map from all files
        for file_path_str in self.file_paths:
            try:
                with h5py.File(file_path_str, 'r') as h5_file:
                    num_graphs = h5_file['metadata'].attrs['num_graphs']
                    for graph_idx in range(num_graphs):
                        self.graph_map.append((file_path_str, graph_idx))
            except Exception as e:
                print(f"Warning: Could not read {file_path_str}: {e}", flush=True)
        
        if not self.graph_map:
            raise ValueError(f"No graphs found in provided files: {self.file_paths}")
        
        # Load metadata from first file
        first_file = self.graph_map[0][0]
        try:
            with h5py.File(first_file, 'r') as h5_file:
                for key, val in h5_file['metadata'].attrs.items():
                    if isinstance(val, bytes):
                        try:
                            val = val.decode('utf-8')
                        except UnicodeDecodeError:
                            pass
                    try:
                        import ast
                        self.metadata[key] = ast.literal_eval(val) if isinstance(val, str) else val
                    except (ValueError, SyntaxError):
                        self.metadata[key] = val
                self.metadata["num_graphs_total"] = len(self.graph_map)
        except Exception as e:
            print(f"Warning: Could not read metadata from {first_file}: {e}", flush=True)

    def __len__(self):
        return len(self.graph_map)
    
    def _get_h5_handle(self, file_path_str):
        """Gets a cached HDF5 file handle or opens a new one."""
        if file_path_str not in _worker_h5_handles or not _worker_h5_handles[file_path_str].id.valid:
            if file_path_str in _worker_h5_handles:
                try:
                    _worker_h5_handles[file_path_str].close()
                except Exception:
                    pass
            _worker_h5_handles[file_path_str] = h5py.File(file_path_str, 'r')
        return _worker_h5_handles[file_path_str]

    def __getitem__(self, idx):
        if not (0 <= idx < len(self.graph_map)):
            raise IndexError(f"Index {idx} out of bounds for dataset length {len(self.graph_map)}")
        
        file_path_str, graph_idx_in_file = self.graph_map[idx]
        
        # Extract file identifier from filename (e.g., 'combined_123.h5' -> '123')
        filename = os.path.basename(file_path_str)
        match = re.search(r'combined_(\d+)\.h5', filename)
        if match:
            file_numeric_identifier = match.group(1)
        else:
            # Fallback: use file index
            file_numeric_identifier = str(self.file_paths.index(file_path_str))
        
        graph_group_name = f'graph_{file_numeric_identifier}_{graph_idx_in_file}'
        
        h5_file = self._get_h5_handle(file_path_str)
        
        if graph_group_name not in h5_file:
            raise KeyError(f"Graph '{graph_group_name}' not found in {file_path_str}")
        
        graph_group = h5_file[graph_group_name]
        
        pos_data = graph_group['pos'][:]
        x_data = graph_group['x'][:]
        
        pos = torch.from_numpy(pos_data).float()
        x = torch.from_numpy(x_data).float()
        
        if self.return_labels and 'y' in graph_group:
            y_data = graph_group['y'][()]
            y = torch.tensor(y_data, dtype=torch.long)
            # Remap label 2 to 1 (if applicable to your data)
            y = torch.where(y == 2, torch.tensor(1, dtype=torch.long), y)
            return Data(x=x, pos=pos, y=y)
        else:
            return Data(x=x, pos=pos)


class SimpleGraphDataset(Dataset):
    """
    Simple dataset for custom data not in HDF5 format.
    
    Use this for raw numpy arrays or tensors.
    
    Args:
        features: List of feature arrays, each [num_nodes, num_features]
        positions: List of position arrays, each [num_nodes, 3]
        labels: Optional list of labels
    """
    
    def __init__(self, features, positions, labels=None):
        assert len(features) == len(positions), "Features and positions must have same length"
        if labels is not None:
            assert len(features) == len(labels), "Features and labels must have same length"
        
        self.features = features
        self.positions = positions
        self.labels = labels
    
    def __len__(self):
        return len(self.features)
    
    def __getitem__(self, idx):
        x = torch.tensor(self.features[idx], dtype=torch.float32)
        pos = torch.tensor(self.positions[idx], dtype=torch.float32)
        
        if self.labels is not None:
            y = torch.tensor(self.labels[idx], dtype=torch.long)
            return Data(x=x, pos=pos, y=y)
        else:
            return Data(x=x, pos=pos)

