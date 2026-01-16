"""
DGCNN Inference Package for LLP Classification

This package provides tools for running inference with a pre-trained
DGCNN (ParticleNet) model for Long-Lived Particle classification.

Usage:
    from inference import ParticleNet, H5GraphDataset
    
    # Or run from command line:
    python -m inference.run_inference --data /path/to/data.h5
"""

from .model import ParticleNet, ParticleDynamicEdgeConv, ParticleStaticEdgeConv
from .data_loader import H5GraphDataset, SimpleGraphDataset, worker_init_fn

__all__ = [
    'ParticleNet',
    'ParticleDynamicEdgeConv', 
    'ParticleStaticEdgeConv',
    'H5GraphDataset',
    'SimpleGraphDataset',
    'worker_init_fn',
]

__version__ = '1.0.0'

