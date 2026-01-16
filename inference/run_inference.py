#!/usr/bin/env python3
"""
DGCNN Inference Script for LLP Classification
==============================================

This script runs inference using the pre-trained DGCNN model.

Usage:
    python run_inference.py --data /path/to/data --output results.csv
    python run_inference.py --data /path/to/data --config custom_config.yaml
    
For help:
    python run_inference.py --help
"""

import argparse
import os
import sys
from pathlib import Path
from typing import Optional, Union, List, Dict, Any
import yaml
import torch
import numpy as np
from tqdm import tqdm

# Import local modules
from model import ParticleNet
from data_loader import H5GraphDataset, SimpleGraphDataset, worker_init_fn

try:
    from torch_geometric.loader import DataLoader
except ImportError:
    from torch_geometric.data import DataLoader


def load_config(config_path: str) -> dict:
    """Load configuration from YAML file."""
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


def get_device(device_config: str) -> torch.device:
    """Determine the device to use for inference."""
    if device_config == 'auto':
        if torch.cuda.is_available():
            return torch.device('cuda')
        else:
            return torch.device('cpu')
    return torch.device(device_config)


def load_model(config: dict, weights_path: str, device: torch.device) -> ParticleNet:
    """
    Load the ParticleNet model with pre-trained weights.
    
    Args:
        config: Model configuration dict
        weights_path: Path to weights file (.pth)
        device: Device to load model onto
        
    Returns:
        Loaded model in eval mode
    """
    # Build model settings
    model_config = config['model']
    settings = {
        'input_features': model_config['input_features'],
        'output_classes': model_config['output_classes'],
        'conv_params': model_config['conv_params'],
        'fc_params': model_config['fc_params'],
    }
    
    # Initialize model
    model = ParticleNet(settings)
    
    # Load weights
    checkpoint = torch.load(weights_path, map_location=device)
    
    # Handle different checkpoint formats
    if 'model_state_dict' in checkpoint:
        state_dict = checkpoint['model_state_dict']
    elif 'state_dict' in checkpoint:
        state_dict = checkpoint['state_dict']
    else:
        state_dict = checkpoint
    
    # Remove 'module.' prefix if present (from DDP training)
    new_state_dict = {}
    for k, v in state_dict.items():
        if k.startswith('module.'):
            new_state_dict[k[7:]] = v
        else:
            new_state_dict[k] = v
    
    model.load_state_dict(new_state_dict)
    model.to(device)
    model.eval()
    
    return model


def run_inference(
    model: ParticleNet,
    data_loader: DataLoader,
    device: torch.device,
    output_mode: str = 'both',
    show_progress: bool = True
) -> Dict[str, np.ndarray]:
    """
    Run inference on a dataset.
    
    Args:
        model: Loaded ParticleNet model
        data_loader: DataLoader with inference data
        device: Device for inference
        output_mode: 'labels', 'probabilities', or 'both'
        show_progress: Whether to show progress bar
        
    Returns:
        Dictionary containing:
            - 'labels': Predicted labels (if output_mode includes labels)
            - 'probabilities': Class probabilities (if output_mode includes probabilities)
            - 'true_labels': Ground truth labels (if available in data)
    """
    all_labels = []
    all_probs = []
    all_true_labels = []
    
    iterator = tqdm(data_loader, desc="Running inference") if show_progress else data_loader
    
    with torch.no_grad():
        for batch in iterator:
            batch = batch.to(device)
            
            # Forward pass
            logits = model(batch)
            probs = torch.nn.functional.softmax(logits, dim=1)
            preds = torch.argmax(logits, dim=1)
            
            all_labels.append(preds.cpu().numpy())
            all_probs.append(probs.cpu().numpy())
            
            # Get true labels if available
            if hasattr(batch, 'y') and batch.y is not None:
                all_true_labels.append(batch.y.cpu().numpy())
    
    results = {}
    
    if output_mode in ['labels', 'both']:
        results['labels'] = np.concatenate(all_labels)
    
    if output_mode in ['probabilities', 'both']:
        results['probabilities'] = np.concatenate(all_probs)
    
    if all_true_labels:
        results['true_labels'] = np.concatenate(all_true_labels)
    
    return results


def save_results(
    results: Dict[str, np.ndarray],
    output_path: str,
    format: str = 'csv'
) -> None:
    """
    Save inference results to file.
    
    Args:
        results: Dictionary of results from run_inference
        output_path: Path to output file
        format: Output format ('csv', 'npz', or 'npy')
    """
    if format == 'csv':
        import csv
        with open(output_path, 'w', newline='') as f:
            writer = csv.writer(f)
            
            # Write header
            header = []
            if 'labels' in results:
                header.append('predicted_label')
            if 'probabilities' in results:
                num_classes = results['probabilities'].shape[1]
                header.extend([f'prob_class_{i}' for i in range(num_classes)])
            if 'true_labels' in results:
                header.append('true_label')
            writer.writerow(header)
            
            # Write data
            n_samples = len(results.get('labels', results.get('probabilities', results['true_labels'])))
            for i in range(n_samples):
                row = []
                if 'labels' in results:
                    row.append(int(results['labels'][i]))
                if 'probabilities' in results:
                    row.extend(results['probabilities'][i].tolist())
                if 'true_labels' in results:
                    row.append(int(results['true_labels'][i]))
                writer.writerow(row)
        print(f"Results saved to {output_path}")
    
    elif format == 'npz':
        np.savez(output_path, **results)
        print(f"Results saved to {output_path}")
    
    elif format == 'npy':
        # For npy, save probabilities as main array
        if 'probabilities' in results:
            np.save(output_path, results['probabilities'])
        else:
            np.save(output_path, results['labels'])
        print(f"Results saved to {output_path}")
    
    else:
        raise ValueError(f"Unknown format: {format}")


def compute_metrics(results: Dict[str, np.ndarray]) -> Dict[str, float]:
    """
    Compute evaluation metrics if true labels are available.
    
    Args:
        results: Dictionary containing 'labels' and 'true_labels'
        
    Returns:
        Dictionary of metrics
    """
    if 'true_labels' not in results or 'labels' not in results:
        return {}
    
    preds = results['labels']
    true = results['true_labels']
    
    metrics = {}
    metrics['accuracy'] = np.mean(preds == true)
    
    # Per-class accuracy
    for c in np.unique(true):
        mask = true == c
        metrics[f'accuracy_class_{c}'] = np.mean(preds[mask] == true[mask])
    
    # Confusion matrix elements
    for true_c in np.unique(true):
        for pred_c in np.unique(preds):
            mask = (true == true_c) & (preds == pred_c)
            metrics[f'count_true{true_c}_pred{pred_c}'] = int(np.sum(mask))
    
    return metrics


def main():
    parser = argparse.ArgumentParser(
        description='Run DGCNN inference for LLP classification',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic usage with HDF5 files
    python run_inference.py --data /path/to/data.h5 --output results.csv
    
    # Process multiple files
    python run_inference.py --data /path/to/data/*.h5 --output results.csv
    
    # Use custom config
    python run_inference.py --data /path/to/data --config custom_config.yaml
    
    # CPU-only inference
    python run_inference.py --data /path/to/data --device cpu
        """
    )
    
    parser.add_argument(
        '--data', '-d',
        type=str,
        required=True,
        help='Path to input data (HDF5 file or directory containing HDF5 files)'
    )
    parser.add_argument(
        '--output', '-o',
        type=str,
        default='inference_results.csv',
        help='Path to output file (default: inference_results.csv)'
    )
    parser.add_argument(
        '--config', '-c',
        type=str,
        default=None,
        help='Path to config file (default: config.yaml in script directory)'
    )
    parser.add_argument(
        '--weights', '-w',
        type=str,
        default=None,
        help='Path to model weights (overrides config)'
    )
    parser.add_argument(
        '--device',
        type=str,
        default=None,
        help='Device for inference: cuda, cuda:0, cpu (overrides config)'
    )
    parser.add_argument(
        '--batch-size', '-b',
        type=int,
        default=None,
        help='Batch size (overrides config)'
    )
    parser.add_argument(
        '--output-mode',
        type=str,
        choices=['labels', 'probabilities', 'both'],
        default=None,
        help='Output mode (overrides config)'
    )
    parser.add_argument(
        '--format', '-f',
        type=str,
        choices=['csv', 'npz', 'npy'],
        default='csv',
        help='Output format (default: csv)'
    )
    parser.add_argument(
        '--no-labels',
        action='store_true',
        help='Data does not contain ground truth labels'
    )
    parser.add_argument(
        '--quiet', '-q',
        action='store_true',
        help='Suppress progress bar'
    )
    
    args = parser.parse_args()
    
    # Determine script directory for default paths
    script_dir = Path(__file__).parent.resolve()
    
    # Load config
    if args.config:
        config_path = Path(args.config)
    else:
        config_path = script_dir / 'config.yaml'
    
    if not config_path.exists():
        print(f"Error: Config file not found: {config_path}")
        sys.exit(1)
    
    config = load_config(str(config_path))
    print(f"Loaded config from: {config_path}")
    
    # Determine weights path
    if args.weights:
        weights_path = Path(args.weights)
    else:
        weights_rel_path = config['weights']['path']
        weights_path = script_dir / weights_rel_path
    
    if not weights_path.exists():
        print(f"Error: Weights file not found: {weights_path}")
        print("Please ensure the model weights are in the 'weights' directory.")
        sys.exit(1)
    
    print(f"Using weights: {weights_path}")
    
    # Get device
    device_str = args.device if args.device else config['inference'].get('device', 'auto')
    device = get_device(device_str)
    print(f"Using device: {device}")
    
    # Load model
    print("Loading model...")
    model = load_model(config, str(weights_path), device)
    print(f"Model loaded: {config['model']['name']}")
    
    # Prepare data
    data_path = Path(args.data)
    if data_path.is_file():
        file_paths = [data_path]
    elif data_path.is_dir():
        file_paths = list(data_path.glob('*.h5'))
        if not file_paths:
            print(f"Error: No HDF5 files found in {data_path}")
            sys.exit(1)
    else:
        # Try glob pattern
        import glob
        file_paths = [Path(p) for p in glob.glob(args.data)]
        if not file_paths:
            print(f"Error: No files found matching {args.data}")
            sys.exit(1)
    
    print(f"Found {len(file_paths)} data file(s)")
    
    # Create dataset
    return_labels = not args.no_labels
    dataset = H5GraphDataset(file_paths, return_labels=return_labels)
    print(f"Dataset loaded: {len(dataset)} samples")
    
    # Create data loader
    batch_size = args.batch_size if args.batch_size else config['inference'].get('batch_size', 64)
    num_workers = config['inference'].get('num_workers', 4)
    
    data_loader = DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=False,
        num_workers=num_workers,
        pin_memory=(device.type == 'cuda'),
        worker_init_fn=worker_init_fn if num_workers > 0 else None
    )
    
    # Run inference
    output_mode = args.output_mode if args.output_mode else config['inference'].get('output_mode', 'both')
    print(f"\nRunning inference (batch_size={batch_size}, output_mode={output_mode})...")
    
    results = run_inference(
        model, 
        data_loader, 
        device, 
        output_mode=output_mode,
        show_progress=not args.quiet
    )
    
    # Compute and display metrics if labels available
    if 'true_labels' in results:
        metrics = compute_metrics(results)
        print("\n" + "="*50)
        print("EVALUATION METRICS")
        print("="*50)
        for name, value in metrics.items():
            if 'count' not in name:
                print(f"  {name}: {value:.4f}")
        print("="*50)
    
    # Save results
    output_path = Path(args.output)
    save_results(results, str(output_path), format=args.format)
    
    print(f"\nInference complete! Results saved to: {output_path}")
    
    # Summary
    n_samples = len(results.get('labels', results.get('probabilities')))
    print(f"\nSummary:")
    print(f"  Total samples: {n_samples}")
    if 'labels' in results:
        unique, counts = np.unique(results['labels'], return_counts=True)
        print(f"  Predictions:")
        for label, count in zip(unique, counts):
            class_name = 'signal' if label == 1 else 'background'
            print(f"    Class {label} ({class_name}): {count} ({100*count/n_samples:.1f}%)")


if __name__ == '__main__':
    main()

