# DGCNN Inference Package for LLP Classification

This package provides everything needed to run inference using our pre-trained DGCNN (Dynamic Graph CNN / ParticleNet) model for Long-Lived Particle (LLP) classification.

## 📁 Package Contents

```
inference/
├── README.md              # This file
├── config.yaml            # Model architecture and inference settings
├── model.py               # ParticleNet model definition
├── data_loader.py         # Data loading utilities
├── run_inference.py       # Main inference script
├── convert_root_to_h5.py  # ROOT to HDF5 converter
├── utils.py               # Utility functions
├── requirements.txt       # Python dependencies
└── weights/
    └── dgcnn_weights.pth  # Pre-trained model weights
```

## 🚀 Quick Start

### 1. Set Up Environment

```bash
# Create a new conda environment (recommended)
conda create -n dgcnn_inference python=3.10
conda activate dgcnn_inference

# Install all dependencies from requirements.txt
pip install -r requirements.txt
```

Or install manually:
```bash
# Install PyTorch with CUDA support
pip install torch --extra-index-url https://download.pytorch.org/whl/cu118

# Install PyTorch Geometric and extensions
pip install torch-geometric
pip install pyg_lib torch_scatter torch_sparse torch_cluster torch_spline_conv \
    -f https://data.pyg.org/whl/torch-2.5.0+cu118.html

# Install other dependencies
pip install h5py numpy PyYAML tqdm uproot awkward
```

### 2. Convert ROOT Files to HDF5 (if needed)

If your data is in ROOT format, first convert it to HDF5:

```bash
# Convert a single ROOT file (background data)
python convert_root_to_h5.py --input data.root --output data.h5

# Convert signal (LLP) data
python convert_root_to_h5.py --input signal.root --output signal.h5 --label signal

# Convert a directory of ROOT files
python convert_root_to_h5.py --input /path/to/root_files/ --output /path/to/h5_output/

# Use multiple workers for faster processing
python convert_root_to_h5.py --input /path/to/root_files/ --output /path/to/h5/ --workers 8
```

### 3. Run Inference

```bash
# Basic usage
python run_inference.py --data /path/to/your/data.h5 --output results.csv

# Process a directory of HDF5 files
python run_inference.py --data /path/to/data_directory/ --output results.csv

# CPU-only inference
python run_inference.py --data /path/to/data.h5 --device cpu --output results.csv

# Get only class probabilities
python run_inference.py --data /path/to/data.h5 --output-mode probabilities
```

## 📊 Input Data Format

### HDF5 Format (Recommended)

Your data should be in HDF5 format with the following structure:

```
your_data.h5
├── metadata/
│   └── attrs:
│       ├── num_graphs: int  # Number of graphs in file
│       └── rechit feature names: list  # Feature names (optional)
└── graph_{file_id}_{idx}/  # For each graph
    ├── pos: float32[N, 3]   # Node positions (eta, phi, energy)
    ├── x: float32[N, F]     # Node features (F=14 features)
    └── y: int64             # Label (optional, for evaluation)
```

### Expected Features (14 total)

The model expects 14 features per rechit in the following order:

| Index | Feature Name | Description |
|-------|-------------|-------------|
| 0 | rechit_eta | Pseudorapidity |
| 1 | rechit_phi | Azimuthal angle |
| 2 | rechit_energy | Energy deposit |
| 3 | rechit_time | Timing information |
| 4 | rechit_layer | Layer number |
| 5 | rechit_x | x-coordinate |
| 6 | rechit_y | y-coordinate |
| 7 | rechit_z | z-coordinate |
| 8 | rechit_rho | Radial distance |
| 9 | rechit_pt | Transverse momentum |
| 10 | rechit_fraction | Energy fraction |
| 11 | rechit_chi2 | Chi-squared |
| 12 | rechit_flags | Quality flags |
| 13 | rechit_size | Cluster size |

**Note:** Features 0-2 (eta, phi, energy) are used as spatial coordinates for k-NN graph construction.

### Using Custom Data

If your data is in a different format, you can use the `SimpleGraphDataset` class:

```python
from data_loader import SimpleGraphDataset
from torch_geometric.loader import DataLoader

# Your data as numpy arrays
features = [...]  # List of [N, 14] arrays
positions = [...]  # List of [N, 3] arrays (eta, phi, energy)

dataset = SimpleGraphDataset(features, positions)
loader = DataLoader(dataset, batch_size=64)
```

## ⚙️ Configuration

The `config.yaml` file contains all model and inference settings:

```yaml
model:
  input_features: 14      # Must match your data
  output_classes: 2       # [background, signal]
  conv_params:            # DO NOT MODIFY (must match trained model)
    - [16, [64, 64, 64]]
    - [16, [128, 128, 128]]
    - [16, [256, 256, 256]]
  fc_params:
    - [0.1, 256]

inference:
  batch_size: 64          # Adjust based on GPU memory
  device: auto            # 'cuda', 'cpu', or 'auto'
  output_mode: both       # 'labels', 'probabilities', or 'both'
```

## 📈 Output Format

### CSV Output (default)
```csv
predicted_label,prob_class_0,prob_class_1,true_label
0,0.923,0.077,0
1,0.124,0.876,1
...
```

### NPZ Output
```python
import numpy as np
results = np.load('results.npz')
labels = results['labels']           # Predicted labels
probs = results['probabilities']     # Class probabilities [N, 2]
true_labels = results['true_labels'] # Ground truth (if available)
```

## 🔧 Command Line Options

```
usage: run_inference.py [-h] --data DATA [--output OUTPUT] [--config CONFIG]
                        [--weights WEIGHTS] [--device DEVICE]
                        [--batch-size BATCH_SIZE]
                        [--output-mode {labels,probabilities,both}]
                        [--format {csv,npz,npy}] [--no-labels] [--quiet]

Options:
  --data, -d        Path to input data (file or directory)
  --output, -o      Output file path (default: inference_results.csv)
  --config, -c      Custom config file (default: config.yaml)
  --weights, -w     Custom weights file (overrides config)
  --device          Device: 'cuda', 'cuda:0', 'cpu' (overrides config)
  --batch-size, -b  Batch size (overrides config)
  --output-mode     Output: 'labels', 'probabilities', 'both'
  --format, -f      Output format: 'csv', 'npz', 'npy'
  --no-labels       Data does not contain ground truth labels
  --quiet, -q       Suppress progress bar
```

## 🐍 Python API

For programmatic use:

```python
import torch
from model import ParticleNet
from data_loader import H5GraphDataset
from torch_geometric.loader import DataLoader

# Load config
import yaml
with open('config.yaml') as f:
    config = yaml.safe_load(f)

# Setup device
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# Load model
settings = {
    'input_features': config['model']['input_features'],
    'output_classes': config['model']['output_classes'],
    'conv_params': config['model']['conv_params'],
    'fc_params': config['model']['fc_params'],
}
model = ParticleNet(settings)

# Load weights
checkpoint = torch.load('weights/dgcnn_weights.pth', map_location=device)
state_dict = checkpoint['model_state_dict']
# Remove 'module.' prefix if present
state_dict = {k.replace('module.', ''): v for k, v in state_dict.items()}
model.load_state_dict(state_dict)
model.to(device)
model.eval()

# Load data
dataset = H5GraphDataset(['/path/to/data.h5'])
loader = DataLoader(dataset, batch_size=64)

# Run inference
with torch.no_grad():
    for batch in loader:
        batch = batch.to(device)
        
        # Get predictions
        labels = model.predict(batch)
        
        # Get probabilities
        probs = model.predict_proba(batch)
        
        # Process results...
```

## 📋 Model Architecture

The model is a Dynamic Graph CNN (DGCNN) / ParticleNet architecture:

```
Input: Point cloud with 14 features per point
    ↓
BatchNorm(14)
    ↓
DynamicEdgeConv(k=16, 14 → 64 → 64 → 64) + Skip Connection
    ↓
DynamicEdgeConv(k=16, 64 → 128 → 128 → 128) + Skip Connection
    ↓
DynamicEdgeConv(k=16, 128 → 256 → 256 → 256) + Skip Connection
    ↓
Global Mean Pooling
    ↓
Linear(256 → 256) + Dropout(0.1) + ReLU
    ↓
Linear(256 → 2)
    ↓
Output: Logits for [background, signal]
```

## 🔄 ROOT to HDF5 Conversion

The `convert_root_to_h5.py` script converts CMS ROOT ntuple files to the HDF5 format expected by the inference pipeline.

### ROOT File Requirements

Your ROOT files should contain the `ntuples/llp` tree with CSC rechit data:
- `cscRechitsX`, `cscRechitsY`, `cscRechitsZ`, etc. (rechit features)
- `cscRechitsClusterId` (cluster assignments)
- `cscRechitClusterSize`, `cscRechitClusterTime`, etc. (cluster features)

For signal (LLP) samples, the following branches are used for truth matching:
- `gLLP_csc` (generator LLP in CSC)
- `cscRechitCluster_match_gParticle_minDeltaR`
- `cscRechitCluster_match_gParticle_id`

### Conversion Options

```
usage: convert_root_to_h5.py [-h] --input INPUT --output OUTPUT
                              [--label {signal,background,ggH,QCD,Gun}]
                              [--size-cutoff SIZE_CUTOFF] [--workers WORKERS]
                              [--no-labels]

Options:
  --input, -i       Input ROOT file or directory
  --output, -o      Output HDF5 file or directory
  --label, -l       Data label: 'signal' for LLP, 'background' for QCD/etc
  --size-cutoff     Minimum cluster size (default: 50 rechits)
  --workers, -w     Parallel workers for directory processing
  --no-labels       Omit ground truth labels from output
```

### Label Mapping

| Label | Numeric | Description |
|-------|---------|-------------|
| signal | 0 | LLP signal (ggH samples) |
| background | 1 | Background (QCD samples) |
| ggH | 0 | Alias for signal |
| QCD | 1 | Alias for background |

## ❓ Troubleshooting

### CUDA Out of Memory
- Reduce `batch_size` in config.yaml or via `--batch-size`
- Use CPU inference: `--device cpu`

### Missing PyG Extensions
```bash
# Check your PyTorch and CUDA versions first
python -c "import torch; print(torch.__version__, torch.cuda.is_available())"

# Install matching PyG extensions
pip install pyg_lib torch_scatter torch_sparse torch_cluster torch_spline_conv \
    -f https://data.pyg.org/whl/torch-{TORCH_VERSION}+{CUDA_VERSION}.html
```

### Data Format Issues
- Ensure your HDF5 files have the correct structure
- Check feature count matches `input_features` in config (should be 14)
- Verify position data is shape [N, 3]

## 📞 Contact

For questions or issues, please contact the MDS-ML team.

## 📄 License

This code is provided for research collaboration purposes.

