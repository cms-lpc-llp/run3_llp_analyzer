"""
DGCNN (ParticleNet) Model for LLP Classification
Adapted from the MDS-ML project for standalone inference.
"""

import torch
import torch_geometric


class ParticleStaticEdgeConv(torch_geometric.nn.MessagePassing):
    """Static edge convolution layer for point cloud processing."""
    
    def __init__(self, in_channels, out_channels):
        super(ParticleStaticEdgeConv, self).__init__(aggr='mean')
        self.mlp = torch.nn.Sequential(
            torch.nn.Linear(2 * in_channels, out_channels[0], bias=False),
            torch_geometric.nn.BatchNorm(out_channels[0]), 
            torch.nn.ReLU(),
            torch.nn.Linear(out_channels[0], out_channels[1], bias=False),
            torch_geometric.nn.BatchNorm(out_channels[1]),
            torch.nn.ReLU(),
            torch.nn.Linear(out_channels[1], out_channels[2], bias=False),
            torch_geometric.nn.BatchNorm(out_channels[2]),
            torch.nn.ReLU()
        )

    def forward(self, x, edge_index, k):
        return self.propagate(edge_index, size=(x.size(0), x.size(0)), x=x)

    def message(self, edge_index, x_i, x_j):
        tmp = torch.cat([x_i, x_j - x_i], dim=1)
        out_mlp = self.mlp(tmp)
        return out_mlp

    def update(self, aggr_out):
        return aggr_out


class ParticleDynamicEdgeConv(ParticleStaticEdgeConv):
    """Dynamic edge convolution layer with k-NN graph construction."""
    
    def __init__(self, in_channels, out_channels, k=7):
        super(ParticleDynamicEdgeConv, self).__init__(in_channels, out_channels)
        self.k = k
        self.skip_mlp = torch.nn.Sequential(
            torch.nn.Linear(in_channels, out_channels[2], bias=False),
            torch_geometric.nn.BatchNorm(out_channels[2]),
        )
        self.act = torch.nn.ReLU()

    def forward(self, pts, fts, batch=None):
        edges = torch_geometric.nn.knn_graph(pts, self.k, batch, loop=False, flow=self.flow)
        aggrg = super(ParticleDynamicEdgeConv, self).forward(fts, edges, self.k)
        x = self.skip_mlp(fts)
        out = torch.add(aggrg, x)
        return self.act(out)


class ParticleNet(torch.nn.Module):
    """
    ParticleNet/DGCNN model for point cloud classification.
    
    Architecture:
        - Input batch normalization
        - Dynamic graph convolution layers with skip connections
        - Global mean pooling
        - Fully connected layers with dropout
        - Output classification layer
    
    Args:
        settings (dict): Model configuration containing:
            - input_features (int): Number of input features per point
            - conv_params (list): List of [k, [c1, c2, c3]] for each conv layer
            - fc_params (list): List of [dropout_rate, units] for each FC layer
            - output_classes (int): Number of output classes
    """

    def __init__(self, settings):
        super().__init__()
        previous_output_shape = settings['input_features']

        self.input_bn = torch_geometric.nn.BatchNorm(settings['input_features'])
        
        self.conv_process = torch.nn.ModuleList()
        for layer_idx, layer_param in enumerate(settings['conv_params']):
            K, channels = layer_param
            self.conv_process.append(ParticleDynamicEdgeConv(previous_output_shape, channels, k=K))            
            previous_output_shape = channels[-1]

        self.fc_process = torch.nn.ModuleList()
        for layer_idx, layer_param in enumerate(settings['fc_params']):
            drop_rate, units = layer_param
            seq = torch.nn.Sequential(
                torch.nn.Linear(previous_output_shape, units),
                torch.nn.Dropout(p=drop_rate),
                torch.nn.ReLU()
            )
            self.fc_process.append(seq)
            previous_output_shape = units

        self.output_mlp_linear = torch.nn.Linear(previous_output_shape, settings['output_classes'])

    def forward(self, batch):
        """
        Forward pass.
        
        Args:
            batch: PyTorch Geometric Data batch with:
                - x: Node features [num_nodes, input_features]
                - pos: Node positions [num_nodes, 3] (used for k-NN graph)
                - batch: Batch assignment vector [num_nodes]
        
        Returns:
            torch.Tensor: Logits [batch_size, output_classes]
        """
        fts = self.input_bn(batch.x)
        pts = batch.pos

        for idx, layer in enumerate(self.conv_process):
            fts = layer(pts, fts, batch.batch)
            pts = fts

        x = torch_geometric.nn.global_mean_pool(fts, batch.batch)

        for layer in self.fc_process:
            x = layer(x)

        x = self.output_mlp_linear(x)
        return x
    
    def predict_proba(self, batch):
        """
        Get class probabilities (softmax of logits).
        
        Args:
            batch: PyTorch Geometric Data batch
            
        Returns:
            torch.Tensor: Probabilities [batch_size, output_classes]
        """
        logits = self.forward(batch)
        return torch.nn.functional.softmax(logits, dim=1)
    
    def predict(self, batch):
        """
        Get predicted class labels.
        
        Args:
            batch: PyTorch Geometric Data batch
            
        Returns:
            torch.Tensor: Predicted labels [batch_size]
        """
        logits = self.forward(batch)
        return torch.argmax(logits, dim=1)

