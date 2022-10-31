#!/usr/bin/env python
# coding: utf-8

import numpy as np
import torch
from torch_geometric.data import Data
import os.path as osp
import os
import torch
import torch.nn.functional as F

from torch import nn
from torch.nn import Linear
from torch.nn import BatchNorm1d
from torch.utils.data import Dataset
from torch_geometric.nn import GCNConv, GATConv
from torch_geometric.data import DataLoader
from torch_geometric.nn import global_add_pool, global_mean_pool

n_features = 396

class Net(torch.nn.Module):
    def __init__(self):
        super(Net, self).__init__()
        self.conv1 = GCNConv(n_features, 4096, cached=False)
        self.bn1 = BatchNorm1d(4096)
        self.conv2 = GCNConv(4096, 2048, cached=False)
        self.bn2 = BatchNorm1d(2048)
        self.conv3 = GCNConv(2048, 1024, cached=False)
        self.bn3 = BatchNorm1d(1024)
        self.conv4 = GCNConv(1024, 1024, cached=False)
        self.bn4 = BatchNorm1d(1024)
        self.conv5 = GCNConv(1024, 2048, cached=False)
        self.bn5 = BatchNorm1d(2048)
        self.conv6 = GCNConv(2048, 256, cached=False)
        self.bn6 = BatchNorm1d(256)
        self.conv7 = GCNConv(256, 1, cached=False)

        self.fc2 = Linear(2048, 1024)
        self.fc3 = Linear(1024, 512)
        self.fc4 = Linear(512, 1)

    def forward(self, data):
        x, edge_index, batch = data.x, data.edge_index, data.batch
        x = F.relu(self.conv1(x, edge_index))
        x = self.bn1(x)
        x = F.relu(self.conv2(x, edge_index))
        x = self.bn2(x)
        x = F.relu(self.conv3(x, edge_index))
        x = self.bn3(x)
        x = F.relu(self.conv4(x, edge_index))
        x = self.bn4(x)
        x = F.relu(self.conv5(x, edge_index))
        x = self.bn5(x)
        x = F.relu(self.conv6(x, edge_index))
        x = self.bn6(x)
        x = self.conv7(x, edge_index)
        x = global_add_pool(x, batch)
        return x

def load_model(device="cpu"):
    root_path = os.path.abspath(os.path.dirname( os.path.dirname(__file__)) )
    
    nmodel= Net().to(device)
    nmodel_file = os.path.join(root_path, "moltaut_weights/neutral.pth")
    nweights = torch.load(nmodel_file, map_location=device)
    nmodel.load_state_dict(nweights, strict=True)
    nmodel.eval()

    imodel= Net().to(device)
    imodel_file = os.path.join(root_path, "moltaut_weights/ionic.pth")
    iweights = torch.load(imodel_file, map_location=device)
    imodel.load_state_dict(iweights, strict=True)
    imodel.eval()
    return nmodel, imodel


