import torch_geometric
from torch_geometric.data import Data
import torch


# define 4 nodes with 2 features
x = torch.tensor([[2, 1], [5, 6], [3, 7], [12, 0]], dtype=torch.float)
# define labels
y = torch.tensor([0, 1, 0, 1], dtype=torch.float)
# define edges
edge_index = torch.tensor([[0, 2, 1, 0, 3],
                           [3, 1, 0, 1, 2]], dtype=torch.long)
# create data
data = Data(x=x, y=y, edge_index=edge_index)