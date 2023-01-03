import numpy as np
import matplotlib.pyplot as plt
from matplotlib import use
from sklearn.datasets import make_blobs


X, y = make_blobs(n_samples=1000, n_features=2, centers=[
                  [-1, -1], [0, 0], [1, 1], [2, 2]], cluster_std=[0.4, 0.2, 0.2, 0.2], random_state=9)
plt.scatter(X[:, 0], X[:, 1], marker='o')
plt.show()
