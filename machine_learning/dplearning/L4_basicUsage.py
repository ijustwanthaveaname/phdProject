import tensorflow as tf
import numpy as np


tf.constant(1)  # similar to convert_to_tensor
tf.constant(1.)
tf.constant("strings")
tf.constant([True, False])

with tf.device("cpu"):
    a = tf.constant(1)
with tf.device("gpu"):
    b = tf.range(4)
a.device
b.device
aa = a.gpu()
bb = b.cpu()
aa.device
bb.device

b.numpy()
b.ndim
b.shape
# return dimension in tensor form
tf.rank(b)
tf.rank(tf.ones([3, 3, 3]))
# determine if type of tensor
isinstance(a ,tf.Tensor)  # bad in some situation
tf.is_tensor(a)
a.dtype == tf.float32
a = np.arange(5)
aa = tf.convert_to_tensor(a, dtype=tf.int32)
# convert dtype
tf.cast(aa, dtype=tf.float32)
# variable, parameters
b = tf.Variable(a, name="input_data")
b.name
b.trainable
isinstance(b, tf.Tensor)  # False
isinstance(b, tf.Variable) # True
tf.is_tensor(b)  # true
a = tf.ones([])  # 1
a
int(a)
float(a)
tf.reshape(aa,[5, -1])
tf.zeros_like(aa)  # ones_like
tf.zeros(aa.shape)  # ones() is same 
tf.fill([2, 4], 9)
tf.random.normal([2, 2], mean=1, stddev=1)
tf.random.truncated_normal([2, 2], mean=0, stddev=1)  # avoid gradient vanish
tf.random.uniform([2, 2], 0, 1)
idx = tf.range(10)
idx = tf.random.shuffle(idx)
a = tf.random.normal([10, 784])
b = tf.random.uniform([10], maxval=10, dtype=tf.int32)
a = tf.gather(a, idx)
b = tf.gather(b, idx)

# calculate loss
out = tf.random.uniform([4, 10])  # 4 images with potential 10 classes
y = tf.range(4)  # 4 classes
y = tf.one_hot(y, depth=10)  # one-hot
loss = tf.keras.losses.mse(y, out)
loss = tf.reduce_mean(loss)
# create a layer
x = tf.random.normal([4, 784])
net = tf.keras.layers.Dense(10)
net.build((4, 784))
net(x)
net.kernel  # w 
net.bias  # b

a = tf.ones([4, 35 ,8])
tf.gather(a, axis=0, indices=[1, 2])
tf.gather_nd(a, [0, 1, 2])
tf.gather_nd(a, [[0, 1, 2]])
tf.gather_nd(a, [0, 1])

a = tf.ones([2, 3, 4])
tf.boolean_mask(a, mask=[[True, False, False], [False, True, True]])
tf.transpose(a, perm=[0, 2, 1])

tf.expand_dims(a, axis=0).shape
tf.expand_dims(a,axis=1).shape # +1 location
tf.expand_dims(a, axis=-2).shape # +0 location
a = tf.zeros([1, 2, 1, 1, 3])
tf.squeeze(a).shape
tf.squeeze(a, axis=2).shape

b = tf.broadcast_to(tf.random.normal([4, 1, 1, 1]), [4, 32, 32, 3])
b.shape

b = tf.tile(tf.random.normal([1, 32, 1]), [4, 1, 32])
b.shape