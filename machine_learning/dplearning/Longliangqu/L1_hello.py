import tensorflow as tf
import time
import numpy as np


# cpu/gpu computation
start = time.time()
cpu_a = tf.random.normal([10000, 1000])
cpu_b = tf.random.normal([1000, 2000])
cpu_c = tf.matmul(cpu_a, cpu_b)
end = time.time()
print(end - start)
start = time.time()
np_a = np.random.normal(0, 1, [10000,1000])
np_b = np.random.normal(0, 1, [1000, 2000])
np_c = np_a @ np_b
end = time.time()
print(end - start)
# Derivative
x = tf.constant(1.)
a = tf.constant(2.)
b = tf.constant(3.)
c = tf.constant(4.)

with tf.GradientTape() as tape:
    tape.watch([a, b, c])
    y = a**2 * x + b * x + c

[dy_da, dy_db, dy_dc] = tape.gradient(y, [a, b, c])
print(dy_da, dy_db, dy_dc)
# test GPU
tf.config.list_physical_devices("GPU")