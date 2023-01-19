import tensorflow as tf


y = tf.constant([1, 2, 3, 0, 2])
y = tf.one_hot(y, depth=4)
y = tf.cast(y, dtype=tf.float32)
out = tf.random.normal([5, 4])
loss1 = tf.reduce_mean(tf.square(y-out))
loss2 = tf.square(tf.norm(y-out))/(5*4)
print(tf.losses.MSE(y, out))
loss3 = tf.reduce_mean(tf.losses.MSE(y, out))
print(loss1, loss2, loss3)