import tensorflow as tf
from sklearn.datasets import load_boston
import matplotlib.pyplot as plt

boston = load_boston()

x = boston.data
y = boston.target
model = tf.keras.Sequential()
model.add(tf.keras.layers.Dense(1, input_shape=(13,)))
model.compile(
        optimizer="adam",
        loss="mse"
)
history = model.fit(x, y, epochs=100)
plt.plot(history.history.epoch, history.history.get("loss"))

y_pred = model.predict(x)
model.evaluate(x, y)