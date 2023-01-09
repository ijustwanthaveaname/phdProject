import tensorflow as tf
import matplotlib.pyplot as plt


(train_image, train_label), (test_image, test_label) = tf.keras.datasets.fashion_mnist.load_data()
train_image.shape
train_label.shape


model = tf.keras.Sequential()
model.add(tf.keras.layers.Flatten(input_shape=(28, 28)))
model.add(tf.keras.layers.Dense(128, activation="relu"))
model.add(tf.keras.layers.Dense(10, activation="softmax"))
model.compile(
    optimizer="adam",  # tf.keras.optimizers.Adam(learning_rate=0.001)
    loss="sparse_categorical_crossentropy",
    metrics=["acc"]
)
model.fit(train_image, train_label, batch_size=int(train_image.shape[0]/100), epochs=100)
model.evaluate(test_image, test_label)


# acc and overfit 
history = model.fit(train_image, train_label, 
                    epochs=10, validation_data=(test_image, test_label))
plt.plot(history.epoch, history.history.get("loss"), label="loss")
plt.plot(history.epoch, history.history.get("val_loss"), label="val_loss")
plt.legend()

# dropout
model2 = tf.keras.Sequential()
model2.add(tf.keras.layers.Flatten(input_shape=(28, 28)))
model2.add(tf.keras.layers.Dense(128, activation="relu"))
model2.add(tf.keras.layers.Dropout(0.5))
model2.add(tf.keras.layers.Dense(10, activation="softmax"))
model2.compile(
    optimizer="adam",  # tf.keras.optimizers.Adam(learning_rate=0.001)
    loss="sparse_categorical_crossentropy",
    metrics=["acc"]
)
model2.fit(train_image, train_label, batch_size=int(train_image.shape[0]/100), epochs=100)
model2.evaluate(test_image, test_label)