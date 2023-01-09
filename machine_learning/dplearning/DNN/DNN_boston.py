from sklearn.datasets import load_boston
import tensorflow as tf
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error


def load_data():
    boston = load_boston()
    boston.feature_names
    return boston.data, boston.target


def dnn(x, y):
    print(x.shape)
    model = tf.keras.Sequential()
    model.add(tf.keras.layers.Dense(15, activation="relu", input_shape=(13, )))
    model.add(tf.keras.layers.Dense(15, activation="relu"))
    model.add(tf.keras.layers.Dense(1))
    model.compile(
        optimizer="adam",
        loss="mse"
    )
    model.fit(x, y, epochs=20)
    y_pred = model.predict(x)
    print(mean_squared_error(y_pred, y))
    print(model.evaluate(x, y))


def main():
    x, y = load_data()
    dnn(x, y)


if __name__ == "__main__":
    main()
