import os
import numpy as np
import matplotlib.pyplot as plt
from tensorflow import keras

os.chdir("/home/weishen/phdProject/machine_learning/dplearning/RNN/")
#### 模拟数据生成函数：多个正弦波+随机噪音
def generate_time_series(batch_size, n_steps,seed=10):
    np.random.seed(seed)
    freq1, freq2, offsets1, offsets2 = np.random.rand(4, batch_size, 1)
    time = np.linspace(0, 1, n_steps)
    series = 0.5 * np.sin((time - offsets1) * (freq1 * 10 + 10)) # wave 1
    series += 0.2 * np.sin((time - offsets2) * (freq2 * 20 + 20)) # + wave 2
    series += 0.1 * (np.random.rand(batch_size, n_steps) - 0.5) # + noise
    return series[..., np.newaxis].astype(np.float32)
#### 生成模拟时间序列以及训练集和测试集
n_steps = 50
series = generate_time_series(10000, n_steps + 1)
X_train, Y_train = series[:7000, :n_steps], series[:7000, -1]
X_valid, Y_valid = series[7000:9000, :n_steps], series[7000:9000, -1]
X_test, Y_test = series[9000:, :n_steps], series[9000:, -1]
#### linear regression
model = keras.models.Sequential([
    keras.layers.Flatten(input_shape=[50, 1]),
    keras.layers.Dense(1)
])
model.compile(loss=keras.losses.mean_squared_error,
              optimizer=keras.optimizers.Adam(0.01))
model.fit(X_train,Y_train,epochs=20,verbose=0)
model.evaluate(X_valid,Y_valid)
#### simple RNN
model = keras.models.Sequential([
    keras.layers.SimpleRNN(1, input_shape=(None, 1))
])
model.compile(loss=keras.losses.mean_squared_error,
              optimizer=keras.optimizers.Adam(0.01))
model.fit(X_train,Y_train,epochs=20,verbose=0)
model.evaluate(X_valid,Y_valid)
#### deep RNN:LSTM
model = keras.models.Sequential([
    keras.layers.LSTM(20, return_sequences=True, input_shape=[None, 1]),
    keras.layers.LSTM(20),
    keras.layers.Dense(1)
])
model.compile(loss=keras.losses.mean_squared_error,
              optimizer=keras.optimizers.Adam(0.01))
model.fit(X_train,Y_train,epochs=20,verbose=0)
model.evaluate(X_valid,Y_valid)
# deep RNN:GRU
model = keras.models.Sequential([
    keras.layers.GRU(20, return_sequences=True, input_shape=[None, 1]),
    keras.layers.GRU(20),
    keras.layers.Dense(1)
])
model.compile(loss=keras.losses.mean_squared_error,
              optimizer=keras.optimizers.Adam(0.01))
model.fit(X_train,Y_train,epochs=20,verbose=0)
model.evaluate(X_valid,Y_valid)
#### 10次预测10个值:Sequence-to-Vector RNN
series = generate_time_series(2000, n_steps + 10)
X_valid, y_valid = series[:, :n_steps], series[:,n_steps:]
X = X_valid
for step_ahead in range(10):
    y_pred_one = model.predict(X[:, step_ahead:])[:, np.newaxis, :]
    X = np.concatenate([X, y_pred_one], axis=1)
y_pred = X[:, n_steps:]
mse=keras.losses.mean_squared_error(y_valid,y_pred)
print(mse.numpy().mean())

# 1次预测10个值:Sequence-to-Vector RNN, 标签为10维
series = generate_time_series(10000, n_steps + 10)
X_train, Y_train = series[:7000, :n_steps], series[:7000, -10:, 0]
X_valid, Y_valid = series[7000:9000, :n_steps], series[7000:9000, -10:, 0]
X_test, Y_test = series[9000:, :n_steps], series[9000:, -10:, 0]
model = keras.models.Sequential([
    keras.layers.SimpleRNN(20, return_sequences=True, input_shape=[None, 1]),
    keras.layers.SimpleRNN(20),
    keras.layers.Dense(10)
])
model.compile(loss=keras.losses.mean_squared_error,
              optimizer=keras.optimizers.Adam(0.01))
model.fit(X_train,Y_train,epochs=20,verbose=0)
model.evaluate(X_valid,Y_valid)

# 1次预测10个值:1-d-CNN+Sequence-to-Sequence RNN，第一层加CNN可能会提高准确率
model = keras.models.Sequential([
    keras.layers.Conv1D(filters=20, kernel_size=4, strides=2, padding="valid",
    input_shape=[None, 1]),
    keras.layers.SimpleRNN(20, return_sequences=True),
    keras.layers.SimpleRNN(20, return_sequences=True),
    keras.layers.TimeDistributed(keras.layers.Dense(10))
])
model.compile(loss="mse", optimizer=keras.optimizers.Adam(0.01),metrics=[last_time_step_mse])
history = model.fit(X_train, Y_train[:, 3::2], epochs=20,verbose=0)
model.evaluate(X_valid,Y_valid[:, 3::2])