from tensorflow.keras.layer import Input, Lambda, Dense, concatenate
from tensorflow.keras.models import Model


# input (i1 ,i2) hidden (h1, h2) output (o1, o2) for just one hidden layer
# h1 = w2*i2, others are fully connected
inp = Input(shape=(2,))
inp2 = Lambda(lambda x: x[:,1:2])(inp)   # get the second neuron, tf.gather(x, [1,2,5], axis=1) 

h1_out = Dense(1, activation='sigmoid')(inp2)  # only connected to the second neuron
h2_out = Dense(1, activation='sigmoid')(inp)  # connected to both neurons
h_out = concatenate([h1_out, h2_out])

out = Dense(2, activation='sigmoid')(h_out)

model = Model(inp, out)

# simply train it using `fit`
model.fit(...)