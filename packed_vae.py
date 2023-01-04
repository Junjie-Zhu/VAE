import os
import random as rn
import sys

import numpy as np
import biobox as bb
import tensorflow as tf
from keras import backend as K
from keras.layers import Input, Dense, Lambda
from keras.models import Model
from scipy import stats
from sklearn.preprocessing import MinMaxScaler

os.environ["CUDA_VISIBLE_DEVICES"] = "1"

np.random.seed(123)
rn.seed(123)
tf.random.set_seed(1)

infile = sys.argv[1]
encoding_dim = 2
BATCH_SIZE = 64
EPOCHS = 50


def main():
    # load data
    train_infile = infile.replace(".pdb", "_train.dat")
    test_infile = infile.replace(".pdb", "_test.dat")

    x_train_orig = np.loadtxt(train_infile)
    x_test_orig = np.loadtxt(test_infile)

    # normalization of input data
    scaler = MinMaxScaler(feature_range=(0, 1))
    x_train = scaler.fit_transform(x_train_orig)
    x_test = scaler.fit_transform(x_test_orig)

    integral_model = vae_structure(x_train)
    encoder = integral_model[0]
    sampling_model = integral_model[1]

    history = sampling_model.fit(x_train, x_train,
                                 epochs=EPOCHS,
                                 batch_size=BATCH_SIZE,
                                 shuffle=True,
                                 validation_data=(x_test, x_test))

    # save reconstructed structures
    x_train = scaler.fit_transform(x_train_orig)
    rc_train = sampling_model.predict(x_train)
    reconstruct_train = scaler.inverse_transform(rc_train)
    decoded_reshaped_train = reconstruct_train.reshape(reconstruct_train.shape[0],
                                                       int(reconstruct_train.shape[1] / 3), 3)

    x_test = scaler.fit_transform(x_test_orig)
    rc_test = sampling_model.predict(x_test)
    reconstruct_test = scaler.inverse_transform(rc_test)
    decoded_reshaped_test = reconstruct_test.reshape(reconstruct_test.shape[0], int(reconstruct_test.shape[1] / 3), 3)

    spr = spearman_corr(x_train_orig, reconstruct_train)
    spr_test = spearman_corr(x_test_orig, reconstruct_test)

    print("> spearman corelation coefficients:\n  train set: %.3f\n  test set: %.3f" % (spr, spr_test))

    M = bb.Molecule()
    M.import_pdb(infile)
    idx = M.atomselect("*", "*", ["CA", "CB", "C", "N", "CG", "CD", "NE2", "CH3", "O", "OE1"], get_index=True)[1]

    M4 = M.get_subset(idxs=idx, conformations=[0])
    M4.coordinates = decoded_reshaped_test
    M4.write_pdb("./data/%s" % infile.replace(".pdb", ".generated.pdb"))


def spearman_corr(y_true, y_pred):
    sprm = 0
    for i in range(y_true.shape[0]):
        sprm += stats.spearmanr(y_true[i], y_pred[i])[0]
    sprm = sprm / y_true.shape[0]
    return sprm


def vae_structure(x_train):
    # encoder
    input_prot = Input(shape=(x_train.shape[1],))

    encoded = Dense(1024, activation='relu')(input_prot)
    encoded = Dense(256, activation='relu')(encoded)
    encoded = Dense(64, activation='relu')(encoded)
    encoded = Dense(16, activation='relu')(encoded)

    encode_mean = Dense(encoding_dim)(encoded)
    encode_log_var = Dense(encoding_dim)(encoded)

    sample = Lambda(sampling, name='sampling')([encode_mean, encode_log_var])

    encoder = Model(input_prot, sample)

    # decoder
    decoded = Dense(16, activation='relu')(sample)
    decoded = Dense(64, activation='relu')(decoded)
    decoded = Dense(256, activation='relu')(decoded)
    decoded = Dense(1024, activation='relu')(decoded)
    decoded = Dense(x_train.shape[1], activation='sigmoid')(decoded)

    # decoder = Model([encode_mean, encode_log_var], decoded)
    sampling_model = Model(input_prot, decoded)

    # define loss, optimizer
    msle = tf.keras.losses.MeanSquaredLogarithmicError()
    xent_loss = K.sum(msle(input_prot, decoded))
    kl_loss = -5e-4 * K.mean(1 + encode_log_var - K.square(encode_mean) - K.exp(encode_log_var))
    vae_loss = K.mean(xent_loss + kl_loss)
    sampling_model.add_loss(vae_loss)

    # from keras import optimizers
    # adam = optimizers.Adam(lr=0.0001, beta_1=0.9, beta_2=0.999, epsilon=None, decay=3e-7, amsgrad=False)
    sampling_model.compile(optimizer='rmsprop', metrics='mean_squared_error')

    return encoder, sampling_model


def sampling(arg):
    mean = arg[0]
    logvar = arg[1]
    epsilon = K.random_normal(shape=K.shape(mean), mean=0.5, stddev=1.)
    return mean + K.exp(0.5 * logvar) * epsilon


if __name__ == "__main__":
    main()
