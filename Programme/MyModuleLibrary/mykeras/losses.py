# -*- coding: utf-8 -*-

"""
This module contains the custom losses or metrics that can be used to
train or to evaluate a neural network.
It is made to work as a usual loss or metric.
"""

try:
    from keras import backend as K
except ModuleNotFoundError:
    from tensorflow.keras import backend as K


def correlate(y_true, y_pred):
    """
    Calculate the correlation between the predictions and the labels.

    :Example:

    >>> model.compile(optimizer = 'adam', losses = correlate)
    >>> load_model('file', custom_objects = {'correlate : correlate})
    """
    X = y_true - K.mean(y_true)
    Y = y_pred - K.mean(y_pred)

    sigma_XY = K.sum(X * Y)
    sigma_X = K.sqrt(K.sum(X * X))
    sigma_Y = K.sqrt(K.sum(Y * Y))

    return sigma_XY / (sigma_X * sigma_Y + K.epsilon())


def mae_cor(y_true, y_pred):
    """
    Calculate the mean absolute error minus the correlation between
    predictions and  labels.

    :Example:

    >>> model.compile(optimizer = 'adam', losses = mae_cor)
    >>> load_model('file', custom_objects = {'mae_cor : mae_cor})
    """
    X = y_true - K.mean(y_true)
    Y = y_pred - K.mean(y_pred)

    sigma_XY = K.sum(X * Y)
    sigma_X = K.sqrt(K.sum(X * X))
    sigma_Y = K.sqrt(K.sum(Y * Y))

    cor = sigma_XY / (sigma_X * sigma_Y + K.epsilon())
    mae = K.mean(K.abs(y_true - y_pred))

    alpha = 1.0
    beta = 1.0

    return beta * mae + alpha * (1 - cor)


def mse_var(y_true, y_pred):
    """
    Calculate the mean squared error between the predictions and the
    labels and add the absolute difference of variance between the
    distribution of labels and the distribution of predictions.

    :Example:

    >>> model.compile(optimizer = 'adam', losses = mse_var)
    >>> load_model('file', custom_objects = {'mse_var' : mse_var})
    """
    X = y_true - y_pred

    Y = K.mean(X**2) + K.abs(K.var(y_true) - K.var(y_pred))

    return Y
