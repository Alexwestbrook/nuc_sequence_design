from pathlib import Path

import numpy as np
import tensorflow as tf
from MyModuleLibrary.mykeras.losses import correlate, mae_cor
from numpy.core.numeric import normalize_axis_tuple
from numpy.lib.stride_tricks import as_strided


def read_fasta(filename):
    with open(filename, "r") as f:
        genome = {}
        n_seqs = 0
        sequence = ""
        for line in f:
            if line[0] == ">":  # First line header, discard this line
                # Save sequence of previous chromosome
                if n_seqs >= 1 and sequence != "":
                    genome[id] = sequence
                    sequence = ""
                # Get new chromosome id
                id, *_ = line.split()
                id = id[1:]
                n_seqs += 1
            else:
                sequence += line.rstrip()
    if n_seqs >= 1 and sequence != "":
        genome[id] = sequence
    return genome


def one_hot_encode(seq, length=None, one_hot_type=bool, order="ACGT"):
    if length is None:
        length = len(seq)
    one_hot = np.zeros((length, 4), dtype=one_hot_type)
    for i, base in enumerate(seq):
        if i >= length:
            break
        base = base.upper()
        if base == order[0]:
            one_hot[i, 0] = 1
        elif base == order[1]:
            one_hot[i, 1] = 1
        elif base == order[2]:
            one_hot[i, 2] = 1
        elif base == order[3]:
            one_hot[i, 3] = 1
    return one_hot


def sliding_window_view(x, window_shape, axis=None, *, subok=False, writeable=False):
    """Function from the numpy library"""
    window_shape = tuple(window_shape) if np.iterable(window_shape) else (window_shape,)
    # first convert input to array, possibly keeping subclass
    x = np.array(x, copy=False, subok=subok)

    window_shape_array = np.array(window_shape)
    if np.any(window_shape_array < 0):
        raise ValueError("`window_shape` cannot contain negative values")

    if axis is None:
        axis = tuple(range(x.ndim))
        if len(window_shape) != len(axis):
            raise ValueError(
                f"Since axis is `None`, must provide "
                f"window_shape for all dimensions of `x`; "
                f"got {len(window_shape)} window_shape elements "
                f"and `x.ndim` is {x.ndim}."
            )
    else:
        axis = normalize_axis_tuple(axis, x.ndim, allow_duplicate=True)
        if len(window_shape) != len(axis):
            raise ValueError(
                f"Must provide matching length window_shape and "
                f"axis; got {len(window_shape)} window_shape "
                f"elements and {len(axis)} axes elements."
            )

    out_strides = x.strides + tuple(x.strides[ax] for ax in axis)

    # note: same axis can be windowed repeatedly
    x_shape_trimmed = list(x.shape)
    for ax, dim in zip(axis, window_shape):
        if x_shape_trimmed[ax] < dim:
            raise ValueError("window shape cannot be larger than input array shape")
        x_shape_trimmed[ax] -= dim - 1
    out_shape = tuple(x_shape_trimmed) + window_shape
    return as_strided(
        x, strides=out_strides, shape=out_shape, subok=subok, writeable=writeable
    )


def full_predict(one_hot_chr, model, reverse=False):
    if reverse:
        one_hot_chr = one_hot_chr[::-1, [1, 0, 3, 2]]
    WINDOW = 2001
    side_arr = np.zeros_like(one_hot_chr, shape=(WINDOW // 2, 4))
    one_hot_chr = np.vstack([side_arr, one_hot_chr, side_arr])
    X = sliding_window_view(one_hot_chr, window_shape=(WINDOW, 4)).reshape(
        -1, WINDOW, 4, 1
    )
    pred = model.predict(X).squeeze()
    if reverse:
        pred = pred[::-1]
    return pred


# Read fasta files
sequences = {}
for filename in [
    "167_7_4kbrf.fa",
    "167_601_7_4kbrf.fa",
    "197b_7_4kbrf.fa",
    "197_601_7_4kbrf.fa",
    "237_7_4kbrf.fa",
    "237_601_7_4kbrf.fa",
]:
    sequences.update(read_fasta(Path("..", "genome", filename)))

# Prepare input for model
one_hot_repeats = {}
start = 4000
for seq_id, seq in sequences.items():
    # Get repeat length
    rlen = int(seq_id[:3])
    # Extract repeat sequence
    rep_seq = seq[start : start + rlen]
    # Convert to one_hot
    rep_one_hot = one_hot_encode(rep_seq, order="ATGC")
    # Tile repeats left and right to get a 2000 + 7*rlen sequence
    # centered on the 7 repeats
    q, r = divmod(1000, rlen)
    one_hot = np.tile(rep_one_hot, ((q + 1) * 2 + 7, 1))
    edge = rlen - r
    if edge != 0:
        one_hot = one_hot[edge:-edge]
    one_hot_repeats[seq_id[:-8]] = one_hot

# Load model
model_name = "weights_with_rev_compl_rep2"
model = tf.keras.models.load_model(
    Path("..", "Results_nucleosome", f"{model_name}.hdf5"),
    custom_objects={"mae_cor": mae_cor, "correlate": correlate},
)

# Predict
preds = {}
for k, v in one_hot_repeats.items():
    preds[k] = full_predict(v, model)

# Save
np.savez_compressed(
    Path("..", "Results_nucleosome", f"preds_{model_name}_on_repeats.npz"), **preds
)
