import matplotlib.pyplot as plt
import numpy as np

'''
Many plots work on a grid of integers of the form:
each row is a read and the columns are the positions in the read.
The values in the grid are the base calls, encoded as integers.
If there is no value present, we typically denote that with a -1.
'''

'''
This function plots the column sums
in a way that makes it easier to determine the major base calls
and what are the major changes.

We can also select a set of positions to highlight.
'''

def plot_base_totals(base_array, bit_filter, ax, small_marker=20, bit_marker=30, marker='.', normed=False):
    """
    base_arrays: list of 2D numpy arrays (one per plot)
    bit_filters: list of boolean arrays (one per plot)
    axes: list/array of matplotlib Axes objects
    """
    ms = small_marker
    ms2 = bit_marker    
    x_bits = np.where(bit_filter)[0]
    ans = np.array([np.bincount(position + 1, minlength=6)[1:-1] for position in base_array.T]) # Exclude -1 values
    
    if normed:
        ans = ans / np.maximum(1, np.sum(ans, axis=1)[:, np.newaxis])
    
    sorted_ans = np.sort(ans, axis=1)
    
    ax.plot(sorted_ans, '-k', alpha=0.5)
    ax.plot(ans, marker, ms=ms, label=list('ACGT'))  # Assuming BASES[:4] means 'A', 'C', 'G', 'T'
    ax.plot(x_bits, ans[bit_filter], marker, ms=ms2, mfc='None', mec="black", mew=1)
    
    ax.legend()
    return()
