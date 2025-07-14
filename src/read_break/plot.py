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
    plot_base_totals_only(base_array, ax, small_marker, marker, draw_sorted_lines=True, normed=normed)
    highlight_base_totals(base_array, bit_filter, ax, bit_marker, marker, normed=normed, mfc='None', mec='black', mew=1)
    ax.legend()
    return()


def plot_base_totals_only(base_array, ax, marker_size, marker='.', draw_sorted_lines=True, normed=False):
    """
    Plot the total values for each base in the given array.
    """
    ans = np.array([np.bincount(position + 1, minlength=6)[1:-1] for position in base_array.T])  # Exclude -1 values

    if normed:
        ans = ans / np.maximum(1, np.sum(ans, axis=1)[:, np.newaxis])
    if draw_sorted_lines:
        sorted_ans = np.sort(ans, axis=1)
        ax.plot(sorted_ans, '-k', alpha=0.5)
    ax.plot(ans, marker, ms=marker_size, label=list('ACGT'))  # Assuming BASES[:4] means 'A', 'C', 'G', 'T'
    return()

def highlight_base_totals(base_array, highlight_filter, ax, marker_size, marker='.', normed=False, *, mfc='None', mec='black', mew=1, label=None):
    """
    Highlight specific bases in the base totals plot.
    
    Args:
        base_array (np.ndarray): 2D array of base calls.
        highlight_filter (np.ndarray): Boolean array indicating which bases to highlight.
        ax (matplotlib.axes.Axes): Axes object to plot on.
        marker_size (int): Size of the markers for highlighted bases.
        marker (str): Marker style for highlighted bases.
        normed (bool): Whether to normalize the values.
    """
    ans = np.array([np.bincount(position + 1, minlength=6)[1:-1] for position in base_array.T])  # Exclude -1 values

    if normed:
        ans = ans / np.maximum(1, np.sum(ans, axis=1)[:, np.newaxis])

    ax.plot(np.where(highlight_filter)[0], ans[highlight_filter], marker, ms=marker_size, mfc=mfc, mec=mec, mew=mew, label=label)
    return()


