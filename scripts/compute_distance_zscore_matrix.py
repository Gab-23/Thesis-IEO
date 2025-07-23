import sys

if (sys.argv[1] == "-h") or (sys.argv[1] == "--help"):
    print("\n")
    print("""
    This script is to be run after merge_and_normalize_cooler.sh
    It takes the ICE normalized cooler and applies Z-score normalization, given a fixed distance
    The idea is to remove distance-decay effect.
            """)
    print("\n")
    print("USAGE:")
    print("\n")
    print("conda activate hic_processing_env")
    print("python compute_distance_zscore_matrix.py [COOL_PATH] [OUTFILE_MATRIX_ZSCORE] [OUTFILE_BINS]")
    print("\n")
    exit(0)

import cooler
import warnings
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

COOL_PATH = sys.argv[1]
OUTFILE_MATRIX_ZSCORE = sys.argv[2]
OUTFILE_BINS = sys.argv[3]

print("loading Cooler object")
cooler_obj = cooler.Cooler(COOL_PATH)

print("extracting bins and matrix")
cooler_bins = cooler_obj.bins()[:] # extract bins
cooler_matrix = cooler_obj.matrix(balance = True, sparse = False) # take the balanced (balance = True) matrix
count_matrix_iced = cooler_matrix[:,:]; count_matrix_iced = count_matrix_iced.astype(float)
count_matrix_iced_symm = np.maximum(count_matrix_iced, count_matrix_iced.T) # symmetrize the matrix

test = count_matrix_iced[10:25,10:25]
test_symm = count_matrix_iced_symm[10:25,10:25]

print(pd.DataFrame(test))
print(np.allclose(test, test.T, equal_nan = True))
print(pd.DataFrame(test_symm))
print(np.allclose(test_symm, test_symm.T, equal_nan = True))

print("computing diagonal mean and std")
n_rows = count_matrix_iced.shape[0]

diag_mean_arr = np.empty(n_rows, dtype = float); diag_std_arr = np.empty(n_rows, dtype = float)

for K in range(n_rows): # iterate over diagonals
    
    diagonal = np.diag(count_matrix_iced_symm, K)
    diag_mean = np.nanmean(diagonal) # take the mean
    diag_std = np.nanstd(diagonal) # take the sd

    diag_mean_arr[K] = diag_mean
    diag_std_arr[K] = diag_std

idx = np.arange(n_rows)
distance_matrix = np.abs(idx[:,None] - idx[None,:])

# the idea is that idx[:,None] makes a row vector with idx, idx[None,:] makes a column vector with idx and the subtraction gives
# a symmetrical matrix where each cell represents the distance between two bins

print("computing z-scores")
count_matrix_zscore = (count_matrix_iced_symm - diag_mean_arr[distance_matrix]) / diag_std_arr[distance_matrix] # compute Z-score normalization

test_zscore = count_matrix_zscore[10:25,10:25]
print(pd.DataFrame(test_zscore))
print(np.allclose(test_zscore, test_zscore.T, equal_nan = True))

print("saving files")

fmt_matrix = "%.4f"
fmt_bins = ["%s","%.0f","%.0f","%.6f"]

np.savetxt(OUTFILE_MATRIX_ZSCORE, count_matrix_zscore, fmt=fmt_matrix, delimiter="\t") # save the output
np.savetxt(OUTFILE_BINS, cooler_bins, fmt=fmt_bins, delimiter="\t") # save the output

print("Done!")
exit(0)

