import sys
import cooler
import warnings
import numpy as np

warnings.filterwarnings("ignore")

if (sys.argv[1] == "-h") or (sys.argv[1] == "--help"):
    print("\n")
    print("USAGE:")
    print("conda activate hic_processing_env")
    print("python compute_distance_zscore_matrix.py [COOL_PATH] [OUTFILE_MATRIX_ZSCORE] [OUTFILE_BINS]")
    print("\n")
    exit(0)

COOL_PATH = sys.argv[1]
OUTFILE_MATRIX_ZSCORE = sys.argv[2]
OUTFILE_BINS = sys.argv[3]

print("loading Cooler object")
cooler_obj = cooler.Cooler(COOL_PATH)

print("extracting bins and matrix")
cooler_bins = cooler_obj.bins()[:]
cooler_matrix = cooler_obj.matrix(balance = True, sparse = False)
count_matrix_iced = cooler_matrix[:,:]; count_matrix_iced = count_matrix_iced.astype(float)
count_matrix_iced_symm = np.maximum(count_matrix_iced, count_matrix_iced.T)

import pandas as pd
test = count_matrix_iced[10:25,10:25]
test_symm = count_matrix_iced_symm[10:25,10:25]

print(pd.DataFrame(test))
print(np.allclose(test, test.T, equal_nan = True))
print(pd.DataFrame(test_symm))
print(np.allclose(test_symm, test_symm.T, equal_nan = True))

print("computing diagonal mean and std")
n_rows = count_matrix_iced.shape[0]

diag_mean_arr = np.empty(n_rows, dtype = float); diag_std_arr = np.empty(n_rows, dtype = float)

for K in range(n_rows):
    
    diagonal = np.diag(count_matrix_iced_symm, K)
    diag_mean = np.nanmean(diagonal)
    diag_std = np.nanstd(diagonal)

    diag_mean_arr[K] = diag_mean
    diag_std_arr[K] = diag_std

idx = np.arange(n_rows)
distance_matrix = np.abs(idx[:,None] - idx[None,:])

print("computing z-scores")
count_matrix_zscore = (count_matrix_iced_symm - diag_mean_arr[distance_matrix]) / diag_std_arr[distance_matrix]

test_zscore = count_matrix_zscore[10:25,10:25]
print(pd.DataFrame(test_zscore))
print(np.allclose(test_zscore, test_zscore.T, equal_nan = True))

print("saving files")

fmt_matrix = "%.4f"
fmt_bins = ["%s","%.0f","%.0f","%.6f"]

np.savetxt(OUTFILE_MATRIX_ZSCORE, count_matrix_zscore, fmt=fmt_matrix, delimiter="\t")
np.savetxt(OUTFILE_BINS, cooler_bins, fmt=fmt_bins, delimiter="\t")

print("Done!")
exit(0)

