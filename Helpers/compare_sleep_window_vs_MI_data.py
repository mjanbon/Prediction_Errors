import scipy.io as sio
import numpy as np
import sys

f1 = r'C:\Users\maxim\OneDrive\Desktop\PhD\Juho Code\Co-I Pipeline\CoI-pipeline\Data\MI_Data\R060721sleepBSLEEP_MI_data.mat'
f2 = r'C:\Users\maxim\OneDrive\Desktop\PhD\Juho Code\Co-I Pipeline\CoI-pipeline\Data\MI_Data\tmp_R060721\sleep_window_scan\R060721_BSLEEP_sleep_window_MI_summary.mat'

print('Loading', f2)
d2 = sio.loadmat(f2, struct_as_record=False, squeeze_me=True)
keys2 = [k for k in d2.keys() if not k.startswith('__')]
print('keys in summary mat:', keys2)
summary = d2.get('summary', None)
if summary is None:
    print('No variable named summary in', f2)
    sys.exit(1)
wr = summary.window_results
means = []
for w in wr:
    si = getattr(w, 'sleep_MI')
    a = np.asarray(si, dtype=float)
    means.append(np.nanmean(a))
slide_mean_mean = float(np.nanmean(means))
print('Sliding windows: mean of mean(sleep_MI) across windows =', slide_mean_mean)

print('\nLoading', f1)
d1 = sio.loadmat(f1, struct_as_record=False, squeeze_me=True)
keys1 = [k for k in d1.keys() if not k.startswith('__')]
print('keys in MI data mat:', keys1)
# try to find numeric variables that look like MI
candidates = []
for k in keys1:
    v = d1[k]
    if isinstance(v, (np.ndarray, list)):
        try:
            arr = np.asarray(v, dtype=float)
            candidates.append((k, arr))
        except:
            pass

if not candidates:
    print('No numeric arrays found to compare in', f1)
    sys.exit(0)

# Choose the largest 1D/2D numeric array as MI candidate
candidates.sort(key=lambda x: x[1].size, reverse=True)
for k, arr in candidates[:5]:
    print('Candidate', k, 'shape', arr.shape, 'mean', np.nanmean(arr))

# Print top candidate mean
k, arr = candidates[0]
print('\nUsing', k, 'for comparison. Mean =', float(np.nanmean(arr)))
print('Sliding windows mean vs file mean: ', slide_mean_mean, 'vs', float(np.nanmean(arr)))

# Save results
out = {'slide_mean_mean': slide_mean_mean, 'file_var': k, 'file_mean': float(np.nanmean(arr))}
np.save('CoI-pipeline/Helpers/compare_sleep_window_vs_MI_data_result.npy', out)
print('Saved summary to CoI-pipeline/Helpers/compare_sleep_window_vs_MI_data_result.npy')
