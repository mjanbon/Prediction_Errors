import h5py
import sys

def inspect(path):
    print('Inspecting', path)
    with h5py.File(path,'r') as f:
        def visitor(name,obj):
            if isinstance(obj, h5py.Dataset):
                print(name, obj.shape, obj.dtype)
        f.visititems(visitor)

if __name__=='__main__':
    inspect(r'C:\Users\maxim\OneDrive\Desktop\PhD\Juho Code\Co-I Pipeline\CoI-pipeline\Data\MI_Data\tmp_R060721\sleep_window_scan\R060721_BSLEEP_sleep_window_MI_summary.mat')
    print('---')
    inspect(r'C:\Users\maxim\OneDrive\Desktop\PhD\Juho Code\Co-I Pipeline\CoI-pipeline\Data\MI_Data\R060721sleepBSLEEP_MI_data.mat')
