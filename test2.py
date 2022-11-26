# %%
import pandas as pd
import numpy as np
import pickle
import matplotlib.pyplot as plt
from utils.binstat import plot_xstat

data = pd.read_csv('data/kEcut20MeV_01e8maxtraceE/2022-08-21_rnne_NC_250_fGScatter_20MeV_KE_01e8_max_trackE_cut.csv.xz')

# %%
def string_to_float_list(string):
    """
    By default, convert all strings to list of floats.
    Otherwise, return empty list.
    """
    if not string or not isinstance(string, str):
        return []
    return [float(s) for s in string.split(",")]

# %%
len(np.concatenate(data['mc.fpdgCode'].apply(string_to_float_list)))
# %%
data['mc.fpdgCode']
# %%
len(set(np.concatenate(data['mc.fpdgCode'].apply(string_to_float_list))))
# %%
set(data['mc.fpdgCode'])