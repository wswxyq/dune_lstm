# %%
import pandas as pd
import numpy as np
import pickle
import matplotlib.pyplot as plt
from utils.binstat import plot_xstat
from utils.fs_cut import get_str_max

data = pd.read_csv(
    "data/dataset_rnne_dune_numu.csv.xz"
)
def string_to_float_list(string):
    """
    By default, convert all strings to list of floats.
    Otherwise, return empty list.
    """
    if not string or not isinstance(string, str):
        return [0]
    return [float(s) for s in string.split(",")]

# %%
c=0
for x in data["particle.dir.z"]:
    if not isinstance(x, str):
        c+=1
print(c)
# %%
newdata=data[data["particle.dir.z"].apply(lambda x: isinstance(x, str))]

# %%
newdata.to_csv("data/new_dataset_rnne_dune_numu.csv.xz", index=False, compression='xz')
# %%
c=0
for x in newdata["particle.dir.z"]:
    if not isinstance(x, str):
        c+=1
print(c)
# %%
data.to_csv("data/test.csv.xz", index=False, compression='xz')

# %%diff 
len(data)
# %%
len(newdata)
# %%
