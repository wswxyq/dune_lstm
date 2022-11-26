# %%
import pandas as pd
import numpy as np
import pickle
import matplotlib.pyplot as plt
from utils.binstat import plot_xstat, plot_2d_hist_count
from utils.fs_cut import get_str_max

data = pd.read_csv(
    "data/kEcut20MeV_01e8maxtraceE/2022-08-21_rnne_NC_250_fGScatter_20MeV_KE_01e8_max_trackE_cut_test_with_pred.csv.xz"
)
# %%
data.columns
hadE=data['mc.nuE']-data['mc.lepE']
# %%
predt = data["lstm_EE_pred.total"]
truet = data["mc.nuE"]

# %%
hadcut=1.5
plot_2d_hist_count(truet[hadE>hadcut], predt[hadE>hadcut], bins=50, range=((0, 5), (0, 5)), name="tmp", ext="pdf", xlabel="True E", ylabel="Reco E")
# %%
plt.hist(hadE)
# %%
plot_2d_hist_count(hadE, predt, bins=50, range=((0, 5), (0, 5)), name="tmp", ext="pdf", xlabel="True hadE", ylabel="Reco E")
# %%
