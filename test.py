# %%
import pandas as pd
import numpy as np
import pickle
import matplotlib.pyplot as plt
from utils.binstat import plot_xstat
from utils.fs_cut import get_str_max

data = pd.read_csv(
    "data/kEcut20MeV_01e8maxtraceE/2022-08-21_rnne_NC_250_fGScatter_20MeV_KE_01e8_max_trackE_cut_test_with_pred.csv.xz"
)

# %%
predt = data["lstm_EE_pred.total"]
truet = data["mc.nuE"]

# %%
trans = np.load("flat.npz")

# %%
basename = "tmp"
d = plot_xstat(
    truet,
    (predt - truet) / truet,
    bins=50,
    range=(0, 5),
    name=basename + "_xstat",
    ext="pdf",
    xlabel="True E",
    ylabel="(Reco E - True E) / True E",
)
# %%
# %%
basename = "tmp"
dt = plot_xstat(
    trans["trueval"][:, 0],
    (trans["prediction"][:, 0] - trans["trueval"][:, 0]) / trans["trueval"][:, 0],
    bins=50,
    range=(0, 5),
    name=basename + "_xstat",
    ext="pdf",
    xlabel="True E",
    ylabel="(Reco E - True E) / True E",
)
# %%
plt.close("all")
plt.fill_between(d['bin_x'], d['bin_y'] - d['bin_yerr'], d['bin_y'] + d['bin_yerr'], alpha=0.5, label='lstm')
plt.fill_between(dt['bin_x'], dt['bin_y'] - dt['bin_yerr'], dt['bin_y'] + dt['bin_yerr'], alpha=0.5, label='transformer')
plt.legend()
plt.xlabel("True E")
plt.ylabel("(Reco E - True E) / True E")
plt.show()
# %%
