# %%
import pandas as pd
import matplotlib.pyplot as plt
from utils.fs_cut import FSIdic
import os

# %%
for i in [1, 2, 3, 4, 6, 7, 8, 9]:
    filename = "result"
    data = pd.read_csv(
        "data/kEcut20MeV_01e8maxtraceE/testset_with_FSIcut/FSI_" + str(i) + ".csv.xz"
    )
    basename = filename + str(i)

    _fig, _ax = plt.subplots(
        1, 1, figsize=(8, 6), dpi=80
    )

    _ax.hist(
        data['mc.nuE'],
        bins=50,
        range=(0, 5),
        histtype="step",
    )
    _ax.set_xlabel("True Energy", fontsize=14)
    _ax.set_title(FSIdic[i], fontsize=18)
    _ax.text(
        0.95,
        0.95,
        f"NumEvents : {len(data['mc.nuE']):.3e}\n",
        transform=_ax.transAxes,
        ha="right",
        va="top",
    )
    plt.savefig(
        os.path.join("plot/FSIplot01e8", filename + str(i) + "_trueE_hist" + "." + "pdf")
    )

# %%
