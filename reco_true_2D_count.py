# %%
import pandas as pd
from utils.fs_cut import FSIdic
from utils.binstat import plot_2d_hist_count


# %%
for i in [1, 2, 3, 4, 6, 7, 8, 9]:
    filename = "result"
    data = pd.read_csv(
        "data/kEcut20MeV_01e8maxtraceE/testset_with_FSIcut/FSI_" + str(i) + ".csv.xz"
    )
    basename = filename + str(i)
    plot_2d_hist_count(
        data["mc.nuE"],
        data["lstm_EE_pred.total"],
        xlabel="True Energy",
        ylabel="Reconstructed Energy",
        name=basename + "_reco_true_2D_count",
        title=FSIdic[i],
        outdir="plot/FSIplot01e8",
    )

# %%
data = pd.read_csv(
    "data/kEcut20MeV_01e8maxtraceE/2022-08-21_rnne_NC_250_fGScatter_20MeV_KE_01e8_max_trackE_cut_test_with_pred.csv.xz"
)
plot_2d_hist_count(
    data["mc.nuE"],
    data["lstm_EE_pred.total"],
    xlabel="True Energy",
    ylabel="Reconstructed Energy",
    name="test_reco_true_2D_count",
    title="testset (all FS)",
    outdir="plot/FSIplot01e8",
)

# %%
