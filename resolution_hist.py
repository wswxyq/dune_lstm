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
pred_lstm = data["lstm_EE_pred.total"]
true_lstm = data["mc.nuE"]

# %%
trans = np.load("transformer_ee_test.npz")
pred_trans = trans["prediction"][:, 0]
true_trans = trans["trueval"][:, 0]
# %%
def rms(x):
    return np.sqrt(np.mean(x ** 2))
plt.figure(dpi=200, figsize=(6, 4))
cut=1.5
lstm_resolution = (pred_lstm - true_lstm) / true_lstm
trans_resolution = (pred_trans - true_trans) / true_trans
reduced_trans_resolution = (pred_trans[true_trans>cut] - true_trans[true_trans>cut]) / true_trans[true_trans>cut]
reduced_lstm_resolution = (pred_lstm[true_lstm>cut] - true_lstm[true_lstm>cut]) / true_lstm[true_lstm>cut]
reduced_trans_resolution1 = (pred_trans[true_trans<cut] - true_trans[true_trans<cut]) / true_trans[true_trans<cut]
reduced_lstm_resolution1 = (pred_lstm[true_lstm<cut] - true_lstm[true_lstm<cut]) / true_lstm[true_lstm<cut]
plt.hist(
    trans_resolution,
    bins=np.linspace(-2, 2, 100),
    histtype="step",
    label="Transformer",
)
plt.hist(lstm_resolution, bins=np.linspace(-2, 2, 100), histtype="step", label="LSTM")
plt.hist(reduced_trans_resolution, bins=np.linspace(-2, 2, 100), histtype="step", label="Transformer (trueE>"+str(cut)+")")
plt.hist(reduced_lstm_resolution, bins=np.linspace(-2, 2, 100), histtype="step", label="LSTM (trueE>"+str(cut)+")")
plt.title("Resolution NC")
plt.legend()
plt.savefig("resolution_high.png")
plt.show()
plt.hist(reduced_trans_resolution1, bins=np.linspace(-3, 3, 100), histtype="step", label="Transformer (trueE<"+str(cut)+")")
plt.hist(reduced_lstm_resolution1, bins=np.linspace(-3, 3, 100), histtype="step", label="LSTM (trueE<"+str(cut)+")")
print("LSTM_resolution\t", "Transformer_resolution\t", "reduced_Transformer_resolution\t", "reduced_LSTM_resolution")
print("mean\t", np.mean(lstm_resolution), "\t", np.mean(trans_resolution), "\t", np.mean(reduced_trans_resolution), "\t", np.mean(reduced_lstm_resolution))
print("std\t", np.std(lstm_resolution), "\t", np.std(trans_resolution), "\t", np.std(reduced_trans_resolution), "\t", np.std(reduced_lstm_resolution))
print("rms\t", rms(lstm_resolution), "\t", rms(trans_resolution), "\t", rms(reduced_trans_resolution), "\t", rms(reduced_lstm_resolution))
plt.legend()
#plt.xlim(-1, 1)
plt.title("Resolution NC")
plt.savefig("resolution_low.png")
plt.show()
# %%
plt.figure()
plt.hist(true_lstm, bins=100, histtype="step", label="TrueE")
#plt.show()
# %%
#plt.hist(true_trans, bins=100, histtype="step", label="TrueE")
# %%
