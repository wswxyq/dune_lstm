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
    return np.sqrt(np.mean(x**2))

cuts=[]
lstm_means = []
trans_means = []
lstm_rms = []
trans_rms = []

for x in np.linspace(0., 1.5, 20):
    cut = x
    lstm_resolution = (
        pred_lstm[true_lstm > cut] - true_lstm[true_lstm > cut]
    ) / true_lstm[true_lstm > cut]
    trans_resolution = (
        pred_trans[true_trans > cut] - true_trans[true_trans > cut]
    ) / true_trans[true_trans > cut]
    cuts.append(cut)
    lstm_means.append(np.mean(lstm_resolution))
    trans_means.append(np.mean(trans_resolution))
    lstm_rms.append(rms(lstm_resolution))
    trans_rms.append(rms(trans_resolution))

# %%
plt.figure(dpi=200, figsize=(6, 4))
plt.plot(cuts, lstm_means, label="LSTM")
plt.plot(cuts, trans_means, label="Transformer")
plt.xlabel("True Energy Cut [MeV]")
plt.ylabel("Mean Resolution")
plt.hlines(0, 0, 1.5, linestyles="dashed")
#plt.ylim(-0.5, 0.5)
plt.legend()

plt.show()
# %%
plt.figure(dpi=200, figsize=(6, 4))
plt.plot(cuts, lstm_rms, label="LSTM")
plt.plot(cuts, trans_rms, label="Transformer")
plt.xlabel("True Energy Cut [MeV]")
plt.ylabel("Resolution RMS")
plt.ylim(-0., 0.5)
plt.legend()
plt.show()
# %%
cuts=[]
lstm_means = []
trans_means = []
lstm_rms = []
trans_rms = []

for x in np.linspace(0.5, 5, 50):
    cut = x
    lstm_resolution = (
        pred_lstm[true_lstm < cut] - true_lstm[true_lstm < cut]
    ) / true_lstm[true_lstm < cut]
    trans_resolution = (
        pred_trans[true_trans < cut] - true_trans[true_trans < cut]
    ) / true_trans[true_trans < cut]
    cuts.append(cut)
    lstm_means.append(np.mean(lstm_resolution))
    trans_means.append(np.mean(trans_resolution))
    lstm_rms.append(rms(lstm_resolution))
    trans_rms.append(rms(trans_resolution))

# %%
plt.figure(dpi=200, figsize=(6, 4))
plt.plot(cuts, lstm_means, label="LSTM")
plt.plot(cuts, trans_means, label="Transformer")
plt.xlabel("True Energy Cut [MeV]")
plt.ylabel("Mean Resolution")
plt.hlines(0, 0.5, 5, linestyles="dashed")
#plt.ylim(-0.5, 0.5)
plt.legend()
plt.show()

# %%
plt.figure(dpi=200, figsize=(6, 4))
plt.plot(cuts, lstm_rms, label="LSTM")
plt.plot(cuts, trans_rms, label="Transformer")
plt.xlabel("True Energy Cut [MeV]")
plt.ylabel("Resolution RMS")
plt.ylim(0., 5)
plt.legend()
plt.show()

# %%
