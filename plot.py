import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --- load your data ---
df = pd.read_csv("PCA.csv")  # change path as needed

# pick the PCs to plot
xcol = "PCA|All|PC1 (41.6%)"
ycol = "PCA|All|PC2 (18.2%)"

x = df[xcol].to_numpy()
y = df[ycol].to_numpy()

# if you want to derive labels/groups from the file instead of hardcoding:
labels = df["General|All|label"].tolist()
groups = df["General|All|group"].tolist()

colors = {"IGF": "tab:green", "Rapamycin": "tab:red", "Control": "tab:blue", "QC": "tab:gray"}

fig, ax = plt.subplots(figsize=(7, 6))

# plot each group with its own color
for g in sorted(set(groups)):
    idx = [i for i, gg in enumerate(groups) if gg == g]
    ax.scatter(x[idx], y[idx], label=g, s=60, color=colors.get(g, "tab:cyan"))

    # centroid circle (radius = mean distance to centroid)
    cx = float(np.mean(x[idx]))
    cy = float(np.mean(y[idx]))
    r = float(np.mean(np.hypot(x[idx] - cx, y[idx] - cy)))
    circ = plt.Circle((cx, cy), r, fill=False, lw=1.5, color=colors.get(g, "tab:cyan"))
    ax.add_patch(circ)

# point labels (optional; comment out if cluttered)
for xi, yi, lab in zip(x, y, labels):
    ax.annotate(lab, (xi, yi), fontsize=8, xytext=(4, 3), textcoords="offset points")

ax.set_xlabel(xcol)
ax.set_ylabel(ycol)
ax.set_title("PCA: PC1 vs PC2")
ax.legend(title="Group", frameon=False)
ax.grid(True, alpha=0.3)

# Show a window (interactive) OR save to file for git/CI
# plt.show()
plt.tight_layout()
plt.savefig("fig_pca_PC1_PC2.png", dpi=300)
print("Saved plot to fig_pca_PC1_PC2.png")