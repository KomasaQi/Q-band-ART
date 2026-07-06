from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


OUT_DIR = Path(__file__).resolve().parent
NL = "\n"


plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "mathtext.fontset": "dejavusans",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


fig, ax = plt.subplots(figsize=(15.8, 5.85), dpi=260)
ax.set_xlim(0, 15.8)
ax.set_ylim(0, 5.85)
ax.axis("off")

blue = "#DDECFB"
blue_edge = "#2C6FA6"
green = "#E5F2E8"
green_edge = "#3A7E51"
gray = "#F3F5F7"
gray_edge = "#6B747D"
yellow = "#F8EDC8"
yellow_edge = "#A97816"
ink = "#1F2933"
muted = "#5F6974"


def lane(y, h, face, edge, text):
    ax.add_patch(
        FancyBboxPatch(
            (0.25, y),
            15.3,
            h,
            boxstyle="round,pad=0.02,rounding_size=0.07",
            linewidth=1.0,
            edgecolor=edge,
            facecolor=face,
            alpha=0.45,
            zorder=0,
        )
    )
    ax.text(
        0.47,
        y + h - 0.15,
        text,
        fontsize=10.5,
        fontweight="bold",
        color=edge,
        va="top",
        zorder=1,
    )


def box(x, y, w, h, text, face, edge, lw=1.35, fs=9.5, weight="normal", dashed=False):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.075,rounding_size=0.07",
        linewidth=lw,
        edgecolor=edge,
        facecolor=face,
        linestyle="--" if dashed else "-",
        zorder=3,
    )
    ax.add_patch(patch)
    ax.text(
        x + w / 2,
        y + h / 2,
        text,
        ha="center",
        va="center",
        fontsize=fs,
        color=ink,
        fontweight=weight,
        linespacing=1.18,
        zorder=4,
    )
    return patch


def arrow(xy1, xy2, color=gray_edge, lw=1.4, style="-|>", rad=0.0, dashed=False, z=2):
    arr = FancyArrowPatch(
        xy1,
        xy2,
        arrowstyle=style,
        mutation_scale=12,
        linewidth=lw,
        color=color,
        connectionstyle=f"arc3,rad={rad}",
        linestyle="--" if dashed else "-",
        zorder=z,
    )
    ax.add_patch(arr)
    return arr


def label(x, y, text, fs=8.2, color=muted, ha="center"):
    ax.text(x, y, text, fontsize=fs, color=color, ha=ha, va="center", zorder=5)


lane(4.35, 1.18, green, green_edge, "Offline calibration and theory")
lane(0.48, 3.32, blue, blue_edge, "Online Q-band ART controller")

box(0.86, 4.62, 2.35, 0.54, "Reduced roll-slosh" + NL + "model", green, green_edge, fs=9.0)
box(3.62, 4.62, 2.35, 0.54, "Gain map" + NL + r"$G_{za}(j\omega)$", green, green_edge, fs=9.0)
box(
    6.38,
    4.54,
    2.75,
    0.70,
    "Critical bands" + NL + r"$\mathcal{B}_i,\ \bar E_i,\ \bar G_i$",
    green,
    green_edge,
    fs=9.0,
    weight="bold",
)
box(9.70, 4.54, 2.55, 0.70, "Q-band matrices" + NL + r"$Q_i$", green, green_edge, fs=9.0, weight="bold")
box(12.80, 4.54, 2.30, 0.70, "Safety margins" + NL + r"$\bar E_i,\ s_i$", green, green_edge, fs=9.0)
arrow((3.21, 4.89), (3.62, 4.89), green_edge)
arrow((5.97, 4.89), (6.38, 4.89), green_edge)
arrow((9.13, 4.89), (9.70, 4.89), green_edge)
arrow((12.25, 4.89), (12.80, 4.89), green_edge)

box(
    0.70,
    2.57,
    3.05,
    0.78,
    "Reference preview" + NL + "tractor sensors" + NL + r"$(v_x,r,a_y,\delta)$",
    blue,
    blue_edge,
    fs=9.0,
    weight="bold",
)
box(
    0.70,
    1.55,
    3.05,
    0.56,
    "Hitch-angle estimate" + NL + r"$\hat\gamma \rightarrow m_{eq}, I_{z,eq}$",
    blue,
    blue_edge,
    fs=8.6,
)
box(0.70, 0.92, 3.05, 0.46, r"ESO residuals $\hat d_y,\hat d_r$", blue, blue_edge, fs=8.8)

box(
    4.28,
    1.80,
    2.58,
    1.00,
    "Tractor-side" + NL + "prediction model" + NL + r"$x_{k+1}=f(x_k,u_k,\hat d_k)$",
    blue,
    blue_edge,
    fs=8.55,
    weight="bold",
)
box(7.30, 1.96, 2.28, 0.68, "Predicted" + NL + r"$\mathbf{A}_y$ sequence", yellow, yellow_edge, fs=9.1, weight="bold")
box(
    10.05,
    1.80,
    2.78,
    1.00,
    "Q-band projection" + NL + "and scheduler" + NL + r"$E_i=\mathbf{A}_y^TQ_i\mathbf{A}_y$" + NL + r"$w_i(t),\ s_i$",
    yellow,
    yellow_edge,
    fs=8.15,
    weight="bold",
)
box(
    13.20,
    1.80,
    2.15,
    1.00,
    "NMPC" + NL + "tracking + smoothing" + NL + "band soft constraints",
    blue,
    blue_edge,
    fs=8.55,
    weight="bold",
)
box(13.20, 0.72, 2.15, 0.60, "Semi-trailer" + NL + "tank truck", gray, gray_edge, fs=8.9, weight="bold")

arrow((3.75, 2.95), (4.28, 2.45), blue_edge)
arrow((3.75, 1.83), (4.28, 2.20), blue_edge)
arrow((3.75, 1.15), (4.28, 1.95), blue_edge)
arrow((6.86, 2.30), (7.30, 2.30), yellow_edge)
arrow((9.58, 2.30), (10.05, 2.30), yellow_edge)
arrow((12.83, 2.30), (13.20, 2.30), blue_edge)
arrow((14.28, 1.80), (14.28, 1.32), blue_edge)
label(14.86, 1.52, r"$\delta_{cmd}$", fs=8.8, color=blue_edge)

arrow((10.98, 4.54), (11.25, 3.00), green_edge, lw=1.15, dashed=True, rad=-0.03)
arrow((13.75, 4.54), (12.35, 3.00), green_edge, lw=1.15, dashed=True, rad=0.08)
label(12.28, 3.58, "offline parameters", fs=8.0, color=green_edge)
box(
    7.15,
    3.30,
    3.90,
    0.38,
    r"Sufficient bound: $\|z_{\mathcal{B}}\|_2 \leq \bar G_i\sqrt{\bar E_i}$",
    "#FFFFFF",
    green_edge,
    lw=1.0,
    fs=8.3,
)
arrow((7.75, 4.54), (8.55, 3.68), green_edge, lw=1.0, dashed=True, rad=-0.10)

arrow((13.20, 1.02), (9.75, 1.02), gray_edge, lw=1.2)
arrow((9.75, 1.02), (3.78, 1.02), gray_edge, lw=1.2)
label(8.45, 0.78, "tractor-side feedback", fs=8.2, color=gray_edge)

ax.text(0.28, 5.70, "Q-band ART architecture", fontsize=13.2, fontweight="bold", color=ink, va="center")
ax.text(4.20, 5.70, "input-side critical-band suppression with tractor-side sensing", fontsize=9.3, color=muted, va="center")

for suffix in ("png", "pdf", "svg"):
    fig.savefig(OUT_DIR / f"intro_architecture.{suffix}", bbox_inches="tight", pad_inches=0.04)

plt.close(fig)
