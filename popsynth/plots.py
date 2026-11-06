"""Standard figures for the prototype comparison."""
import matplotlib.pyplot as plt
import numpy as np

from .analysis import COLORS, LABELS


def period_distributions(res, obs=None):
    """Left: intrinsic present-day periods. Right: detected dark-companion orbits."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    ax = axes[0]
    bins = np.linspace(-1, 6, 57)
    for t, r in res.items():
        lp = np.log10(r["pop"].porb)
        ax.hist(lp, bins=bins, histtype="step", lw=2, color=COLORS.get(t), label=LABELS.get(t, t),
                weights=np.full(len(lp), r["mock"].attrs["weight"]))
    ax.axvspan(2, 3, color="0.85", zorder=0, label="DR3 astrometric window (100-1000 d)")
    ax.set_yscale("log")
    ax.set_xlabel("log$_{10}$ P$_{orb}$ [d]")
    ax.set_ylabel("WD + MS binaries within 2 kpc (all, intrinsic)")
    ax.set_title("Intrinsic present-day WD + MS periods")
    ax.legend(fontsize=8, loc="upper left")

    ax = axes[1]
    bins2 = np.linspace(1, 4, 31)
    for t, r in res.items():
        for rel, ls in (("DR3", "-"), ("DR4", "--")):
            d = r[rel][r[rel].wd_class]
            ax.hist(np.log10(d.porb), bins=bins2, histtype="step", lw=2, ls=ls, color=COLORS.get(t),
                    label=f"{LABELS.get(t, t)}, {rel}")
    if obs is not None:
        ax.hist(np.log10(obs.period), bins=bins2, histtype="stepfilled", color="0.6", alpha=0.6,
                label=f"Gaia DR3 observed (triage), N = {len(obs)}")
    ax.axvspan(2, 3, color="0.9", zorder=0)
    ax.set_yscale("log")
    ax.set_xlabel("log$_{10}$ P$_{orb}$ [d]")
    ax.set_ylabel("detected dark-companion orbits")
    ax.set_title("After simplified Gaia selection + AMRF triage")
    ax.legend(fontsize=7, loc="upper left")
    return fig


def period_vs_wd_mass(res):
    fig, axes = plt.subplots(1, len(res), figsize=(5.2 * len(res), 4.2), constrained_layout=True,
                             squeeze=False)
    for ax, (t, r) in zip(axes[0], res.items()):
        for rel, c in (("DR4", "0.75"), ("DR3", COLORS.get(t))):
            d = r[rel][r[rel].wd_class]
            ax.scatter(np.log10(d.porb), d.m_wd, s=3, alpha=0.4, color=c, label=f"{rel} (N = {len(d)})",
                       rasterized=True)
        ax.axvspan(2, 3, color="0.92", zorder=0)
        ax.set_xlim(0.5, 4)
        ax.set_ylim(0.2, 1.4)
        ax.set_xlabel("log$_{10}$ P$_{orb}$ [d]")
        ax.set_ylabel("M$_{WD}$ [M$_\\odot$]")
        ax.set_title(LABELS.get(t, t))
        ax.legend(fontsize=8, markerscale=4)
    return fig
