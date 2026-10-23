#!/usr/bin/env python3
"""orfarena benchmark: one metrics table and ROC curves per species.

Input: scores/<species>.scores.tsv.gz, one row per transcript (key, label, <tool>_score):
  * label 1, annotated coding transcript (key chr:cds_start-cds_end@transcript): the tool's
    score for its prediction of exactly the annotated CDS (genomic span), 0 without one;
  * label 0, lncRNA without repeats (key transcript): the tool's best score on the
    transcript, 0 without a prediction.

Output: arena.tsv (TP/TN/FP/FN at each tool's threshold, sensitivity, specificity,
precision, F1, balanced accuracy, ROC AUC, average precision) and roc/<species>.{png,svg},
roc/all_species.{png,svg}.

usage: python arena.py [--scores scores] [--out .] [--font-path Arial.ttf]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib import font_manager  # noqa: E402
from matplotlib.ticker import MultipleLocator  # noqa: E402

SPECIES = {  # display order
    "HLcanLupFam7": "Dog (HLcanLupFam7)",
    "HLbosTau10": "Cow (HLbosTau10)",
    "HLratNor8": "Rat (HLratNor8)†",
    "HLeleMaxInd3A": "Asian elephant (HLeleMaxInd3A)",
    "HLornAna4": "Platypus (HLornAna4)",
}
# tool -> (label, threshold, color); colors are categorical slots 1-8 in fixed order,
# so a tool keeps its color in every figure
TOOLS = {
    "xorf": ("xORF v0.2.0", 0.29426348209381104, "#2a78d6"),
    "translationai": ("TranslationAI", 0.5, "#eb6834"),
    "td2_strict": ("TD2 strict", 0.5, "#1baf7a"),
    "transaid": ("TRANSAID", 0.635, "#eda100"),
    "td2": ("TD2", 0.5, "#e87ba4"),
    "psauron": ("PSAURON", 0.5, "#008300"),
    "transdecoder": ("TransDecoder", 0.5, "#6250d6"),
    "genemark": ("GeneMark", 0.5, "#e34948"),
}
SURFACE, INK, INK_2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
ZOOM = ((0.0, 0.2), (0.8, 1.0))  # FPR, TPR
FOOTNOTE = "† xORF v0.2.0 was trained with rat: its rat scores are gene-grouped 5-fold out-of-fold predictions."


def roc(y: np.ndarray, s: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Step ROC over the distinct scores (ties move together), trapezoid AUC."""
    order = np.argsort(-s, kind="stable")
    y, s = y[order], s[order]
    cut = np.r_[np.where(np.diff(s))[0], len(s) - 1]
    tps = np.cumsum(y)[cut]
    fps = 1 + cut - tps
    tpr, fpr = np.r_[0, tps / y.sum()], np.r_[0, fps / (len(y) - y.sum())]
    return fpr, tpr, float(np.trapezoid(tpr, fpr))


def average_precision(y: np.ndarray, s: np.ndarray) -> float:
    order = np.argsort(-s, kind="stable")
    y, s = y[order], s[order]
    cut = np.r_[np.where(np.diff(s))[0], len(s) - 1]
    tps = np.cumsum(y)[cut]
    fps = 1 + cut - tps
    recall, precision = tps / y.sum(), tps / (tps + fps)
    return float(np.sum((recall - np.r_[0, recall[:-1]]) * precision))


def metrics(y: np.ndarray, s: np.ndarray, threshold: float) -> dict:
    p = s >= threshold
    tp, fp = int((y & p).sum()), int((~y & p).sum())
    tn, fn = int((~y & ~p).sum()), int((y & ~p).sum())
    sens, spec = tp / (tp + fn), tn / (tn + fp)
    prec = tp / (tp + fp) if tp + fp else np.nan
    return {"threshold": threshold, "TP": tp, "TN": tn, "FP": fp, "FN": fn, "sensitivity": sens, "specificity": spec,
            "precision": prec, "F1": 2 * prec * sens / (prec + sens), "balanced_accuracy": (sens + spec) / 2,
            "ROC_AUC": roc(y.astype(int), s)[2], "average_precision": average_precision(y.astype(int), s)}


def style(ax, zoom: bool) -> None:
    ax.set_facecolor(SURFACE)
    (x0, x1), (y0, y1) = ZOOM if zoom else ((0, 1), (0, 1))
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.xaxis.set_major_locator(MultipleLocator(0.05 if zoom else 0.2))
    ax.yaxis.set_major_locator(MultipleLocator(0.05 if zoom else 0.2))
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(colors=MUTED, labelcolor=INK_2, labelsize=8, width=0.8, length=3)
    ax.set_xlabel("False positive rate (1 − specificity)", color=INK_2, fontsize=9)
    ax.set_ylabel("True positive rate (sensitivity)", color=INK_2, fontsize=9)
    ax.plot([0, 1], [0, 1], color=AXIS, linewidth=0.8, zorder=1)  # chance


def curves(ax, df: pd.DataFrame, aucs: dict) -> dict:
    """Draw every tool's ROC; returns tool -> line, in descending AUC order."""
    y = df["label"].to_numpy().astype(int)
    handles = {}
    # draw xORF last so it sits on top; the legend keeps the AUC ranking
    for tool in sorted(TOOLS, key=lambda t: t == "xorf"):
        label, _, color = TOOLS[tool]
        fpr, tpr, _ = roc(y, df[f"{tool}_score"].to_numpy())
        (h,) = ax.plot(fpr, tpr, color=color, linewidth=1.5, solid_joinstyle="round", solid_capstyle="round",
                       zorder=3 if tool == "xorf" else 2, label=f"{label}  {aucs[tool]:.3f}")
        handles[tool] = h
    return dict(sorted(handles.items(), key=lambda kv: -aucs[kv[0]]))


def legend(fig_or_ax, handles, **kw):
    leg = fig_or_ax.legend(handles=handles, frameon=False, fontsize=8, labelcolor=INK, handlelength=1.6, **kw)
    return leg


def main() -> None:
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scores", type=Path, default=here / "scores")
    ap.add_argument("--out", type=Path, default=here)
    ap.add_argument("--font-path", type=Path, help="Arial TTF, if matplotlib does not find Arial by name")
    args = ap.parse_args()

    if args.font_path:
        font_manager.fontManager.addfont(str(args.font_path))
    plt.rcParams.update({"font.family": "Arial", "svg.fonttype": "none", "figure.facecolor": SURFACE,
                         "savefig.facecolor": SURFACE})
    (args.out / "roc").mkdir(parents=True, exist_ok=True)

    rows, data = [], {}
    for sp, name in SPECIES.items():
        df = pd.read_csv(args.scores / f"{sp}.scores.tsv.gz", sep="\t")
        data[sp] = df
        y = df["label"].to_numpy().astype(bool)
        for tool, (label, thr, _) in TOOLS.items():
            rows.append({"species": sp, "tool": label, "positives": int(y.sum()), "negatives": int((~y).sum()),
                         **metrics(y, df[f"{tool}_score"].to_numpy(), thr)})
    table = pd.DataFrame(rows)
    table.to_csv(args.out / "arena.tsv", sep="\t", index=False, float_format="%.4f")

    aucs = {sp: dict(zip(TOOLS, table[table.species == sp]["ROC_AUC"])) for sp in SPECIES}
    for sp, name in SPECIES.items():
        fig, axes = plt.subplots(1, 2, figsize=(8.0, 4.3), gridspec_kw={"wspace": 0.28})
        for ax, zoom in zip(axes, (False, True)):
            style(ax, zoom)
            handles = list(curves(ax, data[sp], aucs[sp]).values())
        axes[0].set_title("ROC", color=INK, fontsize=10, loc="left")
        axes[1].set_title("Zoom: FPR ≤ 0.2, TPR ≥ 0.8", color=INK, fontsize=10, loc="left")
        n = data[sp]["label"].value_counts()
        fig.suptitle(f"{name}: {n.get(1, 0):,} coding transcripts, {n.get(0, 0):,} lncRNAs",
                     color=INK, fontsize=11, x=0.06, ha="left", y=0.99)
        legend(fig, handles, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.02), title="ROC AUC",
               title_fontsize=8).get_title().set_color(INK_2)
        if name.endswith("†"):  # under the title: the legend fills the bottom band
            fig.text(0.06, 0.935, FOOTNOTE, color=INK_2, fontsize=8, ha="left", va="top")
        fig.subplots_adjust(left=0.07, right=0.98, top=0.86, bottom=0.27)
        for ext in ("png", "svg"):
            fig.savefig(args.out / "roc" / f"{sp}.{ext}", dpi=300)
        plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(10.5, 7.2), gridspec_kw={"wspace": 0.3, "hspace": 0.42})
    for ax, (sp, name) in zip(axes.flat, SPECIES.items()):
        style(ax, True)
        handles = curves(ax, data[sp], aucs[sp])
        ax.set_title(name, color=INK, fontsize=10, loc="left")
    last = axes.flat[-1]
    last.axis("off")
    fixed = [handles[t] for t in TOOLS]  # the shared legend names the tools only; AUCs differ per species
    for t, h in zip(TOOLS, fixed):
        h.set_label(TOOLS[t][0])
    legend(last, fixed, loc="center left", title="Zoom: FPR ≤ 0.2, TPR ≥ 0.8\nTransDecoder and GeneMark stay below\nTPR 0.8 here; ROC AUC per species: arena.tsv",
           title_fontsize=8).get_title().set_color(INK_2)
    fig.text(0.06, 0.012, FOOTNOTE, color=INK_2, fontsize=8, ha="left", va="bottom")
    fig.subplots_adjust(left=0.06, right=0.98, top=0.95, bottom=0.09)
    for ext in ("png", "svg"):
        fig.savefig(args.out / "roc" / f"all_species.{ext}", dpi=300)
    plt.close(fig)
    print(table[["species", "tool", "threshold", "TP", "TN", "FP", "FN", "F1", "ROC_AUC", "average_precision"]]
          .round(4).to_string(index=False))


if __name__ == "__main__":
    main()
