"""
Shared look for all cfCNV figures. EDIT THIS FILE to change colours, fonts and sizes;
every figure script imports it, so one change restyles all figures.
"""
import matplotlib.pyplot as plt

# ---- caller colours (one per caller; change freely) ---------------------------------
TOOL_COLORS = {
    "QDNAseq":   "#1b9e77",
    "ichorCNA":  "#d95f02",
    "BIC-seq2":  "#7570b3",
    "CNVkit":    "#e7298a",
    "Spectre":   "#66a61e",
    "CNV-MEANN": "#e6ab02",
    "dudeML":    "#a6761d",
    "LYCEUM":    "#666666",
    "ACE":       "#999999",
    "Fragle":    "#1f78b4",
}
ML_TOOLS = {"CNV-MEANN", "dudeML", "LYCEUM"}          # drawn with dashed lines
TOOL_ORDER = ["QDNAseq", "ichorCNA", "BIC-seq2", "CNVkit", "Spectre", "CNV-MEANN", "dudeML", "LYCEUM"]
GAIN_COLOR = "#c0392b"            # copy-number gain
LOSS_COLOR = "#2c6fbb"            # copy-number loss
BACKGROUND_COLORS = {"IH01": "#444444", "IH02": "#888888", "IH03": "#bbbbbb"}
HEATMAP_CMAP = "viridis"           # sequential colour map for heatmaps (e.g. "magma", "Blues")
SCORE_CMAP = ["#b2182b", "#f7f7f7", "#2166ac"]   # worse -> better, summary scorecard

# ---- page geometry and fonts ----------------------------------------------------------
MM = 1 / 25.4                       # millimetres -> inches
FULL_WIDTH = 180 * MM               # double-column figure width
FONT_FAMILY = "Arial"               # falls back to DejaVu Sans if Arial is not installed
BASE_SIZE = 7
DPI = 600                           # raster output resolution

def apply():
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": [FONT_FAMILY, "Helvetica", "DejaVu Sans"],
        "font.size": BASE_SIZE, "axes.titlesize": BASE_SIZE + 1, "axes.labelsize": BASE_SIZE,
        "xtick.labelsize": BASE_SIZE - 0.5, "ytick.labelsize": BASE_SIZE - 0.5, "legend.fontsize": BASE_SIZE - 1,
        "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6, "lines.linewidth": 1.2,
        "legend.frameon": False, "pdf.fonttype": 42, "svg.fonttype": "none",
        "savefig.bbox": "tight", "savefig.dpi": DPI,
    })

def color(tool):
    return TOOL_COLORS.get(tool, "#333333")

def panel(ax, letter, x=-0.14, y=1.06):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=BASE_SIZE + 3, fontweight="bold", va="bottom", ha="left")

def save(fig, name, outdir="output"):
    import os
    os.makedirs(outdir, exist_ok=True)
    for ext in ("pdf", "png", "svg"):
        fig.savefig(f"{outdir}/{name}.{ext}")
    print("wrote", f"{outdir}/{name}.[pdf|png|svg]")
