"""Output helpers: tagged tables (CSV + Markdown) and a consistent figure style."""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .constants import OUT, FIG

# Colour-blind-safe qualitative palette (Okabe-Ito)
CB = ["#0072B2", "#D55E00", "#009E73", "#CC79A7",
      "#E69F00", "#56B4E9", "#F0E442", "#000000"]

plt.rcParams.update({
    "figure.figsize": (7.2, 4.6),
    "figure.dpi": 130,
    "font.size": 9.5,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "grid.linewidth": 0.6,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.prop_cycle": plt.cycler(color=CB),
    "legend.frameon": False,
    "savefig.bbox": "tight",
})

_MANIFEST = {}


def _cell(x, sig=4):
    if isinstance(x, (bool, np.bool_)):
        return "yes" if x else "no"
    if isinstance(x, (int, np.integer)):
        return str(int(x))
    if isinstance(x, (float, np.floating)):
        if not np.isfinite(x):
            return "∞" if x > 0 else "−∞"
        return f"{x:.{sig}g}"
    return str(x)


def to_markdown(df: pd.DataFrame, sig=4) -> str:
    """Minimal GitHub-flavoured Markdown table (avoids the tabulate dep)."""
    cols = list(df.columns)
    rows = [[_cell(v, sig) for v in rec] for rec in df.itertuples(index=False)]
    widths = [max(len(str(c)), *(len(r[i]) for r in rows)) if rows else len(str(c))
              for i, c in enumerate(cols)]
    out = ["| " + " | ".join(str(c).ljust(w) for c, w in zip(cols, widths)) + " |",
           "|" + "|".join("-" * (w + 2) for w in widths) + "|"]
    for r in rows:
        out.append("| " + " | ".join(v.ljust(w) for v, w in zip(r, widths)) + " |")
    return "\n".join(out)


def save_table(df: pd.DataFrame, name: str, caption: str = "", tags: dict = None,
               float_fmt="%.6g", md_rows=40):
    """Write a table to outputs/ as CSV plus a Markdown preview."""
    csv = OUT / f"{name}.csv"
    df.to_csv(csv, index=False, float_format=float_fmt)
    md = OUT / f"{name}.md"
    with open(md, "w", encoding="utf-8") as fh:
        if caption:
            fh.write(f"**{name}** — {caption}\n\n")
        if tags:
            fh.write("Traceability: " +
                     ", ".join(f"`{k}`={v}" for k, v in tags.items()) + "\n\n")
        head = df if len(df) <= md_rows else df.head(md_rows)
        fh.write(to_markdown(head))
        if len(df) > md_rows:
            fh.write(f"\n\n_({len(df)} rows total; full data in {name}.csv)_\n")
    _MANIFEST[name] = {"kind": "table", "rows": int(len(df)),
                       "caption": caption, "tags": tags or {}}
    return df


def save_fig(fig, name: str, caption: str = ""):
    """Save vector (SVG + PDF) and a PNG preview."""
    for ext in ("svg", "pdf", "png"):
        fig.savefig(FIG / f"{name}.{ext}")
    plt.close(fig)
    _MANIFEST[name] = {"kind": "figure", "caption": caption}
    return name


def save_json(obj, name: str):
    path = OUT / f"{name}.json"
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=2, default=_np_default)
    _MANIFEST[name] = {"kind": "json"}
    return path


def load_json(name: str):
    with open(OUT / f"{name}.json", encoding="utf-8") as fh:
        return json.load(fh)


def _np_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(f"not serialisable: {type(o)}")


def write_manifest():
    save_json(_MANIFEST, "_manifest")


def banner(text):
    line = "=" * max(60, len(text) + 4)
    print(f"\n{line}\n  {text}\n{line}")


def fmt_range(lo, hi, unit="", sig=3):
    def f(x):
        if not np.isfinite(x):
            return "∞"
        return f"{x:.{sig}g}"
    return f"{f(lo)}–{f(hi)} {unit}".strip()
