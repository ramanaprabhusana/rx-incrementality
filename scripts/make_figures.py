"""Render README figures from results/estimates.json, in light and dark themes.

Palette and mark specs follow a validated categorical palette: series blue and
orange pass colour-vision-deficiency separation in both themes (worst adjacent
Delta E 24.7 light, 26.8 dark). Text uses ink tokens, never series colours.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Rectangle  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "figures"
DPI = 200

THEMES = {
    "light": dict(surface="#fcfcfb", ink="#0b0b0b", ink2="#52514e", muted="#898781",
                  grid="#e1e0d9", base="#c3c2b7", s1="#2a78d6", s2="#eb6834", dim="#b4b2aa"),
    "dark": dict(surface="#1a1a19", ink="#ffffff", ink2="#c3c2b7", muted="#898781",
                 grid="#2c2c2a", base="#383835", s1="#3987e5", s2="#d95926", dim="#5c5b57"),
}


def pct(b: float) -> float:
    return 100.0 * (math.exp(b) - 1.0)


def style(ax, t):
    ax.set_facecolor(t["surface"])
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(t["base"])
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(colors=t["muted"], labelcolor=t["ink2"], length=0, labelsize=9)
    ax.grid(axis="y", color=t["grid"], linewidth=0.6, linestyle="-")
    ax.set_axisbelow(True)


def header(fig, t, title, subtitle):
    fig.text(0.06, 0.95, title, color=t["ink"], fontsize=13, fontweight="bold", va="top")
    fig.text(0.06, 0.885, subtitle, color=t["ink2"], fontsize=9.5, va="top")


def _points(ax, t, rel, coefs, ses, colour, scale, dx=0.0):
    xs = [r + dx for r in rel]
    for x, b, se in zip(xs, coefs, ses, strict=True):
        ax.plot([x, x], [scale(b - 1.96 * se), scale(b + 1.96 * se)], color=colour,
                linewidth=1.6, solid_capstyle="round", zorder=2)
    ax.plot(xs, [scale(b) for b in coefs], "o", color=colour, markersize=6, zorder=3,
            markeredgecolor=t["surface"], markeredgewidth=1.4)


def _event_axes(ax, t, ylim, ticks, labels, ylabel):
    style(ax, t)
    ax.axhline(0, color=t["base"], linewidth=0.8, zorder=1)
    ax.axvline(-0.5, color=t["base"], linewidth=0.8, zorder=1)
    ax.plot([-1], [0], "o", color=t["muted"], markersize=4.5, zorder=3)
    ax.set_xticks(range(-5, 5))
    ax.set_xticklabels([f"{r:+d}" if r else "0" for r in range(-5, 5)])
    ax.set_xlim(-5.6, 4.6)
    ax.set_ylim(*ylim)
    # Ticks at clean values. A rounding formatter on 2.5% steps once printed 7.5% as "+8%".
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels)
    ax.set_xlabel("Years relative to first payment", color=t["ink2"], fontsize=8.5)
    ax.set_ylabel(ylabel, color=t["ink2"], fontsize=8.5)


def event_study(events, theme):
    t = THEMES[theme]
    fig = plt.figure(figsize=(8, 4.6), dpi=DPI, facecolor=t["surface"])
    left = fig.add_axes([0.08, 0.13, 0.40, 0.58])
    right = fig.add_axes([0.58, 0.13, 0.40, 0.58])

    intensive = events["physician"]["cohort (Sun-Abraham)"]
    extensive = events["physician_extensive"]["cohort (Sun-Abraham)"]

    _event_axes(left, t, (-11, 12), [-10, -5, 0, 5, 10], ["-10%", "-5%", "0", "+5%", "+10%"],
                "Change in claims")
    _points(left, t, intensive["rel_periods"], intensive["coefs"], intensive["ses"], t["s1"], pct)
    left.text(-0.42, 11.3, "first payment", color=t["muted"], fontsize=7.5, va="top")
    p_int = intensive["pretrend_p"]
    left.set_title(f"How much established prescribers use the drug\npre-trend test passes (p = {p_int:.2f})",
                   color=t["ink2"], fontsize=8.5, loc="left", pad=6)

    to_pts = lambda b: 100.0 * b  # noqa: E731
    _event_axes(right, t, (-7, 10.5), [-5, 0, 5, 10], ["-5", "0", "+5", "+10"], "Change in probability (points)")
    _points(right, t, extensive["rel_periods"], extensive["coefs"], extensive["ses"], t["s1"], to_pts, dx=-0.1)
    dt = extensive["detrended"]
    post = [i for i, r in enumerate(extensive["rel_periods"]) if r >= 0]
    _points(right, t, [extensive["rel_periods"][i] for i in post], [dt["coefs"][i] for i in post],
            [dt["ses"][i] for i in post], t["s2"], to_pts, dx=0.1)
    right.set_title("Whether physicians prescribe the drug at all\npre-trend test fails (p < 0.001)",
                    color=t["ink2"], fontsize=8.5, loc="left", pad=6)
    handles = [
        Line2D([], [], marker="o", linestyle="", color=t["s1"], markersize=6,
               markeredgecolor=t["surface"], markeredgewidth=1.4, label="Estimate"),
        Line2D([], [], marker="o", linestyle="", color=t["s2"], markersize=6,
               markeredgecolor=t["surface"], markeredgewidth=1.4, label="Net of the pre-trend"),
    ]
    leg = right.legend(handles=handles, loc="lower right", frameon=False, fontsize=8,
                       handletextpad=0.3, borderaxespad=0.2)
    for txt in leg.get_texts():
        txt.set_color(t["ink2"])

    header(fig, t, "Prescribing rises after the first payment, but adoption was already rising",
           "71,869 physicians, diabetes drugs, Medicare Part D 2019 to 2024. "
           "Cohort event-study estimates with 95% intervals.")
    path = OUT / f"event-study-{theme}.png"
    fig.savefig(path, dpi=DPI, facecolor=t["surface"])
    plt.close(fig)
    return path


def rounded_hbar(ax, y, value, height, colour, radius_px, surface):
    """Horizontal bar: square at the zero baseline, rounded at the data end."""
    to_px = ax.transData
    px_x = abs(to_px.transform((1, 0))[0] - to_px.transform((0, 0))[0])
    px_y = abs(to_px.transform((0, 1))[1] - to_px.transform((0, 0))[1])
    rx, ry = radius_px / px_x, radius_px / px_y
    x0, x1 = 0.0, value
    if abs(x1) <= 2 * rx:
        ax.add_patch(Rectangle((min(x0, x1), y - height / 2), abs(x1), height, color=colour, lw=0))
        return
    ax.add_patch(FancyBboxPatch((min(x0, x1), y - height / 2), abs(x1), height,
                                boxstyle=f"round,pad=0,rounding_size={rx}",
                                mutation_aspect=ry / rx, color=colour, lw=0))
    sq_start = x0 if x1 > 0 else x1 + rx
    ax.add_patch(Rectangle((sq_start, y - height / 2), abs(x1) - rx, height, color=colour, lw=0))


def specification(results, theme):
    t = THEMES[theme]
    want = [("2-way: physician-year + drug-year", "Physician-year +\ndrug-year"),
            ("published: physician-drug + drug-year", "Physician-drug +\ndrug-year (published)"),
            ("3-way: all three", "All three\n(this repository)")]
    by = {r["spec"]: r for r in results if "coef" in r}
    fig = plt.figure(figsize=(8, 3.9), dpi=DPI, facecolor=t["surface"])
    axes = [fig.add_axes([0.25, 0.17, 0.33, 0.56]), fig.add_axes([0.64, 0.17, 0.33, 0.56])]
    panels = [("Effect of this year's payment", lambda s: by[s], "coef", "se", 30, [0, 10, 20, 30]),
              ("Apparent effect of next year's payment\n(should be zero)",
               lambda s: by[f"{s} + next year's payment"], "lead", "lead_se", 18, [0, 5, 10, 15])]
    ys = [2, 1, 0]
    for ax, (ptitle, getter, coef_key, se_key, xmax, ticks) in zip(axes, panels, strict=True):
        style(ax, t)
        ax.grid(axis="y", visible=False)
        ax.grid(axis="x", color=t["grid"], linewidth=0.6)
        ax.set_xlim(0, xmax)
        ax.set_ylim(-0.6, 2.6)
        fig.canvas.draw()
        for y, (spec, _) in zip(ys, want, strict=True):
            r = getter(spec)
            b = r["coef"] if coef_key == "coef" else r["extra"]["pay_any_lead"]
            s = r["se"] if se_key == "se" else r["extra"]["pay_any_lead_se"]
            colour = t["s1"] if spec.startswith("3-way") else t["dim"]
            rounded_hbar(ax, y, pct(b), 0.42, colour, 6, t["surface"])
            ax.plot([pct(b - 1.96 * s), pct(b + 1.96 * s)], [y, y], color=t["ink2"],
                    linewidth=1.0, zorder=3)
            ax.text(pct(b + 1.96 * s) + xmax * 0.02, y, f"{pct(b):+.1f}%", va="center",
                    color=t["ink"] if spec.startswith("3-way") else t["ink2"], fontsize=9,
                    fontweight="bold" if spec.startswith("3-way") else "normal")
        ax.set_title(ptitle, color=t["ink2"], fontsize=9, loc="left", pad=6)
        ax.set_xticks(ticks)
        ax.set_xticklabels([f"{v}%" for v in ticks])
        ax.set_yticks(ys)
        ax.set_yticklabels([lab for _, lab in want] if ax is axes[0] else [])
    header(fig, t, "Most of the naive effect is targeting",
           "Physicians, diabetes drugs, 2019 to 2024. A payment next year cannot cause prescribing this year.")
    path = OUT / f"specification-{theme}.png"
    fig.savefig(path, dpi=DPI, facecolor=t["surface"])
    plt.close(fig)
    return path


def main() -> None:
    plt.rcParams["font.family"] = ["Helvetica Neue", "Arial", "DejaVu Sans"]
    data = json.loads((ROOT / "results" / "estimates.json").read_text())
    OUT.mkdir(parents=True, exist_ok=True)
    for theme in THEMES:
        print(event_study(data["event_studies"], theme).relative_to(ROOT))
        print(specification(data["estimates"], theme).relative_to(ROOT))


if __name__ == "__main__":
    main()
