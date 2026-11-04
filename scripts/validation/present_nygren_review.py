"""Generate TeX numerical tables from the Julia analytical-budget HDF5."""
from pathlib import Path
import h5py
ROOT = Path(__file__).resolve().parents[2]
out = ROOT / "latex_generated/nygren_review"
out.mkdir(parents=True, exist_ok=True)
with h5py.File(ROOT / "results/nygren_review/budgets.h5") as h:
    # Julia matrices have reversed dimensions when read through h5py.
    res = h["resolution_fwhm_fraction"][:].T
    limits = h["maximum_independent_channels"][:].T
    times = h["arrival_duration_s"][:].T
    def table(headers, rows):
        return (r"\begin{tabular}{" + "l" + "r"*(len(headers)-1) + r"}\toprule" + "\n"
                + " & ".join(headers) + r"\\\midrule" + "\n"
                + "\n".join(" & ".join(map(str, row)) + r"\\" for row in rows)
                + "\n" + r"\bottomrule\end{tabular}" + "\n")
    (out / "energy_table.tex").write_text(table(
        [r"Independent channels $n$", *[f"{e:g} electrons ENC" for e in h["enc_e"][:]]],
        [[n, *[f"{100*x:.3f}\\%" for x in row]] for n, row in zip(h["independent_channels"][:], res)]))
    (out / "limits_table.tex").write_text(table(
        ["Target FWHM", "Total noise budget (electrons)", r"Max. $n$, ENC 50", r"Max. $n$, ENC 120"],
        [[f"{100*r:.1f}\\%", f"{b:.1f}", f"{row[1]:.1f}", f"{row[2]:.1f}"]
         for r,b,row in zip(h["target_fwhm_fraction"][:],h["target_noise_budget_e"][:],limits)]))
    (out / "timing_table.tex").write_text(table(
        [r"Extent $\Delta z$ (mm)", r"$v=0.1$ m/s: time (s)", r"$v=0.2$ m/s: time (s)"],
        [[f"{z*1000:g}", *[f"{x:g}" for x in row]] for z,row in zip(h["z_extent_m"][:],times)]))
    (out / "values.tex").write_text("\n".join([
        rf"\newcommand{{\PrimaryPairs}}{{{h['primary_pairs'][()]:.0f}}}",
        rf"\newcommand{{\IntrinsicFWHM}}{{{100*h['intrinsic_fwhm_fraction'][()]:.3f}\%}}",
        rf"\newcommand{{\ThousandENCResolution}}{{{100*h['electronic_1000e_fwhm_fraction'][()]:.2f}\%}}",
        rf"\newcommand{{\SpeedLow}}{{{h['mobility_speed_m_s'][0]:.3f}}}",
        rf"\newcommand{{\SpeedHigh}}{{{h['mobility_speed_m_s'][1]:.3f}}}",
        rf"\newcommand{{\LeakageShotOneSecond}}{{{h['shot_noise_e'][2]:.1f}}}",
    ]) + "\n")
print(out)
