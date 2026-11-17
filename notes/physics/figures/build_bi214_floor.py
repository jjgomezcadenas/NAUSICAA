"""Figures, tables and macros of notes/physics/bi214_compton_floor.tex.

  python notes/physics/figures/build_bi214_floor.py [floor.h5]

Reads results/iitaca/II_L1000_D500_bi214_floor/floor.h5 (written by
scripts/iitaca/prepare_bi214_floor.jl) and writes, without repeating any
physical calculation:
  notes/physics/figures/bi214_floor_displays.pdf    event displays
  notes/physics/figures/bi214_floor_distance.pdf    acceptance against the distance of the second electron
  notes/physics/figures/bi214_floor_thirds.pdf      single-electron acceptance in thirds of the drift
  notes/physics/figures/bi214_floor_attenuation.pdf attenuation length and scattered-photon energy
  notes/physics/figures/bi214_floor_centre.pdf      distance from the Compton vertex to the second electron
  notes/physics/generated/bi214_floor_values.tex    numerical macros
"""
import sys, os, hashlib
import numpy as np, h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), "../../.."))
SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "results/iitaca/II_L1000_D500_bi214_floor/floor.h5")
FIG = os.path.join(ROOT, "notes/physics/figures")
GEN = os.path.join(ROOT, "notes/physics/generated")
BLUE, ORANGE, AQUA, YELLOW, GRAY = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#8a8a85"
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": "#e4e3df", "grid.linewidth": 0.6, "savefig.bbox": "tight"})

f = h5py.File(SRC, "r")
r = lambda k: f[k][()]
macros = {}

def pct(x, d=1): return f"{100 * x:.{d}f}"
def num(n): return f"{int(n):,}".replace(",", "{,}")
def put(name, value): macros[name] = value

# ---- cut flow, E2 cut, network -------------------------------------------
flow = dict(zip([s.decode() for s in r("flow/step")], r("flow/events")))
put("BiGenerated", num(flow["generated"])); put("BiSaved", num(flow["read"])); put("BiTrigger", num(flow["trigger"]))
put("BiROI", num(flow["roi"])); put("BiSelected", num(flow["track"]))
put("BiTrackFraction", pct(flow["track"] / flow["roi"]))
tc = {k: f[f"flow/trigger_counts/{k}"][()] for k in f["flow/trigger_counts"].keys()}
put("VetoGamma", num(tc["outside_frame_gamma_events"])); put("VetoGammaFrac", pct(tc["outside_frame_gamma_events"] / tc["input_events"], 0))
put("VetoEL", num(tc["main_track_near_el_events"])); put("VetoFrame", num(tc["frame_outside_detector_events"]))
put("TriggerFrac", pct(flow["trigger"] / flow["read"], 1)); put("RoiFrac", pct(flow["roi"] / flow["trigger"], 1))
eb, en = r("e2/bi214_passed"), r("e2/bi214_events"); sb, sn = r("e2/single_electrons_passed"), r("e2/single_electrons_events")
put("EtwoThreshold", f"{r('e2/threshold_keV'):.1f}")
put("BiEtwo", pct(eb / en)); put("BiEtwoErr", pct(np.sqrt(eb / en * (1 - eb / en) / en)))
put("SeEtwo", pct(sb / sn)); put("SeEtwoErr", pct(np.sqrt(sb / sn * (1 - sb / sn) / sn)))
put("SeEtwoEvents", num(sn))
nb, nn = r("network/bi214_passed"), r("network/bi214_events"); tb, tn = r("network/single_electrons_passed"), r("network/single_electrons_events")
put("BiNet", pct(nb / nn, 2)); put("BiNetErr", pct(np.sqrt(nb / nn * (1 - nb / nn) / nn), 2)); put("BiNetPassed", num(nb))
put("SeNet", pct(tb / tn, 2)); put("SeNetErr", pct(np.sqrt(tb / tn * (1 - tb / tn) / tn), 2)); put("SeNetPassed", num(tb)); put("SeNetEvents", num(tn))
put("BiNetRatio", f"{(nb / nn) / (tb / tn):.1f}")

# ---- truth -----------------------------------------------------------------
H = {k: r(f"truth/bi214/by_history/{k}") for k in ("events", "passed", "fraction", "error")}
put("PhotEvents", num(H["events"][1])); put("PhotPassed", num(H["passed"][1])); put("PhotNet", pct(H["fraction"][1], 2)); put("PhotNetErr", pct(H["error"][1], 2))
put("ComptEvents", num(H["events"][2])); put("ComptPassed", num(H["passed"][2])); put("ComptNet", pct(H["fraction"][2])); put("ComptNetErr", pct(H["error"][2]))
put("ComptShareSelected", pct(H["events"][2] / nn)); put("ComptSharePassed", pct(H["passed"][2] / nb, 0))
D = {k: r(f"truth/bi214/by_second_distance/{k}") for k in ("events", "passed", "fraction", "error", "e2_passed")}
for i, tag in enumerate(("A", "B", "C")):
    put(f"Dist{tag}Events", num(D["events"][i])); put(f"Dist{tag}Passed", num(D["passed"][i]))
    put(f"Dist{tag}Net", pct(D["fraction"][i])); put(f"Dist{tag}NetErr", pct(D["error"][i]))
    put(f"Dist{tag}Etwo", pct(D["e2_passed"][i] / D["events"][i], 0))
put("DistAShareSelected", pct(D["events"][0] / nn)); put("DistAContribution", pct(D["passed"][0] / nn, 2))
put("DistBContribution", pct(D["passed"][1] / nn, 2)); put("DistASharePassed", pct(D["passed"][0] / nb, 0))
put("DistARejection", f"{0.8 / D['fraction'][0]:.1f}")
X = {k: r(f"truth/bi214/by_xray/{k}") for k in ("events", "passed", "fraction", "error")}
put("XrayNoneNet", pct(X["fraction"][0], 2)); put("XrayNoneEvents", num(X["events"][0]))
put("XrayNearNet", pct(X["fraction"][2], 2)); put("XrayNearEvents", num(X["events"][2]))
put("ExpectedNet", pct(r("truth/expected_bi214_fraction"), 2)); put("ExpectedNetErr", pct(r("truth/expected_bi214_fraction_error"), 2))
excess = nb - tb / tn * nn
put("ExcessEvents", f"{excess:.0f}")
put("ComptExcess", f"{H['passed'][2] - tb / tn * H['events'][2]:.0f}")
cls, sd, ps = r("truth/class"), r("truth/second_distance_mm"), r("truth/passed")
lp, cs = r("truth/leading_process"), r("truth/second_centre_start_mm")
comp = cls == 3
put("ComptLeadCompton", num(np.sum(comp & (ps == 1) & (lp == 2))))
put("ComptCentreMedianPassed", f"{np.median(cs[comp & (ps == 1)]):.0f}")
put("ComptCentreMedianAll", f"{np.median(cs[comp]):.0f}")
put("ComptCentreNinetyPassed", f"{np.quantile(cs[comp & (ps == 1)], 0.9):.0f}")
put("ComptLtFour", pct(np.mean(sd[comp] < 4), 0))

fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9), gridspec_kw={"wspace": 0.3})
centres = np.array([2.0, 6.0, 10.0]); labels = ["< 4 mm", "4–8 mm", "> 8 mm"]
ax[0].errorbar(centres, 100 * D["fraction"], yerr=100 * D["error"], fmt="o", color=BLUE, ms=6, capsize=3, label="network, 80 % signal")
ax[0].errorbar(centres, 100 * D["e2_passed"] / D["events"], yerr=100 * np.sqrt(D["e2_passed"] / D["events"] * (1 - D["e2_passed"] / D["events"]) / D["events"]),
               fmt="s", color=ORANGE, ms=6, capsize=3, label="E$_2$ cut, 80 % signal")
ax[0].axhline(100 * tb / tn, color=BLUE, lw=1, ls="--"); ax[0].axhline(100 * sb / sn, color=ORANGE, lw=1, ls="--")
ax[0].text(4.4, 100 * tb / tn + 1.0, "single e$^-$, network", ha="left", fontsize=6.5, color="#444")
ax[0].text(0.2, 100 * sb / sn + 1.0, "single e$^-$, E$_2$ cut", ha="left", fontsize=6.5, color="#444")
for c, n, k in zip(centres, D["events"], D["passed"]):
    ax[0].annotate(f"{k}/{n}", (c, 100 * k / n), textcoords="offset points", xytext=(8, -3), fontsize=7)
ax[0].set_xticks(centres, labels); ax[0].set_xlim(0, 12); ax[0].set_ylim(0, 105)
ax[0].set_xlabel("second electron: distance to main track"); ax[0].set_ylabel("acceptance (%)")
ax[0].set_title("Compton then photoelectric, selected $^{214}$Bi", fontsize=8.5)
ax[0].legend(fontsize=7, frameon=False, loc="upper right")
bins = np.arange(0, 30.5, 1.0)
ax[1].hist(np.clip(sd[comp], 0, 30), bins, color=GRAY, alpha=0.5, label=f"all ({comp.sum()})")
ax[1].hist(np.clip(sd[comp & (ps == 1)], 0, 30), bins, color=BLUE, label=f"passing the network ({np.sum(comp & (ps == 1))})")
ax[1].axvline(4, color="#444", lw=0.8, ls=":"); ax[1].axvline(8, color="#444", lw=0.8, ls=":")
ax[1].set_xlabel("distance to main track (mm)"); ax[1].set_ylabel("events / mm")
ax[1].set_title("same events; last bin includes > 30 mm", fontsize=8.5); ax[1].legend(fontsize=7, frameon=False)
fig.savefig(os.path.join(FIG, "bi214_floor_distance.pdf")); plt.close(fig)

fig, ax = plt.subplots(figsize=(3.6, 2.7))
b2 = np.arange(0, 61, 2.0)
ax.hist(np.clip(cs[comp], 0, 60), b2, color=GRAY, alpha=0.5, label="all")
ax.hist(np.clip(cs[comp & (ps == 1)], 0, 60), b2, color=BLUE, label="passing")
ax.set_xlabel("start of leading electron to energy centre\nof second electron (mm)"); ax.set_ylabel("events / 2 mm")
ax.set_title(f"selected $^{{214}}$Bi, Compton class ({comp.sum()}); last bin > 60 mm", fontsize=8)
ax.legend(fontsize=7, frameon=False)
fig.savefig(os.path.join(FIG, "bi214_floor_centre.pdf")); plt.close(fig)

# ---- drift thirds ---------------------------------------------------------------
T = {}
for tag in ("II_L1000_D500", "xehe_itaca"):
    g = f[f"thirds/{tag}"]; T[tag] = {k: g[k][()] for k in g.keys()}
fig, ax = plt.subplots(figsize=(3.8, 2.7))
for tag, col, mk, lab in (("II_L1000_D500", BLUE, "o", "II_L1000_D500, pure Xe"), ("xehe_itaca", ORANGE, "s", "ITACA, Xe–15 % He")):
    t = T[tag]
    ax.errorbar(t["sigma_median_mm"], 100 * t["fraction"], yerr=100 * t["error"],
                xerr=[t["sigma_median_mm"] - t["sigma_lower_mm"], t["sigma_upper_mm"] - t["sigma_median_mm"]],
                fmt=mk, color=col, ms=6, capsize=2, lw=1, label=lab)
    ax.axhline(100 * t["overall_fraction"], color=col, lw=0.8, ls="--")
ax.set_xlabel("thermal ion spread $\\sigma$ in the third (mm)"); ax.set_ylabel("single e$^-$ accepted (%)")
ax.set_ylim(0, 2.5); ax.set_xlim(0, 1.45)
ax.set_title("network at 80 % signal; thirds of the drift length;\nbars: spread at the edges, dashed: all test events", fontsize=8)
ax.legend(fontsize=7, frameon=False, loc="lower right")
fig.savefig(os.path.join(FIG, "bi214_floor_thirds.pdf")); plt.close(fig)
for tag, p in (("II_L1000_D500", "Pure"), ("xehe_itaca", "Hybrid")):
    t = T[tag]
    for i, w in enumerate(("One", "Two", "Three")):
        put(f"Third{p}{w}", pct(t["fraction"][i], 2)); put(f"Third{p}{w}Err", pct(t["error"][i], 2))
        put(f"Third{p}{w}Sigma", f"{t['sigma_median_mm'][i]:.2f}"); put(f"Third{p}{w}Events", num(t["events"][i]))
        put(f"Third{p}{w}Lo", f"{t['edges_cm'][i]:.0f}"); put(f"Third{p}{w}Hi", f"{t['edges_cm'][i + 1]:.0f}")
        put(f"Third{p}{w}Sig", pct(t["signal_efficiency"][i], 1))
    put(f"Third{p}All", pct(t["overall_fraction"], 2)); put(f"Third{p}AllErr", pct(t["overall_error"], 2))
    put(f"Third{p}FullSigma", f"{t['full_drift_sigma_mm']:.2f}")

# ---- attenuation --------------------------------------------------------------------
A = lambda k: r(f"attenuation/{k}")
put("GasDensity", f"{A('density_kg_m3'):.1f}")
put("BackscatterKeV", f"{A('backscatter_keV'):.0f}")
put("NinetyKeV", f"{A('compton_scattered_keV')[np.argmin(np.abs(A('compton_angle_rad') - np.pi / 2))]:.0f}"); put("BackscatterLambda", f"{A('backscatter_length_mm') / 10:.0f}")
for e in (300, 400, 500):
    put(f"Lambda{['Three', 'Four', 'Five'][e // 100 - 3]}", f"{A(f'at_{e}keV_length_mm') / 10:.0f}")
    put(f"Pten{['Three', 'Four', 'Five'][e // 100 - 3]}", pct(A(f"at_{e}keV_p10mm"), 1))
put("LineLambda", f"{A('line_length_mm') / 10:.0f}")
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.6))
ax[0].plot(np.degrees(A("compton_angle_rad")), A("compton_scattered_keV"), color=BLUE, lw=2)
ax[0].set_xlabel("Compton scattering angle (deg)"); ax[0].set_ylabel("scattered photon energy (keV)")
ax[0].set_title("2447.86 keV line", fontsize=8.5); ax[0].set_xlim(0, 180); ax[0].set_ylim(0, 2600)
ax[1].plot(A("energy_keV"), A("length_mm") / 10, color=BLUE, lw=2)
ax[1].set_xlabel("photon energy (keV)"); ax[1].set_ylabel("attenuation length (cm)")
ax[1].set_title(f"xenon, {A('density_kg_m3'):.1f} kg/m$^3$ (16 bar, 290 K); NIST SRD 126", fontsize=8.5)
ax[1].set_xlim(200, 2500); ax[1].set_ylim(0, 280)
fig.savefig(os.path.join(FIG, "bi214_floor_attenuation.pdf")); plt.close(fig)

# ---- voxelised Monte Carlo ------------------------------------------------------------------
V = {k: f[f"voxelised/{k}"][()] for k in f["voxelised"].keys()}
vs, rs = V["selected"] == 1, V["readout_selected"] == 1
put("VoxEvents", num(len(vs))); put("VoxSelected", num(vs.sum())); put("VoxFraction", pct(vs.mean()))
put("VoxBoth", num((vs & rs).sum())); put("VoxOnly", num((vs & ~rs).sum())); put("ReadoutOnly", num((~vs & rs).sum()))
put("VoxReadoutSelected", num(rs.sum())); put("VoxBoxDiff", f"{V['max_box_energy_difference_keV']:.3f}")
ev_v = dict(zip(V["event_number"], range(len(vs))))
idx = np.array([ev_v[e] for e in r("truth/event_number")])
cv = vs[idx]
put("VoxComptKept", num(cv[comp].sum())); put("VoxComptKeptFrac", pct(cv[comp].mean(), 0))
put("VoxComptPassKept", num(cv[comp & (ps == 1)].sum())); put("VoxComptPass", num((comp & (ps == 1)).sum()))
for i, (lo, hi) in enumerate(((0, 4), (4, 8), (8, 1e9))):
    m = comp & (sd >= lo) & (sd < hi)
    put(f"VoxDist{'ABC'[i]}Kept", pct(cv[m].mean(), 0)); put(f"VoxDist{'ABC'[i]}Events", num(m.sum()))
put("VoxPhotKept", pct(cv[cls != 3].mean(), 0))
put("VoxPassKept", num(cv[ps == 1].sum())); put("VoxNonComptPassKept", num(cv[(ps == 1) & ~comp].sum()))
# reasons: deposits and segments
put("VoxRejDeposit", num((rs & ~vs & (V["deposits"] > 0)).sum()))
put("ReadoutRejDeposit", num((vs & ~rs & (V["readout_deposits"] > 0)).sum()))

# ---- event displays ------------------------------------------------------------------------
order = [("compton_lt4mm", "Compton; second electron {dmin} mm from the main track, photon absorbed {dist} mm from the vertex; passes the network"),
         ("compton_4to8mm", "Compton; second electron {dmin} mm from the main track, photon absorbed {dist} mm from the vertex; passes the network"),
         ("photoelectric", "photoelectric absorption only; passes the network"),
         ("compton_rejected", "Compton; photon absorbed {dist} mm from the vertex; rejected by the track selection")]
lab_names = {1: "main track", 2: "second electron", 3: "fluorescence photon", 4: "other (bremsstrahlung, relaxation)"}
lab_cols = {1: BLUE, 2: ORANGE, 3: AQUA, 4: GRAY}
fig = plt.figure(figsize=(7.4, 9.4))
gs = fig.add_gridspec(4, 4, left=0.07, right=0.99, top=0.965, bottom=0.07, hspace=0.42, wspace=0.42)
handles = {}
for row, (tag, title) in enumerate(order):
    g = f[f"displays/{tag}"]; d = {k: g[k][()] for k in g.keys()}
    mx, my, mz, me, ml = d["mc_x_mm"], d["mc_y_mm"], d["mc_z_mm"], d["mc_energy_keV"], d["mc_label"]
    off = d["vertex_offset_mm"]
    rx, ry, rz, re = d["reco_x_mm"] + off[0], d["reco_y_mm"] + off[1], d["reco_z_mm"] + off[2], d["reco_energy_keV"]
    keep = re > 5.0
    c = [(mx.max() + mx.min()) / 2, (my.max() + my.min()) / 2, (mz.max() + mz.min()) / 2]
    half = max(np.ptp(mx), np.ptp(my), np.ptp(mz)) / 2 + 8
    dist = d.get("second_start_distance_mm", None)
    if "second_keV" in d:
        put({"compton_lt4mm": "DispA", "compton_4to8mm": "DispB", "compton_rejected": "DispR"}[tag] + "SecondKeV", f"{d['second_keV']:.0f}")
        put({"compton_lt4mm": "DispA", "compton_4to8mm": "DispB", "compton_rejected": "DispR"}[tag] + "PhotonPath", f"{dist:.0f}")
    ev = int(d["event_number"])
    i_t = np.where(r("truth/event_number") == ev)[0]
    dmin = sd[i_t[0]] if len(i_t) else np.nan
    axs = [fig.add_subplot(gs[row, k]) for k in range(4)]
    for colm, (a, b, vert) in enumerate(((mx, my, "y"), (mx, mz, "z"))):
        ax = axs[colm]
        for Lb in (1, 4, 3, 2):
            m = ml == Lb
            if m.any(): handles[Lb] = ax.scatter(a[m], b[m], s=np.clip(me[m], 1.0, 40) * 0.7, color=lab_cols[Lb], lw=0, zorder=2 + (Lb != 1))
        if "second_start_mm" in d:
            s = d["second_start_mm"]; ax.plot(s[0], s[1] if vert == "y" else s[2], "x", color="black", ms=7, mew=1.3, zorder=5)
        ax.plot(0, 0, "o", mfc="none", mec="black", ms=7, mew=1, zorder=5)
    for colm, (a, b, vert) in enumerate(((rx, ry, "y"), (rx, rz, "z"))):
        ax = axs[2 + colm]
        ax.scatter(a[keep], b[keep], c=re[keep], s=10, marker="s", cmap="Blues", vmin=-0.2 * np.quantile(re[keep], 0.95), vmax=np.quantile(re[keep], 0.95), lw=0)
        ax.plot(0, 0, "o", mfc="none", mec="black", ms=7, mew=1, zorder=5)
    for colm, ax in enumerate(axs):
        vert = "y" if colm % 2 == 0 else "z"
        ax.set_xlim(c[0] - half, c[0] + half); ax.set_ylim((c[1] if vert == "y" else c[2]) - half, (c[1] if vert == "y" else c[2]) + half)
        ax.set_aspect("equal"); ax.tick_params(labelsize=6.5); ax.grid(False)
        ax.set_xlabel("x (mm)", fontsize=7, labelpad=1); ax.set_ylabel(f"{vert} (mm)", fontsize=7, labelpad=0)
        ax.set_title(("Geant4 hits" if colm < 2 else "reconstructed hits > 5 keV") + f", x{vert}", fontsize=6.8, pad=2)
    tt = f"event {ev}: " + title.format(dmin=f"{dmin:.1f}", dist=f"{dist:.0f}" if dist is not None else "")
    bb = axs[0].get_position()
    fig.text(0.07, bb.y1 + 0.022, tt, fontsize=7.5, ha="left", weight="bold")
from matplotlib.lines import Line2D
hl = [Line2D([], [], marker="o", ls="", color=lab_cols[k], ms=5) for k in (1, 2, 3, 4)] + \
     [Line2D([], [], marker="o", ls="", mfc="none", mec="black", ms=6), Line2D([], [], marker="x", ls="", color="black", ms=6)]
fig.legend(hl, [lab_names[k] for k in (1, 2, 3, 4)] + ["start of the main track", "start of the second electron"],
           loc="lower center", ncol=3, fontsize=7, frameon=False, bbox_to_anchor=(0.5, 0.0))
fig.savefig(os.path.join(FIG, "bi214_floor_displays.pdf")); plt.close(fig)

# ---- 214Bi in drift thirds, pitch estimate -------------------------------------------------
B = {k: f[f"bi214_thirds/{k}"][()] for k in f["bi214_thirds"].keys()}
for i, w in enumerate(("One", "Two", "Three")):
    put(f"BiThird{w}Lo", f"{B['edges_cm'][i]:.0f}"); put(f"BiThird{w}Hi", f"{B['edges_cm'][i + 1]:.0f}")
    put(f"BiThird{w}Sel", pct(B["selected"][i] / B["roi"][i])); put(f"BiThird{w}Vox", pct(B["voxelised_selected"][i] / B["voxelised_roi"][i]))
    a = B["passed"][i] / B["selected"][i]; put(f"BiThird{w}Net", pct(a, 2)); put(f"BiThird{w}NetErr", pct(np.sqrt(a * (1 - a) / B["selected"][i]), 2))
    put(f"BiThird{w}Compt", num(B["compton_selected"][i])); put(f"BiThird{w}ComptShare", pct(B["compton_selected"][i] / B["selected"][i]))
    c = B["compton_passed"][i] / B["compton_selected"][i]; put(f"BiThird{w}ComptNet", pct(c, 0)); put(f"BiThird{w}ComptNetErr", pct(np.sqrt(c * (1 - c) / B["compton_selected"][i]), 0))
    put(f"BiThird{w}PerRoi", f"{1000 * B['passed'][i] / B['roi'][i]:.2f}")
put("BiThirdRoi", num(B["roi"][0]))
for i, w in enumerate(("One", "Two", "Three")):
    put(f"BiThird{w}ReadoutOnly", pct(B["readout_only"][i] / B["roi"][i])); put(f"BiThird{w}VoxOnly", pct(B["voxelised_only"][i] / B["roi"][i]))
    put(f"BiThird{w}SelN", num(B["selected"][i]))
P = {k: f[f"pitch_estimate/{k}"][()] for k in f["pitch_estimate"].keys()}
put("PitchEdgeOne", f"{P['edges_mm'][1]:.1f}"); put("PitchEdgeTwo", f"{P['edges_mm'][2]:.1f}")
for i, w in enumerate("ABC"): put(f"PitchDist{w}Events", num(P["events"][i]))
put("PitchComptPassed", f"{P['compton_passed']:.0f}"); put("PitchComptNet", pct(P["compton_fraction"], 0))
put("PitchTotalNet", pct(P["total_fraction"], 2)); put("PitchOthers", num(P["others_passed"]))
put("PitchIrreducible", f"{P['events'][0] * P['acceptance_3mm'][0]:.0f}")

# ---- macros ------------------------------------------------------------------------------------
sha = hashlib.sha256(open(SRC, "rb").read()).hexdigest()
with open(os.path.join(GEN, "bi214_floor_values.tex"), "w") as out:
    out.write(f"% Generated by notes/physics/figures/build_bi214_floor.py from {os.path.relpath(SRC, ROOT)}\n% sha256 {sha}\n")
    for k, v in macros.items():
        out.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")
print("wrote", len(macros), "macros")
