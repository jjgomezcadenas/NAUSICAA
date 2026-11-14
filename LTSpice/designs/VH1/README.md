# VH1: reset and two-capacitor CSA prototype

Prototype of the architecture proposed in `notes/physics/nausicaa_VH.tex`
(Section 7), built on the amplifier core of the present circuit so that the
two can be compared one to one.

- Core: M1 to M6 with their sizes and biases, the 1.25 pF pad capacitance
  and the 3 fF feedback value of `LTSpice/releases/csa_v3`; model
  `corrected.inc` copied from that release. The discharge transistor Mf is
  removed.
- Feedback: two capacitors Cf1, Cf2 of 3 fF between the input node `x` and
  nodes `a1`, `a2`. Ideal switches (`SW`, Ron 1 kohm, Roff 1e18 ohm) connect
  each of `a1`, `a2` either to `vref` (parked) or to `out` (feedback).
  Odd frames integrate on Cf1, even frames on Cf2, break before make with
  1 us gaps.
- Reset: a switch closes the loop from `out` to the gate through an ideal
  level shift `Vshift` of 0.885 V. In the present circuit the gate of M1
  sits at 1.3 mV and the output at 0.886 V, so a direct short of output to
  input would force a different operating point; with the short the solver
  finds out = in = 56 mV, a latched state with no gain. A real implementation
  needs a reset that returns the gate to its bias voltage, which the level
  shift emulates.
- `.options gmin=1e-18`: the default 1e-12 S across every off switch would
  discharge 3 fF with a 3 ms time constant.
- `Voff` between the gate (`in`) and the capacitor node (`x`) emulates an
  amplifier input offset for the cancellation test.

Decks are regenerated and run by
`scripts/electronics/design/vh1_reset_csa.py`; numerical results and plots
are in `results/ltspice/design/VH1/` (`summary.md`, `summary.json`, CSV).

| Deck | Purpose |
|---|---|
| `ac_cf3f.cir`, `ac_cf30f.cir` | Open-loop gain A = -V(out)/V(in), input impedance and transimpedance, closed loop with Cf |
| `ac_cf3f_nopad.cir`, `ac_cf3f_pad.cir` | Input capacitance of the amplifier alone and with the pad |
| `noise_cf3f.cir` | Output noise spectral density of the core in feedback, Mf removed |
| `swap_vref0886_voff0.cir` | Reset, 100 ions in frame 1, 200 in frame 2, none in frame 3 |
| `swap_vref0886_voff50mV.cir` | Same with a 50 mV input offset: cancellation test |
| `swap_vref1000_voff0.cir` | Same with Vref = 1.0 V: Vref independence and output range |
| `swap_packet100.cir` | Realistic 100-ion packet (`csa_v3` input current), frames aligned at 9 ms |
| `swap_range2000.cir` | 2000 ions per frame over six frames: linearity and range |

Switches are ideal: no charge injection, clock feedthrough or switch
capacitance is modelled. Reset kTC noise is not accessible to `.noise`; its
cancellation by the swap is tested deterministically through the offset run.
