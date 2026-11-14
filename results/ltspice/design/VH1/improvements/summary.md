# VH1 improvements

Every variant re-tuned with the reset loop closed: Vnload (M3) for I(M1) = 210.3 uA, Vnmos (M2) for V(out) = 0.886 V.
ENC for 20 ms frames, frame difference of 5 ms averages (as the present calibration) and of 20 ms averages.

| variant | Vnmos (V) | Vnload (V) | A0 | C_amp (pF) | charge fraction | uV/ion | ENC 5 ms avg | ENC 20 ms avg |
|---|---|---|---|---|---|---|---|---|
| base | 1.600 | 0.538 | 7605 | 4.34 | 0.803 | 42.9 | 40.2 | 27.5 |
| L2_1u | 1.613 | 0.538 | 6816 | 3.92 | 0.798 | 42.6 | 34.5 | 23.6 |
| L2_2u | 1.650 | 0.538 | 5720 | 3.34 | 0.789 | 42.1 | 34.0 | 23.3 |
| L2_4u | 1.726 | 0.538 | 4482 | 2.69 | 0.773 | 41.3 | 37.6 | 25.7 |
| L3_4u | 1.600 | 0.582 | 3669 | 2.45 | 0.749 | 40.0 | 38.4 | 26.2 |
| L3_8u | 1.600 | 0.661 | 1363 | 1.34 | 0.612 | 32.7 | 43.1 | 29.4 |
| L3_16u | 1.600 | 0.819 | 755 | 1.05 | 0.496 | 26.5 | 53.2 | 36.3 |
| W4_5.76u | 1.600 | 0.538 | 11016 | 6.05 | 0.819 | 43.7 | 39.7 | 27.1 |
| W4_11.5u | 1.600 | 0.538 | 14281 | 7.70 | 0.827 | 44.2 | 39.6 | 27.1 |
| W4_23u | 1.600 | 0.538 | 17213 | 9.15 | 0.832 | 44.5 | 39.6 | 27.0 |
| W4_46u | 1.600 | 0.538 | 19745 | 10.31 | 0.837 | 44.7 | 39.6 | 27.1 |
| W1_110u | 1.654 | 0.538 | 6353 | 1.91 | 0.858 | 45.8 | 39.4 | 27.1 |
| W1_55u | 1.729 | 0.538 | 4977 | 0.81 | 0.879 | 46.9 | 47.4 | 33.0 |
| combined | 1.705 | 0.582 | 2379 | 0.90 | 0.768 | 41.0 | 30.3 | 21.0 |
| combined_W4_23u | 1.705 | 0.580 | 9269 | 2.63 | 0.878 | 46.9 | 28.3 | 19.7 |

Combined variant choices: {'M2 length': 'L2_2u', 'M3 length': 'L3_4u', 'M4 width': 'base', 'M1 width': 'W1_110u'}

## Practical reset (8.85 Mohm and 0.1 uA level shift) on the combined circuit

Output during reset 0.8853 V, gate 0.0003 V.
Frame-end drops: 4.464, 12.975, 12.975 mV for 100, 300, 300 accumulated ions: 44.64, 43.25, 43.25 uV per ion; frame 2 alone 42.56 uV per ion.

Simulator runs: 170.
