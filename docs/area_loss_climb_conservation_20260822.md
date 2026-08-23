# GB-area-loss climb source and conservation derivation

## Source term

For a two-dimensional calculation with unit out-of-plane thickness, the
grain-boundary area is represented by the physical GB length.  At saved solver
state \(n\),

\[
V_{\mathrm{ex}}^n
=\alpha\sum_\beta v_{\mathrm{ex},\beta}^{\mathrm{GB}}L_\beta^n,
\qquad
v_{\mathrm{ex}}^{\mathrm{GB}}
=\delta_{\mathrm{GB}}\left(1-\frac{\rho_{\mathrm{lattice}}}{\rho_{\mathrm{GB}}}\right).
\]

The implementation also accepts a directly calibrated
`excess_volume_per_area`; this takes precedence over the density expression.
The released volume and signed vacancy-equivalent quota are

\[
\Delta V_{\mathrm{ex}}^{n+1}
=\max\!\left(V_{\mathrm{ex}}^n-V_{\mathrm{ex}}^{n+1},0\right),
\qquad
\Delta N_{\mathrm{req}}^{n+1}
=\frac{\Delta V_{\mathrm{ex}}^{n+1}}{\Omega_{\mathrm{def}}}.
\]

Thus a rigidly translating straight GB has \(L^{n+1}=L^n\) and creates no
climb demand.  Shortening creates positive demand, while lengthening does not
create an equal-and-opposite positive source.  Neither swept area
\(L|\Delta x_n|\) nor the absolute GB-length change appears in this source.

In the executed 900 K campaigns, the direct calibration is
`excess_volume_per_area = 0.01`, `point_defect_formation_volume = 0.02`, and
`alpha = 1`.  Therefore each unit decrease of total physical GB length releases
0.01 units of excess volume and requires 0.5 vacancy-equivalent units.  No
`delta_gb`, `rho_lattice`, or `rho_gb` values are silently substituted when the
direct calibration is present.

## Material-wide signed ledger

Vacancy and interstitial inventories are stored separately; a positive signed
quota denotes vacancy-equivalent content and a negative quota denotes
interstitial-equivalent content.  For either species \(s\), the invariant is

\[
N_{\mathrm{required},s}
=N_{\mathrm{GB},s}+N_{\mathrm{TJ},s}+N_{\mathrm{external},s}
+q_{\mathrm{active},s}+q_{\mathrm{retired},s}.
\]

Summing the two nonnegative species ledgers gives the reported global balance

\[
N_{\mathrm{required}}
=N_{\mathrm{accommodated}}^{\mathrm{GB}}
+N_{\mathrm{accommodated}}^{\mathrm{TJ}}
+N_{\mathrm{external}}+q_{\mathrm{stored}}.
\]

The invariant is inductive under every ledger operation:

1. A source increment adds the same amount to `required` and to one active GB
   domain (or the material reservoir when a disappearing domain cannot be
   attributed locally).
2. A split partitions a parent's inventory among children with normalized
   nonnegative weights, so the stored sum is unchanged.
3. A merge removes the parent entries and adds their exact sum to the child.
4. Topological retirement transfers active inventory to the retired material
   reservoir instead of deleting it.
5. A sink moves the accepted amount from active/retired storage to exactly one
   of the GB, TJ, or external accommodated counters.  Accepted quota is capped
   by available same-sign inventory, and competing completions are resolved in
   continuous completion-time order, so it cannot be consumed twice.
6. Checkpoint state serializes the complete signed ledger and reloads it only
   after verifying the invariant.

The runtime assertion uses a relative tolerance of \(10^{-10}\) against the
total required inventory.  Persisted `defect_inventory.csv` records the total,
species residuals, maximum absolute residual, active and retired storage, and
the GB/TJ/external partition at every output checkpoint.

## Sink variants

`C_GB` permits only a signed climbing-disconnection/connected-GB sink.  A
completion selects a disconnection whose point-defect sign matches the
available material inventory; TJs cannot consume quota.

`C_GBTJ` retains that GB sink as a competing alternative and adds a TJ sink.
The TJ candidate must have same-sign flux, be adjacent to a demand-bearing GB
(or the retired material reservoir), admit a least-squares common TJ velocity
within the configured compatibility tolerance, and pass either strict Burgers
closure or the configured finite residual-energy penalty.  A stochastic stage
clock only creates a candidate; it cannot by itself authorize TJ climb.

The executed reference parameters are:

| Parameter | Value |
|---|---:|
| `climb_trigger_quota` | 0.25 |
| `climb_release_quota` | 1.0 |
| `tj_sink_release_quota` | 1.0 |
| GB nucleation / exchange / transport barriers | 0.45 / 0.55 / 0.65 eV |
| TJ nucleation / exchange / transport barriers | 0.45 / 0.55 / 0.65 eV |
| all stage prefactors | \(10^5\) |
| `tj_sink_compatibility_tolerance` | 0.25 |
| `tj_sink_step_length` | 1.0 |
| `tj_sink_burgers` | 0.25 |
| `tj_residual_stiffness_ev` | 1.0 eV |
| `free_volume_stiffness` | 0.05 |

The fast/slow discriminating cases change the three serial barriers together
to 0.35/0.45/0.55 eV and 0.55/0.65/0.75 eV, respectively; they do not change
the area-loss source calibration.
