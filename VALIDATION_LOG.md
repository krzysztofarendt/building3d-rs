# Validation Progress Log

BESTEST Cases 600 & 900 vs OpenStudio/EnergyPlus reference (Boston Logan TMY3 EPW).

Reference annual values:
- Case 600: heating = 4324.8 kWh, cooling = 6044.1 kWh
- Case 900: heating = 1661.2 kWh, cooling = 2498.2 kWh

## Results by Commit

| Commit | Description | 600 Heat | 600 Cool | 900 Heat | 900 Cool | Notes |
|--------|-------------|----------|----------|----------|----------|-------|
| 536a19f | Merge heat into fvm (dynamic convection infra + MRT fix + pure air cap) | -3.8% | -53.2% | +87.6% | -43.7% | Convection defaults to Fixed; infra only, not enabled |

### Ablation: enabling dynamic convection + solar-to-air (on top of 536a19f)

| Change | 600 Heat | 600 Cool | 900 Heat | 900 Cool |
|--------|----------|----------|----------|----------|
| Baseline (Fixed h_in=3.0, Fixed h_out) | -3.8% | -53.2% | +87.6% | -43.7% |
| + TARP interior convection | +12.5% | -33.0% | +101.9% | -24.5% |
| + DOE-2 exterior convection | +13.0% | -32.9% | +102.3% | -24.8% |
| + transmitted_solar_to_air_fraction=0.4 | +13.1% | -25.0% | +105.8% | -5.3% |

**Observations:**
- TARP is the biggest single lever: +20pp improvement on cooling, but worsens heating by ~15pp
- DOE-2 exterior convection is negligible (~0.3pp)
- Solar-to-air=0.4 improves cooling further (Case 900 cooling from -24.5% to -5.3%)
- Heating over-prediction is the dominant remaining gap (especially Case 900: +106%)

| bf4ae52 | Enable TARP + DOE-2 + solar-to-air=0.4 | +13.1% | -25.0% | +105.8% | -5.3% | All tests pass; Case 900 cooling nearly matched |
| e9c8993 | Fix TARP h_min_iso clamp + insulation thickness + TARP radiative bypass | -0.7% | -7.2% | +42.7% | +26.8% | Case 600 nearly perfect; Case 900 much improved |

### Bug fixes applied in latest commit

**Bug 1: TARP clamped to h_min_iso = 7.69 W/m²K**
- `interior_convection_h()` in `convection.rs` applied `h.max(h_min_iso)` where h_min_iso = 1/R_si = 7.69
- ISO R_si = 0.13 is a COMBINED (convective + radiative) coefficient; TARP returns convection-only values (typically 1.5-4.0 W/m²K)
- The clamp effectively made TARP return 7.69 always, identical to pre-TARP behavior
- Fix: removed h_min_iso clamp for TARP, keeping only H_MIN_FLOOR = 0.1 as physical minimum

**Bug 2: Case 600 fiberglass over-insulated by 16%**
- Fiberglass thickness was 0.0776m (RSI 1.94 * k=0.04) instead of correct 0.066m
- BESTEST-GSR reference (NREL OpenStudio) specifies RSI 1.65 for fiberglass layer → 0.066m at k=0.04
- This artificially reduced envelope losses, suppressing heating demand

**Bug 3: Linearized radiation through MRT loses energy**
- First attempt to fix TARP added h_rad = 4εσT³ ≈ 5.1 W/m²K as separate radiative component
- This routed heat through MRT, but simplified area-weighted MRT model cannot properly return it to air
- Results worsened: 600H +21.6%, 600C -38.1%, 900H +134.6%, 900C -22.1%
- Fix: TARP mode bypasses radiative exchange entirely (h_conv = h_total = TARP, t_eff = T_air)
- This is physically justified: TARP is designed as a total interior surface coefficient for the convective path; adding a separate radiation model requires proper view-factor calculations

### View-factor interior longwave radiation

| Change | 600 Heat | 600 Cool | 900 Heat | 900 Cool |
|--------|----------|----------|----------|----------|
| Baseline (e9c8993, TARP bypass) | -0.7% | -7.2% | +42.7% | +26.8% |
| + View factors (VF on, solar-to-walls=false) | +9.9% | +4.1% | +63.7% | +50.3% |
| + Solar to FVM walls (VF on, solar-to-walls=true) | +13.3% | -0.8% | +74.6% | +41.8% |
| + Solar-to-air=0.0 (all solar to mass, walls=false) | +9.2% | +3.6% | +53.5% | +42.0% |
| + Solar-to-air=0.0, solar-to-walls=true | +13.0% | -5.7% | +63.5% | +20.7% |

**View-factor model:**
- Monte Carlo ray-cast computation of geometric view factors (F_ij) between interior zone surfaces
- Uniform linearized h_rad = eps * 4 * sigma * T_mean^3 (~5.1 W/m²K at 20C)
- Per-surface MRT: T_mrt_i = sum_j(F_ij * T_j) (instead of area-weighted MRT)
- Energy conservation guaranteed by reciprocity: sum(A_i * h_r * (T_i - T_mrt_i)) = 0
- Interior BC: h_total = h_conv(TARP) + h_rad, t_eff = (h_conv * T_air + h_rad * T_mrt_i) / h_total
- Combined h_in ~ 7 W/m²K (TARP ~2.5 + radiation ~4.6), vs old TARP-only ~2.5

**Key bugs fixed during implementation:**
- Ray normals: building polygons have outward-facing normals; VF code must negate them for interior rays
- Energy conservation: air gain must use h_total * (T_surf - t_eff), not h_conv * (T_surf - T_air), because radiation nets to zero by reciprocity but omitting it from zone balance loses energy

**Analysis of view-factor impact:**
- Case 600 (lightweight): results remain good (+10% heat, +4% cool)
- Case 900 (heavyweight): worsened (+64% heat, +50% cool) vs baseline (+43%, +27%)
- Root cause: the old TARP-only code (h_in ~ 2.5 W/m²K) created "phantom insulation" (R_si = 0.40 vs correct 0.13 m²K/W) that accidentally compensated for other model deficiencies
- With correct h_in ~ 7 W/m²K: R_si drops to 0.14, increasing effective U-value ~7%, which disproportionately affects Case 900 because its reference demand is smaller (1661 vs 4325 kWh)

**Solar-to-air fraction (changed from 0.4 to 0.0):**
- Previous: 40% of transmitted solar went directly to zone air, 60% to mass surfaces
- EnergyPlus FullInteriorAndExterior distributes ALL solar to surfaces (0% to air)
- With solar-to-air=0.4: solar bypasses mass storage, causing weaker thermal lag
- With solar-to-air=0.0: all solar flows through mass, improving heavyweight buffering by ~10pp
- Case 600 barely affected (lightweight mass releases heat quickly anyway)
- Case 900 improved: +64% → +54% heating, +50% → +42% cooling

**Solar to FVM walls (rejected):**
- `ConvectiveWithFluxToDomain` injects 100% of solar flux into wall domain
- Solar deposited on exterior walls partially conducts through to outside = heat loss
- Heating worsened for both cases; cooling improved slightly (solar leakage reduces cooling load)
- Reverted: `distribute_transmitted_solar_to_fvm_walls = false`

| 5d9645f | View-factor interior radiation (VF on, solar-to-walls=false) | +9.9% | +4.1% | +63.7% | +50.3% | Physics improved; Case 900 worse due to removed error compensation |
| (pending) | Solar-to-air=0.0 (all solar to mass surfaces) | +9.2% | +3.6% | +53.5% | +42.0% | Matches EnergyPlus FullInteriorAndExterior solar distribution |
