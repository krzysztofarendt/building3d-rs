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
