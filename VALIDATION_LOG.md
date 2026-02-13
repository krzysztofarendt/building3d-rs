# Validation Progress Log

BESTEST Cases 600 & 900 vs OpenStudio/EnergyPlus reference (Boston Logan TMY3 EPW).

Reference annual values:
- Case 600: heating = 4324.8 kWh, cooling = 6044.1 kWh
- Case 900: heating = 1661.2 kWh, cooling = 2498.2 kWh

## Results by Commit

| Commit | Description | 600 Heat | 600 Cool | 900 Heat | 900 Cool | Notes |
|--------|-------------|----------|----------|----------|----------|-------|
| 536a19f | Merge heat into fvm (dynamic convection infra + MRT fix + pure air cap) | -3.8% | -53.2% | +87.6% | -43.7% | Convection defaults to Fixed; infra only, not enabled |
