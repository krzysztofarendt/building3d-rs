# BESTEST: Building Energy Simulation Test

## What is BESTEST?

BESTEST (Building Energy Simulation Test) is a standardized diagnostic methodology
for testing, validating, and improving building energy simulation software. It was
developed at the National Renewable Energy Laboratory (NREL) by Ron Judkoff and
Joel Neymark, originally published in 1995 under the IEA Solar Heating and Cooling
Programme Task 12. It was later codified as **ANSI/ASHRAE Standard 140** (first
published 2001, with updates in 2004, 2007, 2011, 2014, 2017, 2020, and 2023).

### Three-tier validation framework

1. **Analytical verification** -- test cases with exact mathematical solutions
   (e.g., steady-state conduction through a known wall).
2. **Empirical validation** -- comparison to measured data from instrumented test
   cells (e.g., ETNA test cells in France, LBNL's FLEXLAB, ORNL's Flexible
   Research Platform).
3. **Comparative testing** -- the dominant method: multiple reference programs run
   the same cases, and their result spread defines an acceptance range.

The acceptance range reflects the consensus of state-of-the-art programs, not
ground truth. If all programs share a systematic error, BESTEST will not detect it.
BESTEST is best understood as a diagnostic and quality-assurance tool -- it catches
bugs and outlier behavior, but does not guarantee a match to physical reality.

### Test case structure

BESTEST defines a series of carefully specified single-zone buildings that
progress from simple to realistic. Each case isolates specific physical phenomena
so that if a program disagrees with the reference range, a diagnostic logic flow
can identify which algorithm is responsible.

- **Case 600 series** -- lightweight (wood-frame) single-zone 8 m x 6 m x 2.7 m
  building with south-facing windows. Denver (DRYCOLD) TMY weather. Thermostat
  deadband: 20 C heating / 27 C cooling. Variants: 610 (overhang shading),
  620 (east/west windows), 630 (east/west shading), 640 (night setback),
  650 (night ventilation), 600FF (free-float, no HVAC).

- **Case 900 series** -- identical geometry but heavyweight construction (high
  thermal mass). Key diagnostic: Case 900 heating is 34--38% of Case 600,
  demonstrating the thermal mass effect. The free-float variant (900FF) shows
  temperature swings of +/-6 C vs. +/-25 C for 600FF.

---

## Programs tested

### 1995 original (IEA BESTEST, 8 programs)

| Program        | Organization                     |
|----------------|----------------------------------|
| BLAST 3.0      | US Army CERL                     |
| DOE-2.1D       | NREL                             |
| ESP            | De Montfort University (DMU)     |
| SRES/SUN       | NREL                             |
| SERIRES 1.2    | Building Research Establishment  |
| S3PAS          | University of Seville            |
| TASE           | VTT Finland                      |
| TRNSYS 13.1    | BEL / BRE                        |

### Early 2000s additions (EnergyPlus validation report)

| Program               | Organization |
|-----------------------|--------------|
| BLAST 3.0-334         | US-IT        |
| DOE-2.1E              | LBNL         |
| DOE-2.1E (RevWindow)  | LBNL         |
| EnergyPlus 1.2.0      | US DOE       |

### 2020 update (ASHRAE 140-2020, 7 programs)

| Program          | Version   | Organization          |
|------------------|-----------|-----------------------|
| BSIMAC           | 9.0.74    | --                    |
| CSE              | 0.861.1   | California Energy Commission |
| DeST             | 2.0       | Tsinghua University   |
| EnergyPlus       | 9.0.1     | US DOE / NREL         |
| ESP-r            | 13.3      | University of Strathclyde |
| NewHASP          | --        | (excluded from standard due to input deviations) |
| TRNSYS           | 18.01     | University of Wisconsin / TESS |

---

## Results: Case 600 (lightweight construction)

### 1995 original results (Denver TMY weather)

Individual per-program values are available in the NREL/TP-472-6231 report and
ASHRAE Standard 140 companion spreadsheet (Std140_TF_Results.xlsx). The freely
available summary statistics are:

| Metric               | Min     | Max     | Mean    | EnergyPlus |
|----------------------|---------|---------|---------|------------|
| Annual heating (MWh) | 4.296   | 5.709   | 5.046   | 4.673      |
| Annual cooling (MWh) | 6.137   | 8.448   | 7.053   | 6.792      |
| Peak heating (kW)    | 3.437   | 4.354   | 3.952   | 3.838      |
| Peak cooling (kW)    | 5.965   | 7.188   | 6.535   | 6.664      |

**Spread (max-min as % of min):**
- Heating: 32.9%
- Cooling: 29.8% (37.6% including outliers up to 8.448 MWh)

### 2020 updated results (Denver TMY weather, 6 programs)

| Metric               | Min     | Max     | Mean    |
|----------------------|---------|---------|---------|
| Annual heating (MWh) | 3.993   | 4.504   | 4.213   |
| Annual cooling (MWh) | 5.432   | 6.162   | 5.856   |
| Peak heating (kW)    | 3.020   | 3.359   | 3.184   |
| Peak cooling (kW)    | 5.422   | 6.481   | 6.024   |

**Spread (max-min as % of min):**
- Heating: 12.8%
- Cooling: 13.4%

---

## Results: Case 900 (heavyweight construction)

### 1995 original results (Denver TMY weather)

| Metric               | Min     | Max     | Mean    | EnergyPlus |
|----------------------|---------|---------|---------|------------|
| Annual heating (MWh) | 1.170   | 2.041   | 1.649   | 1.350      |
| Annual cooling (MWh) | 2.132   | 3.669   | 2.826   | 2.424      |
| Peak heating (kW)    | 2.850   | 3.797   | 3.452   | 3.240      |
| Peak cooling (kW)    | 2.888   | 3.932   | 3.460   | 3.359      |

**Spread (max-min as % of min):**
- Heating: 74.4%
- Cooling: 60.2% (72.1% including outliers up to 3.669 MWh)

### 2020 updated results (Denver TMY weather, 6 programs)

| Metric               | Min     | Max     | Mean    |
|----------------------|---------|---------|---------|
| Annual heating (MWh) | 1.379   | 1.814   | 1.626   |
| Annual cooling (MWh) | 2.267   | 2.714   | 2.467   |
| Peak heating (kW)    | 2.443   | 2.778   | 2.591   |
| Peak cooling (kW)    | 2.556   | 3.376   | 2.975   |

**Spread (max-min as % of min):**
- Heating: 31.5%
- Cooling: 19.7%

---

## Summary: spread across programs

| Case | Metric          | 1995 Min | 1995 Max | 1995 Spread | 2020 Min | 2020 Max | 2020 Spread |
|------|-----------------|----------|----------|-------------|----------|----------|-------------|
| 600  | Heating (MWh)   | 4.296    | 5.709    | 32.9%       | 3.993    | 4.504    | 12.8%       |
| 600  | Cooling (MWh)   | 6.137    | 8.448    | 37.6%       | 5.432    | 6.162    | 13.4%       |
| 900  | Heating (MWh)   | 1.170    | 2.041    | 74.4%       | 1.379    | 1.814    | 31.5%       |
| 900  | Cooling (MWh)   | 2.132    | 3.669    | 72.1%       | 2.267    | 2.714    | 19.7%       |

The 2020 update significantly tightened results, reflecting 25 years of improved
test specifications and corrected software errors.

## BESTEST vs. ASHRAE 140 conformance

BESTEST and ASHRAE Standard 140 are related but distinct:

- **BESTEST** (1995) is the diagnostic methodology -- the test cases, building
  specs, and the idea of running multiple programs and comparing results. The
  original report had no pass/fail criteria. It showed the spread across programs
  and provided diagnostic flowcharts to help developers find bugs. It is
  educational and diagnostic in nature.

- **ASHRAE Standard 140** (2001+) codified BESTEST into a formal standard with
  normative **acceptance ranges**. These ranges are deliberately wider than the
  actual program spread to allow for legitimate modeling differences. A program
  must fall within these ranges to "conform." ASHRAE 140 conformance is required
  by energy codes such as ASHRAE 90.1 and IECC -- if a simulation tool does not
  pass, it cannot be used for building permit compliance in the US.

The acceptance ranges compared to the actual program spread (2020 round):

| Case | Metric  | Program spread      | ASHRAE acceptance range |
|------|---------|---------------------|-------------------------|
| 600  | Heating | 3.99 -- 4.50 MWh    | 3.75 -- 4.98 MWh        |
| 600  | Cooling | 5.43 -- 6.16 MWh    | 5.00 -- 6.83 MWh        |
| 900  | Heating | 1.38 -- 1.81 MWh    | 1.04 -- 2.28 MWh        |
| 900  | Cooling | 2.27 -- 2.71 MWh    | 2.35 -- 2.60 MWh        |

## building3d test tolerances

The building3d test suite (`tests/bestest_energy_suite.rs`) compares against
EnergyPlus/OpenStudio reference values from the NREL BESTEST-GSR workflow, using
a **Boston Logan TMY3 EPW** (not the standard Denver weather). The tolerances are:

| Case | Heating tolerance | Cooling tolerance | Notes                          |
|------|-------------------|-------------------|--------------------------------|
| 600  | +/- 10%           | +/- 15%           | 1R1C steady-state model        |
| 900  | +/- 25%           | +/- 25%           | 2R2C transient, approximate    |

These tolerances are comparable to the 2020 inter-program spread (13--32%) and
are tighter than the original 1995 spread (30--74%). Since the comparison is
against a single EnergyPlus run with the same weather file, these are regression
guards rather than ASHRAE conformance tests.

---

## Sources

- Judkoff, R. and Neymark, J. (1995). International Energy Agency Building
  Energy Simulation Test (IEA BESTEST) and Diagnostic Method. NREL/TP-472-6231.
  https://www.osti.gov/biblio/90674
- EnergyPlus Testing with ANSI/ASHRAE Standard 140-2001. LBNL.
  https://simulationresearch.lbl.gov/dirpubs/epl_bestest_ash.pdf
- Neymark, J. et al. (2020). Update of ASHRAE Standard 140 Section 5.2 and
  Related Sections. ANL.
  https://publications.anl.gov/anlpubs/2020/05/158451.pdf
- ASHRAE Standard 140-2023 Resource Files.
  https://data.ashrae.org/standard140/AccompanyingSection7.html
- Modelica Buildings Library BESTEST validation.
  https://simulationresearch.lbl.gov/modelica/releases/v7.0.0/help/Buildings_ThermalZones_Detailed_Validation_BESTEST_Cases6xx.html
