# Flood-Resilient Power System Optimisation

**Mitigation of Flood-induced Network Fragmentation and Cascade Failure in Power Systems**

Kinuri Rathnayake | Supervised by Professor Hongjian Sun | Durham University, Department of Engineering

---

## Overview

This repository contains the MATLAB implementation for a probabilistic, multi-objective framework for flood-resilient power system operation. The framework jointly minimises total economic loss and socially weighted criticality loss through a bi-objective DC Optimal Power Flow (DC-OPF) formulation, integrated with spatiotemporal flood hazard modelling and Monte Carlo scenario generation.

---

## Repository Structure

```
flood_resilience/
├── src/
│   ├── main_simulation.m          % Main simulation loop (entry point)
│   ├── solve_DC_OPF.m             % Bi-objective DC-OPF via epsilon-constraint
│   ├── compute_criticality_spatial.m  % Spatial criticality index computation
│   ├── data_loaders.m             % Demand, price, generation mix loaders
│   ├── plot_results.m             % All results figures
│   ├── sensitivity_analysis.m     % OAT sensitivity analysis
│   └── solver_NR.m                % Newton-Raphson AC power flow (external)
├── data/
│   ├── busdata.xlsx               % IEEE 30-bus system bus data
│   ├── linedata.xlsx              % IEEE 30-bus system line data
│   ├── York_Flood_Depth.xlsx      % [30 x 48] flood depth matrix (m)
│   ├── demand2015.xlsx            % UK National Grid demand, 26 Dec 2015
│   ├── marketindexprice2015.xlsx  % APX half-hourly market prices
│   ├── generationmix2015.xlsx     % National Grid generation mix
│   ├── export.geojson             % OSM building centroids, York
│   └── critical.geojson           % Hospital/emergency service locations
└──results/                       % Output figures and tables (generated)
 
```

---

## Requirements

- MATLAB R2022b or later
- Optimisation Toolbox (`linprog`, `fmincon`)


---

## Methods

### Flood Fragility Model
Substation failure probabilities follow a lognormal fragility function:

$$P^{\text{fail}}_i(d_i) = \Phi\left(\frac{\ln d_i - \ln d_{50}}{\beta}\right), \quad d_i > d_0$$

with parameters $d_0 = 0.107$ m, $d_{50} = 1.144$ m, $\beta = 0.4381$ derived from the 2015 York Boxing Day flood event.

### Composite Criticality Index
A time-varying composite criticality index combines three spatial components:

$C_i(t) = 0.4 \cdot \frac{d_{i,t}}{d_{\max}} + 0.3 \cdot Y_i + 0.3 \cdot \frac{\rho_i}{\rho_{\max}}$

where $d_{i,t}/d_{\max}$ is normalised flood exposure, $Y_i \in \{0,1\}$ is proximity to critical infrastructure (hospitals, police, fire stations) within 500 m, and $\rho_i/\rho_{\max}$ is normalised building density within 700 m.

### Bi-objective DC-OPF
Two competing objectives are minimised simultaneously:

**Economic objective** (all buses, GBP):
$f_{econ}$ = \sum_{g} C_g P_g + \sum_{i} \text{VoLL}_i (1-\alpha_i) P_{D,i}

**Criticality objective** (quadratic penalty, GBP):
$f_{crit}$ = \sum_{i} \text{VoLL}_i \cdot C_i \cdot (1-\alpha_i)^2 \cdot P_{D,i}

The $\epsilon$-constraint method traces the Pareto front; the knee-point solution minimises the normalised Euclidean distance to the ideal point.

### VoLL Assignments

| Sector | Bus IDs | VoLL (£/MWh) |
|--------|---------|--------------|
| Emergency | 1, 12, 22 | 50,000 |
| Commercial | 2,3,6,10,15,18,19,24,25 | 15,000 |
| Residential | 4,5,7,11,17,20,23,27,28,29,30 | 8,000 |
| Public | 8,9,13,14,16,21,26 | 2,000 |

---
