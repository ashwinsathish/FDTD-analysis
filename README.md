# FDTD Analysis of EM Waves in Corroded Reinforced Concrete

## Overview

This project implements a Finite-Difference Time-Domain (FDTD) simulation to analyze electromagnetic wave propagation in corroded reinforced concrete structures. The simulation aims to provide a non-destructive testing method for assessing the integrity and durability of concrete-rebar infrastructures.

## Methodology

### FDTD Algorithm

The simulation uses the Yee algorithm for spatial and temporal interleaving of electric and magnetic fields. The algorithm provides second-order accuracy in finite difference method (FDM) calculations.

Key equations for lossy media:

- Electric field update:
  $E_{x,y,z}^{n+1} = \frac{2\varepsilon - \sigma\Delta t}{2\varepsilon + \sigma\Delta t}E_{x,y,z}^{n} - \frac{2\Delta t}{2\varepsilon + \sigma\Delta t}(\nabla \times H)_{x,y,z}^{n+1/2}$

- Magnetic field update:
  $H_{x,y,z}^{n+1/2} = H_{x,y,z}^{n-1/2} - \frac{\Delta t}{\mu}(\nabla \times E)_{x,y,z}^n$

### Absorbing Boundary Conditions

Mur's 2nd order ABCs are implemented to simulate infinite space propagation:

$E_z \big|_{x=0} = E_z^{\text{past}} \big|_{x=1} + c_1 \left(E_z \big|_{x=1} - E_z^{\text{past}} \big|_{x=0}\right) + c_2 (...) + c_3 (...)$

Where $c_1$, $c_2$, and $c_3$ are constants derived from the wave equation.

### Corrosion Modeling

Corrosion is modeled using a multi-layer dielectric approach:
- Conductive steel layer (rebar)
- Rust coating with specific relative permittivity
- Surrounding concrete

The degree of corrosion is estimated by the percentage loss in rebar weight, with rust occupying up to four times the volume lost by the rebar.

## Key Features

1. 3D FDTD simulation of EM wave propagation
2. Modulated Gaussian pulse excitation
3. Mur's 2nd order Absorbing Boundary Conditions
4. Circular and cuboidal rebar modeling
5. Variable corrosion levels (0% to 30%)
6. FFT analysis for frequency-domain characterization

## Results

### Time-Domain Analysis

- Demonstrated wave propagation through concrete with varying levels of rebar corrosion
- Observed increased delay in wave propagation with increasing corrosion levels

### Frequency-Domain Analysis

- FFT analysis revealed changes in magnitude and phase spectra for different corrosion levels
- Observed shifts in center frequency and variations in 3dB bandwidth with increasing corrosion

### Corrosion-Delay Relationship

- Established an approximately linear relationship between rust thickness and output wave delay
- 3D simulations showed more significant delays compared to 2D simulations, better representing real-world interactions

## References

1. "Analysis of the Electromagnetic Signature of Reinforced Concrete Structures for Non-destructive Evaluation of Corrosion Damage." IEEE TRANSACTIONS ON INSTRUMENTATION AND MEASUREMENT, VOL. 61 (NO. 4), April 2012.

2. "Total-Field Absorbing Boundary Conditions for the Time-Domain Electromagnetic Field Equations." IEEE TRANSACTIONS ON ELECTROMAGNETIC COMPATIBILITY, VOL. 40 (NO. 2), May 1998.
