# Nutrient-Phytoplankton-Zooplankton-Detritus model (Nitrogen and Carbon units)

This model is developed under **Blue-Cloud2026 VLab 3: Carbon-Plankton Dynamics**, simulating the dynamics of carbon and plankton in marine ecosystems using a mechanistic Nutrient-Phytoplankton-Zooplankton-Detritus (**NPZD**) model.

## Table of Contents
- [Overview](#overview)
- [Features](#features)
- [Software](#Software)
- [Installation](#installation)
- [Usage](#usage)
- [Model Details](#model-details)
- [Contributing](#contributing)
- [Scientific output](#Scientific-output)
- [License](#license)
- [Acknowledgements](#Acknowledgements)

## Overview
The **Blue-Cloud2026 VLab 3 Carbon-Plankton Dynamic** project models interactions between nutrients, phytoplankton, zooplankton, and detritus in marine environments. The NPZD model captures essential ecological processes, providing insights into carbon cycling and plankton population dynamics.

## Features
- Mechanistic NPZD model with customizable parameters.
- Simulates carbon flow through different trophic levels.
- Outputs include nutrient concentrations, plankton biomass, and detritus levels over time.
- Configurable for different marine ecosystems and scenarios.

## Software
- Jupyer Notebook
- R

## Installation

Via Blue-Cloud2026 Vlab 3 — The model and dependencies are pre-installed in the virtual environment of the Vlab 3 – Cabon-Plankton dynamics.

## Usage
Via Blue-Cloud2026 Vlab 3
   1. **Open the Jupyter Notebook to start:**
      NPZD.ipynb

   2. **Follow the step-by-step guide in the Jupyter Notebook**

   3. **Copy your results to your workspace**

## Model Details
The NPZD model simulates four key compartments:

   1. **Nutrients (N):** Essential elements required for phytoplankton growth.
   2. **Phytoplankton (P):** Primary producers converting nutrients and carbon into biomass through photosynthesis.
   3. **Zooplankton (Z):** Consumers feeding on phytoplankton.
   4. **Detritus (D):** Non-living organic material resulting from dead organisms and waste products.

Key processes:
   - Nutrient and carbon uptake by phytoplankton.
   - Grazing of phytoplankton by zooplankton.
   - Mortality and decomposition contributing to detritus.
   - Recycling of detritus into nutrients.

Input data:
   - Sea Surface Temperature (SST) in °C
   - Sea Surface Salinity (SSS) in PSU
   - Nutrients (DIN, PO4, SiO4) in mmol N m-3 (or P or Si equivalent)
   - Carbon data (pCO2 atmosphere, windspeed, pH)
   - Threshold values for each parameter per region

   All these data are required.

Calibration and validation data:
   - Chlorophyll-a (µg Chl m-3)
   - Zooplankton abundances (ind m3-p
   - partial pressure of CO2 (pCO2 seawater; µatm)

   Chlorophyll-a data is required, the other are optional. 
   
## Contributing
Please follow these steps:
   1. Fork the repository.
   2. Create a new branch for your feature or bugfix.
   3. Submit a pull request with a clear description of your changes.

## Scientific output
Otero, V.; Pint, S.; Deneudt, K.; De Rijcke, M.; Mortelmans, J.; Schepers, L.; Martin-Cabrera, P.; Sabbe, K.; Vyverman,W.; Vandegehuchte, M.; et al. Pronounced Seasonal and Spatial Variability in Determinants of Phytoplankton Biomass Dynamics along a Near–Offshore Gradient in the Southern North Sea. J. Mar. Sci. Eng. 2023, 11, 1510. https://doi.org/10.3390/jmse11081510

## License
This project is licensed under the Creative Commons CC-BY 4.0 license.

## Acknowledgments
The production of this work has been supported by the Blue Cloud working environment via the Blue Cloud (https://blue-cloud.d4science.org, accessed on 15 March 2023) operated by D4Science.org (www.d4science.org, accessed on 15 March 2023). This work makes use of the LifeWatch data and infrastructure funded by Research Foundation—Flanders (FWO) as part of the Belgian contribution to LifeWatch ESFRI.

The processed data and scripts have also been published in the Blue-Cloud Catalogue (https://data.d4science.org/ctlg/Zoo-Phytoplankton_EOV/spatiotemporal_analysis_of_plankton_drivers_in_the_belgian_part_of_the_north_sea_data, accessed on 10 March 2023), as well as deposited in Zenodo (https://doi.org/10.5281/zenodo.6794084, accessed on 10 March 2023) under a Creative Commons CC-BY 4.0 license. This allows for the use of the data and code under the condition of providing the reference to the original source: Steven Pint and Viviana Otero (2021). Data, scripts, and model output to perform spatiotemporal analysis of the plankton drivers in the Belgian part of the North Sea (dataset). Zenodo: https://doi.org/10.5281/zenodo.6794084, accessed on 10 March 2023. The raw input data J. Mar. Sci. Eng. 2023, 11, 1510 20 of 36 that was obtained from LifeWatch can be accessed through https://doi.org/10.14284/441 and https://doi.org/10.14284/445, accessed on 10 March 2023, or via https://rshiny.lifewatch.be/, accessed on 10 March 2023.

---

For questions or further information, please contact the project maintainer at [steven.pint@vliz.be].

---






