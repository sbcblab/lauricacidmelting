# **Phase Change Material Analysis**
This repository contains the code, data, and documentation for analyzing the melting behavior of phase change materials (PCMs) using the Lattice Boltzmann Method (LBM). PCMs absorb or release thermal energy during phase transitions and are widely used in thermal energy storage systems.

## **Overview**
The purpose of this project is to simulate and analyze the melting process of PCMs inside a controlled environment using LBM. By leveraging the OpenLB library, this project aims to model heat transfer, phase transitions, and other critical phenomena to optimize PCM performance for energy storage applications.

### **Key Features**
- Simulation of the melting process of PCMs.
- Heat transfer modeling using LBM with enthalpy-based methods.
- Grid-based analysis for temperature distribution and phase states.
- Analysis of PCM melting behavior in rectangular cavities with heating elements.

## **Simulation Setup**
### **Geometry**
The domain is a rectangular cavity with the following specifications:
- **Height**: 120 mm
- **Width**: 50 mm
- **Configuration**: 
  - A horizontal fin is attached to the middle of the right-hand wall.
  - The right wall is heated, while the other walls are thermally insulated.

### **Phase Change Material**
- Material: Lauric Acid
- Melting range: 42–47 °C
- Thermal properties:
  - Latent heat: 187210 kJ/kg
  - Thermal conductivity: 0.15 W/mK

### **Simulation Techniques**
- **Numerical Method**: Lattice Boltzmann Method (LBM).
- **Library**: [OpenLB](https://www.openlb.net/).
- **Grid Type**: Uniform structured grid for simulation accuracy.
- **Boundary Conditions**:
  - Heated wall: Constant heat flux or temperature.
  - Insulated walls: No heat transfer.

## **Results**
The analysis provides:
- **Temperature Distribution**: A grid-based map showing temperature variations.
- **Phase Change Progress**: Visualization of liquid and solid fractions over time.

## **Usage**
### **Prerequisites**
- OpenLB installed on your system ([Installation Guide](https://www.openlb.net/installation/)).
- A C++ compiler compatible with OpenLB.
- Python 3.x (optional, for data post-processing and visualization).
- Required Python libraries: `numpy`, `matplotlib`, `pandas` (optional).

### **Steps to Run the Simulation**
1. Clone this repository:
   ```bash
   git clone https://github.com/username/pcm-analysis.git
   cd pcm-analysis
