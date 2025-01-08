# **Phase Change Material Analysis**
This repository contains the code, data, and documentation for analyzing the behavior of phase change materials (PCMs) under various thermal conditions. PCMs are substances that absorb or release thermal energy during phase transitions, such as melting or solidification, and are widely used in thermal energy storage systems.

## **Overview**
The purpose of this project is to simulate and analyze the thermal behavior of PCMs inside a controlled environment using computational fluid dynamics (CFD). By leveraging advanced simulation techniques, the project aims to model heat transfer, phase transitions, and other critical phenomena to optimize PCM performance for energy storage applications.

### **Key Features**
- Simulation of melting and solidification processes.
- Heat transfer modeling using the enthalpy-porosity approach.
- Grid-based analysis for temperature distribution and phase state.
- Analysis of PCM behavior in rectangular cavities with heating elements.

## **Simulation Setup**
### **Geometry**
The domain is a rectangular cavity with the following specifications:
- **Height**: 120 mm
- **Width**: 50 mm
- **Configuration**: 
  - A horizontal fin is attached to the middle of the right-hand wall.
  - The right wall is heated, while the other walls are thermally insulated.

### **Phase Change Material**
- Material: [Specify PCM used, e.g., paraffin wax, salt hydrates].
- Melting range: [Provide range, e.g., 45–55 °C].
- Thermal properties:
  - Latent heat: [Value, e.g., 200 kJ/kg].
  - Thermal conductivity: [Value, e.g., 0.2 W/mK].

### **Simulation Techniques**
- **Numerical Method**: Finite Volume Method (FVM).
- **Solver**: [Specify solver, e.g., Ansys Fluent, OpenFOAM].
- **Grid Type**: Structured grid for enhanced accuracy.
- **Boundary Conditions**:
  - Heated wall: Constant heat flux or temperature.
  - Insulated walls: No heat transfer.

## **Results**
The analysis provides:
- **Temperature Distribution**: A grid-based map showing temperature variations.
- **Phase Change Progress**: Visualization of liquid and solid fractions over time.
- **Heat Flux Analysis**: Quantification of heat transfer rates.

## **Usage**
### **Prerequisites**
- Python 3.x (for data processing and visualization).
- Required libraries: `numpy`, `matplotlib`, `pandas`, [add any CFD-specific libraries or tools].

### **Steps to Run the Simulation**
1. Clone this repository:
   ```bash
   git clone https://github.com/username/pcm-analysis.git
   cd pcm-analysis
