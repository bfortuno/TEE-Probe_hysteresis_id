# TEE-Probe_hysteresis_id

## Purpose
The TEE-Probe_hysteresis_id project aims to optimize the identification of hysteresis in the bending behavior of an ultrasound TEE Probe. The project involves a MATLAB component for initial identification, followed by a Python component for recreating the identified curve.

## Dependencies
### MATLAB
- Communications Toolbox
- Curve Fitting Toolbox
- Global Optimization Toolbox
- Optimization Toolbox

### Python
- scipy
- numpy
- matplotlib

## Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/your-username/TEE-Probe_hysteresis_id.git
   cd TEE-Probe_hysteresis_id
2. Install MATLAB dependencies.

3. Install Python dependencies:
   ```bash
    pip install -r requirements.txt

## Usage
1. Run the MATLAB identification code to generate the PPModel.mat file.
   ```matlab
    run MATLAB_identification.m

2. Use the generated PPModel.mat file in the Python code to recreate the identified curve.
   ```python
    python PPModel.py

## Contributors
- [Benjamin Fortuno](https://github.com/bfortuno)
