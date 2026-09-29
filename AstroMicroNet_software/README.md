# Software

Integrated MATLAB GUI V11 with its unchanged shared Python detection backend. The real-cell example and its dependent GUI regression fixtures were excluded for this repository. The supplied example is synthetic.

## Python

From the repository root, follow the main README. The input is either `R` (N by N connectivity) or `X` (T by N processed traces) in NPZ or MATLAB v7 MAT format. Optional `roi_ids` must be strictly increasing. Final assignments are in `final`; `before_gate` is diagnostic only. An empty final core is valid. Python array indices are zero-based; MATLAB `coreIndices` are one-based.

Simulation examples explicitly set rho=0.28. Real-data defaults use rho=0.20, support threshold 0.20, 1,024 random starts plus 17 matrix starts, and 499 edge-weight permutations. Candidate capacity is not a forced final core count. The original backend files were copied without algorithm changes.

## MATLAB

Previously tested upstream with MATLAB R2024a on Windows and Python 3.13, NumPy 2.2.5 and SciPy 1.15.3. Image and signal functions require Image Processing Toolbox, Signal Processing Toolbox, Statistics and Machine Learning Toolbox, and Curve Fitting Toolbox (`smooth`).

Install `requirements.txt` into the Python interpreter you will use. Set the MATLAB current folder to this directory:

```matlab
app = AstroMicroNet('C:/path/to/python.exe');
% Optional historical interfaces:
% app = AstroMicroNet('C:/path/to/python.exe', 'v1');
% app = AstroMicroNet('C:/path/to/python.exe', 'guide');
```

The GUI retains the 5 Hz processing convention; local AVI inputs must already be sampled at 5 Hz. Load the local file, delineate structures and microdomains, then select Analyze Signals. The user is responsible for supplying imaging data they are authorized to use.

For an existing matrix:

```matlab
S = load('your_network.mat');
opts = struct('pythonExecutable','C:/path/to/python.exe', ...
              'inputKind','R','rho',0.20,'seedKey','my_cell');
result = astromicronet_analyze(S.R, opts);
```

Both interfaces use the same final labels for display, statistics and export. The historical Distance Threshold control zeros connectivity values below the entered value. Spatial calibration defaults to 0.187 micrometres/pixel and must match the user's recording. This repository does not supply a real imaging fixture or claim a new full GUI regression test.
