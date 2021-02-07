# fortan_dcd_analysis

To create a python importable module:

Requirements:
- Intel compiler
- Bash, linux environment
- Python 2

To use:
- Clone the files into your location of choice

Run

`
  ./script "createmodule"
`

Example python code in file:
- calculations.py
- This file demonstates
  - opening dcd file
  - reading dcd header for #atoms, #frames information
  - skips header when reading frames
  - reading frames
  - doing h-bond and first hydration shell calculation with distance and angle cutoff criteria
 
Test the python file on a N,P=1bar,T=300K simulation of SPC/E water in 5nm cubic box

`
 ./script "createmodule" 
  python calculation.py
`

Conversion for .xtc to .dcd using catdcd (https://www.ks.uiuc.edu/Development/MDTools/catdcd/)

`
  catdcd -o <dcdname> -xtc <xtcname> <optional>
`
