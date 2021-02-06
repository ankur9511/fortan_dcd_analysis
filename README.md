# fortan_dcd_analysis

To create a python importable module:

Requirements:
- Intel compiler
- Bash, linux environment

To use:
- Clone the files into your location of choice

Run

`
  ./script "createmodule"
`

Example python code in file:
- calculations.py

Conversion for .xtc to .dcd using catdcd (https://www.ks.uiuc.edu/Development/MDTools/catdcd/)

`
  catdcd -o <dcdname> -xtc <xtcname> <optional>
`
