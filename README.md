# fortan_dcd_analysis
### Open/close, read, analyze dcd trajectories with python frontend and fortran backend

Requirements:
- Intel compiler
- Bash, linux environment
- Python 2
- catdcd https://www.ks.uiuc.edu/Development/MDTools/catdcd/

#### To use:
- Clone/copy the files into your location of choice/cluster


#### Conversion for .xtc to .dcd using catdcd

`
  catdcd -o <dcdname> -xtc <xtcname> <optional>
`

#### To make a python usable module, run

Add the python2, intel, and gcc environment to the top of file "script"

`
  ./script "createmodule"
`

#### Example usage in python code in file : calculations.py
This file demonstates
  - opening dcd file
  - reading dcd header for #atoms, #frames information
  - skips header when reading frames
  - reading frames
  - doing h-bond and first hydration shell calculation with distance and angle cutoff criteria
 
#### Test the python file on a N,P=1bar,T=300K simulation of SPC/E water in 5nm cubic box

To verify, you should see stuff printed on your terminal as contained in file "output"

Make sure you set the python2, intel, and gcc environment

`
 ./script "createmodule" 
`


`
  python calculations.py
`

-----
#### Quick and dirty use

`
module load python/2.7.15
`

`
module load gcc
`

`
module load intel
`

`
git clone https://github.com/ankur9511/fortan_dcd_analysis.git
`

`
cd fortan_dcd_analysis/
`

`
chmod +x ./script
`

`
./script "createmodule"
`

`
python calculations.py
`
