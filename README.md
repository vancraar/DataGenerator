# Programming Project: Galaxy-Crash - Data Generation

This repository contains a Python3 script to generate simulation data.

## Dependencies

Required Python3 packages:

- [dataclass_csv](https://pypi.org/project/dataclass-csv/)
- [argparse](https://pypi.org/project/argparse/)

They can be installed using:

```shell
pip3 install -r requirements.txt
```

## Usage

The script can be used to generate data sets containing a large and a small galaxy where the small galaxy orbits the larger one and crashing into it during the simulation.
The `RATIO` determines how many of the `NUM_PARTICLES` are associated with the small galaxy.
For example, with `-n 1000` and `-r 0.2`, 800 particles are assigned to the large galaxy and 200 particles to the small galaxy.


```shell
python3 generate_data.py --help
usage: generate_data [-h] -o OUTPUT -n NUM_PARTICLES [-r RATIO]

Generate a data set consisting of two galaxies colliding with each other.

options:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        the output filename containing the Cartesian State Vectors in CSV format
  -n NUM_PARTICLES, --num_particles NUM_PARTICLES
                        the number of particles in the data set
  -r RATIO, --ratio RATIO
                        the ratio of particles between the large and small galaxy (0.0 = all particles are located in the large galaxy)

```

The resulting csv file contains 8 columns:

```
id,mass,pos_x,pos_y,pos_z,vel_x,vel_y,vel_z
```

- the unique ID of the particle to identify it in the ground truth file
- the mass of the particle
- the position of the particle in the 3-dimensional space
- the velocity of the particle in the 3-dimensional space

