# N-Body Data Generation

This repository contains two Python3 script to generate N-body initial conditions for one or more configurable galaxies.

Each galaxy consists of one central black hole and a configurable number of disk particles. Galaxies can either be placed explicitly or orbit around another galaxy. Orbiting galaxies can be configured as stable satellites or as crash/merger candidates.

## Dependencies

The generators only uses Python standard library modules. No additional packages are required.

## Disk Generator 

Generate a data set from a JSON scenario config:

```shell
python3 generate_data.py --config examples/collision_three_galaxies.json
```

The output file can be set in the config via `output` or overridden on the command line:

```shell
python3 generate_data.py --config examples/collision_three_galaxies.json --output data.csv
```

Show CLI help:

```shell
python3 generate_data.py --help
```

## Plummer Generator

`generate_plummer_data.py` generates one or more spherical Plummer systems instead of disk galaxies. It uses a separate config shape based on `systems` and does not create a central black hole particle.

Generate the included Plummer collision example:

```shell
python3 generate_plummer_data.py --config examples/plummer_collision.json
```

Each system requires `name`, `particles`, `total_mass`, and `scale_radius`. Root systems also require `position` and `velocity`; orbiting systems use the same `orbit` modes as the disk generator. `cutoff_radius` is optional and limits the sampled Plummer tail.

By default each system is recentered and momentum-corrected before its global `position` and `velocity` are applied. These defaults can be disabled per system with `correct_center_of_mass: false` or `correct_momentum: false`.

Use the disk generator when you want structured rotating disk galaxies; use the Plummer generator when you want spherical, approximately equilibrium particle systems.

## Disk Generator vs. Plummer Generator

This repository contains two independent generators for different kinds of N-body initial conditions.

`generate_data.py` creates disk-like galaxy setups. Each galaxy contains one central black hole particle and a configurable number of disk particles. Disk particles are placed in a rotating disk with configurable radius, thickness, orientation, particle mass range, and orbit settings. This generator is useful for visual galaxy collision scenarios with clear disk structures.

`generate_plummer_data.py` creates spherical Plummer systems. It does not create a central black hole particle. Instead, all particles are sampled from a self-gravitating Plummer distribution and usually have equal mass. Positions are isotropic, velocities are sampled from the Plummer distribution function, and each system is recentered and momentum-corrected by default. This generator is better suited for approximately equilibrium star clusters or spherical galaxy models.

Both generators write the same CSV output format:

```csv
id,mass,pos_x,pos_y,pos_z,vel_x,vel_y,vel_z
```

The config formats are intentionally different:

- `generate_data.py` uses `galaxies` with disk-specific fields such as `disk_radius`, `disk_thickness`, `disk_orientation`, and `mass.black_hole_factor`.
- `generate_plummer_data.py` uses `systems` with Plummer-specific fields such as `total_mass`, `scale_radius`, `cutoff_radius`, `correct_center_of_mass`, and `correct_momentum`.

## Config Format

Minimal example:

```json
{
  "seed": 12345,
  "output": "data.csv",
  "galaxies": [
    {
      "name": "main",
      "particles": 700,
      "disk_radius": 16.0,
      "disk_thickness": 1.5,
      "position": [0.0, 0.0, 0.0],
      "velocity": [0.0, 0.0, 0.0]
    },
    {
      "name": "satellite",
      "particles": 200,
      "disk_radius": 6.0,
      "disk_thickness": 0.8,
      "orbit": {
        "around": "main",
        "mode": "crash",
        "distance": 30.0,
        "angle": 45.0,
        "tangential_velocity_factor": 0.5,
        "radial_velocity_factor": 0.7
      },
      "disk_orientation": {
        "inclination": 65.0,
        "azimuth": 140.0
      }
    }
  ]
}
```

The `seed` field is required. The same config with the same seed produces the same CSV output.

There is no global particle count. Each galaxy defines its own `particles` count. The first particle of each galaxy is its central black hole, all remaining particles are disk particles.

## Galaxy Fields

Required for every galaxy:

- `name`: unique galaxy name
- `particles`: number of particles in this galaxy, including the central black hole
- `disk_radius`: maximum disk radius
- `disk_thickness`: vertical disk thickness

Root galaxies without `orbit` also require:

- `position`: absolute position `[x, y, z]`
- `velocity`: absolute velocity `[vx, vy, vz]`

Orbiting galaxies require:

- `orbit`: orbit definition

Optional:

- `disk_orientation`: orientation of the disk plane
- `mass`: per-galaxy mass settings

## Disk Orientation

`disk_orientation` controls how the galaxy disk is tilted in 3D space:

```json
"disk_orientation": {
  "inclination": 45.0,
  "azimuth": 120.0
}
```

- `inclination`: tilt away from the XY plane in degrees
- `azimuth`: direction of the tilt in degrees

If omitted, the disk lies in the XY plane:

```json
"disk_orientation": {
  "inclination": 0.0,
  "azimuth": 0.0
}
```

## Orbits

An orbiting galaxy specifies which galaxy it orbits via `orbit.around`:

```json
"orbit": {
  "around": "main",
  "mode": "stable",
  "distance": 40.0,
  "angle": 120.0
}
```

Supported orbit modes:

- `stable`: approximately tangential circular orbit around the host
- `crash`: reduced tangential velocity plus radial velocity toward the host
- `static`: placed relative to the host without additional orbit velocity
- `custom`: custom relative position and velocity

For `stable`, `crash`, and `static`:

- `distance`: distance to the host galaxy
- `angle`: position angle around the host in degrees in the XY plane
- `height`: optional Z offset relative to the host

For `stable` and `crash`:

- `tangential_velocity_factor`: multiplier for tangential orbit velocity
- `radial_velocity_factor`: multiplier for radial velocity toward the host

For `custom`:

```json
"orbit": {
  "around": "main",
  "mode": "custom",
  "position": [30.0, 0.0, 0.0],
  "velocity": [0.0, 4.0, 0.0]
}
```

The custom `position` and `velocity` are relative to the host galaxy.

## Mass Settings

The central black hole mass is calculated automatically by default:

```text
black_hole_mass = average_particle_mass * max(particles - 1, 1) * black_hole_factor
```

Default mass settings:

```json
"mass": {
  "particle_min": 0.03,
  "particle_max": 20.0,
  "black_hole_factor": 10.0
}
```

Defaults can be configured globally:

```json
"defaults": {
  "mass": {
    "particle_min": 0.03,
    "particle_max": 20.0,
    "black_hole_factor": 10.0
  }
}
```

They can also be overridden per galaxy:

```json
"mass": {
  "particle_min": 0.03,
  "particle_max": 20.0,
  "black_hole_factor": 12.0
}
```

For full control, set the black hole mass directly:

```json
"mass": {
  "black_hole_mass": 1000000.0
}
```

## Output Format

The resulting CSV file contains 8 columns:

```csv
id,mass,pos_x,pos_y,pos_z,vel_x,vel_y,vel_z
```

- `id`: unique ID of the particle
- `mass`: mass of the particle
- `pos_x`, `pos_y`, `pos_z`: position in 3D space
- `vel_x`, `vel_y`, `vel_z`: velocity in 3D space

## Examples

Stable multi-galaxy setup:

```shell
python3 generate_data.py --config examples/stable_three_galaxies.json
```

Collision setup:

```shell
python3 generate_data.py --config examples/collision_three_galaxies.json
```
