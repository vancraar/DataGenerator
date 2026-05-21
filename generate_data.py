#!/usr/bin/env python3
# -*- coding: utf-8 -*-

########################################################################################################################
# Authors: Marcel Breyer, Alexander Van Craen                                                                          #
# Copyright (C): 2024 Alexander Van Craen, Marcel Breyer, and Dirk Pflüger                                             #
# License: This file is released under the MIT license. See the LICENSE file in the project root for full information. #
########################################################################################################################

import argparse
import csv
from dataclasses import dataclass
import json
import math
import random


# calculate G using the correct units
gravitational_constant_si = 6.67440e-11  # m^3 / (kg * s^2)
solar_mass_in_kg = 1.988435e30           # kg
parsec_in_m = 3.08567758129e16           # m
year_in_s = 365.25 * 86400.0             # s
G = (solar_mass_in_kg * year_in_s * (gravitational_constant_si / parsec_in_m) *
     (1.0 / parsec_in_m) * (1.0 / parsec_in_m) * year_in_s)


# a class representing a single particle
@dataclass
class Body:
    id: int
    mass: float
    pos_x: float
    pos_y: float
    pos_z: float
    vel_x: float
    vel_y: float
    vel_z: float


@dataclass
class Galaxy:
    name: str
    particles: int
    disk_radius: float
    disk_thickness: float
    particle_min: float
    particle_max: float
    black_hole_mass: float
    position: tuple
    velocity: tuple
    basis_u: tuple
    basis_v: tuple
    basis_n: tuple


def parse_args():
    parser = argparse.ArgumentParser(
        prog="generate_data",
        description="Generate N-body initial conditions for one or more configurable galaxies."
    )

    parser.add_argument("-c", "--config",
                        help="the JSON scenario config file",
                        required=True)
    parser.add_argument("-o", "--output",
                        help="override the output filename from the config")

    return parser.parse_args()


def load_config(path):
    with open(path, "r") as config_file:
        return json.load(config_file)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def get_number(config, key, default=None):
    if key not in config:
        return default

    value = config[key]
    require(isinstance(value, (int, float)), f"'{key}' must be a number")
    return float(value)


def get_vector(config, key):
    require(key in config, f"missing required vector '{key}'")
    value = config[key]
    require(isinstance(value, list) and len(value) == 3, f"'{key}' must be a vector with three numbers")
    require(all(isinstance(component, (int, float)) for component in value), f"'{key}' must contain only numbers")
    return tuple(float(component) for component in value)


def validate_config(config, output_override):
    require(isinstance(config, dict), "config must be a JSON object")
    require("seed" in config, "missing required 'seed'")
    require(isinstance(config["seed"], int), "'seed' must be an integer")
    require(output_override is not None or "output" in config, "missing required 'output' or --output override")
    require("galaxies" in config, "missing required 'galaxies'")
    require(isinstance(config["galaxies"], list) and len(config["galaxies"]) > 0,
            "'galaxies' must be a non-empty list")

    names = set()
    for index, galaxy in enumerate(config["galaxies"]):
        prefix = f"galaxies[{index}]"
        require(isinstance(galaxy, dict), f"{prefix} must be an object")
        require("name" in galaxy and isinstance(galaxy["name"], str) and galaxy["name"],
                f"{prefix}.name must be a non-empty string")
        require(galaxy["name"] not in names, f"duplicate galaxy name '{galaxy['name']}'")
        names.add(galaxy["name"])

        require("particles" in galaxy and isinstance(galaxy["particles"], int) and galaxy["particles"] >= 1,
                f"{prefix}.particles must be an integer >= 1")
        require(get_number(galaxy, "disk_radius") is not None and galaxy["disk_radius"] > 0,
                f"{prefix}.disk_radius must be > 0")
        require(get_number(galaxy, "disk_thickness") is not None and galaxy["disk_thickness"] >= 0,
                f"{prefix}.disk_thickness must be >= 0")

    for galaxy in config["galaxies"]:
        name = galaxy["name"]
        orbit = galaxy.get("orbit")
        if orbit is None:
            get_vector(galaxy, "position")
            get_vector(galaxy, "velocity")
        else:
            require(isinstance(orbit, dict), f"galaxy '{name}' orbit must be an object")
            require("position" not in galaxy and "velocity" not in galaxy,
                    f"galaxy '{name}' must not define root position/velocity when orbit is set")
            require("around" in orbit and orbit["around"] in names,
                    f"galaxy '{name}' orbit.around must reference an existing galaxy")
            require(orbit["around"] != name, f"galaxy '{name}' cannot orbit around itself")
            require(orbit.get("mode") in ("stable", "crash", "static", "custom"),
                    f"galaxy '{name}' orbit.mode must be stable, crash, static, or custom")

            if orbit["mode"] == "custom":
                get_vector(orbit, "position")
                get_vector(orbit, "velocity")
            else:
                require(get_number(orbit, "distance") is not None and orbit["distance"] > 0,
                        f"galaxy '{name}' orbit.distance must be > 0")
                require(get_number(orbit, "angle") is not None,
                        f"galaxy '{name}' orbit.angle must be set")

    resolve_galaxy_order(config["galaxies"])


def resolve_galaxy_order(galaxy_configs):
    by_name = {galaxy["name"]: galaxy for galaxy in galaxy_configs}
    temporary = set()
    permanent = set()
    ordered = []

    def visit(name):
        require(name not in temporary, f"host orbit cycle detected at galaxy '{name}'")
        if name in permanent:
            return

        temporary.add(name)
        orbit = by_name[name].get("orbit")
        if orbit is not None:
            visit(orbit["around"])
        temporary.remove(name)
        permanent.add(name)
        ordered.append(by_name[name])

    for galaxy in galaxy_configs:
        visit(galaxy["name"])

    return ordered


def merge_mass_config(defaults, galaxy_config):
    default_mass = defaults.get("mass", {}) if isinstance(defaults.get("mass", {}), dict) else {}
    galaxy_mass = galaxy_config.get("mass", {}) if isinstance(galaxy_config.get("mass", {}), dict) else {}
    mass = {**default_mass, **galaxy_mass}

    particle_min = float(mass.get("particle_min", 0.03))
    particle_max = float(mass.get("particle_max", 20.0))
    black_hole_factor = float(mass.get("black_hole_factor", 10.0))

    require(particle_min > 0, f"galaxy '{galaxy_config['name']}' mass.particle_min must be > 0")
    require(particle_max >= particle_min, f"galaxy '{galaxy_config['name']}' mass.particle_max must be >= particle_min")

    if "black_hole_mass" in mass:
        black_hole_mass = float(mass["black_hole_mass"])
    else:
        avg_particle_mass = (particle_min + particle_max) / 2.0
        black_hole_mass = avg_particle_mass * max(galaxy_config["particles"] - 1, 1) * black_hole_factor

    require(black_hole_mass > 0, f"galaxy '{galaxy_config['name']}' black hole mass must be > 0")

    return particle_min, particle_max, black_hole_mass


def add(a, b):
    return a[0] + b[0], a[1] + b[1], a[2] + b[2]


def scale(v, factor):
    return v[0] * factor, v[1] * factor, v[2] * factor


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def norm(v):
    return math.sqrt(dot(v, v))


def normalize(v):
    length = norm(v)
    require(length > 0, "cannot normalize zero vector")
    return v[0] / length, v[1] / length, v[2] / length


def disk_basis(orientation):
    inclination = math.radians(float(orientation.get("inclination", 0.0)))
    azimuth = math.radians(float(orientation.get("azimuth", 0.0)))

    normal = normalize((math.sin(inclination) * math.cos(azimuth),
                        math.sin(inclination) * math.sin(azimuth),
                        math.cos(inclination)))

    helper = (0.0, 1.0, 0.0) if abs(normal[2]) > 0.999 else (0.0, 0.0, 1.0)
    u = normalize(cross(helper, normal))
    v = normalize(cross(normal, u))

    return u, v, normal


def transform(local, basis_u, basis_v, basis_n):
    return add(add(scale(basis_u, local[0]), scale(basis_v, local[1])), scale(basis_n, local[2]))


def orbit_state(orbit, host):
    if orbit["mode"] == "custom":
        return add(host.position, get_vector(orbit, "position")), add(host.velocity, get_vector(orbit, "velocity"))

    distance = float(orbit["distance"])
    angle = math.radians(float(orbit["angle"]))
    height = float(orbit.get("height", 0.0))

    relative_position = (distance * math.cos(angle), distance * math.sin(angle), height)
    position = add(host.position, relative_position)

    if orbit["mode"] == "static":
        return position, host.velocity

    radial_direction = normalize(relative_position)
    tangential_direction = (-math.sin(angle), math.cos(angle), 0.0)
    if norm(tangential_direction) == 0.0:
        tangential_direction = (0.0, 1.0, 0.0)
    tangential_direction = normalize(tangential_direction)

    orbital_speed = math.sqrt(G * host.black_hole_mass / norm(relative_position))

    if orbit["mode"] == "stable":
        tangential_factor = float(orbit.get("tangential_velocity_factor", 1.0))
        radial_factor = float(orbit.get("radial_velocity_factor", 0.0))
    else:
        tangential_factor = float(orbit.get("tangential_velocity_factor", 0.5))
        radial_factor = float(orbit.get("radial_velocity_factor", 0.7))

    relative_velocity = add(scale(tangential_direction, orbital_speed * tangential_factor),
                            scale(radial_direction, -orbital_speed * radial_factor))

    return position, add(host.velocity, relative_velocity)


def create_galaxies(config):
    defaults = config.get("defaults", {}) if isinstance(config.get("defaults", {}), dict) else {}
    ordered_configs = resolve_galaxy_order(config["galaxies"])
    galaxies = {}
    result = []

    for galaxy_config in ordered_configs:
        orientation = galaxy_config.get("disk_orientation", {})
        require(isinstance(orientation, dict), f"galaxy '{galaxy_config['name']}' disk_orientation must be an object")

        particle_min, particle_max, black_hole_mass = merge_mass_config(defaults, galaxy_config)
        basis_u, basis_v, basis_n = disk_basis(orientation)

        orbit = galaxy_config.get("orbit")
        if orbit is None:
            position = get_vector(galaxy_config, "position")
            velocity = get_vector(galaxy_config, "velocity")
        else:
            position, velocity = orbit_state(orbit, galaxies[orbit["around"]])

        galaxy = Galaxy(galaxy_config["name"], galaxy_config["particles"], float(galaxy_config["disk_radius"]),
                        float(galaxy_config["disk_thickness"]), particle_min, particle_max, black_hole_mass,
                        position, velocity, basis_u, basis_v, basis_n)
        galaxies[galaxy.name] = galaxy
        result.append(galaxy)

    return result


def generate_galaxy_particles(galaxy, rng, start_id):
    particles = [Body(start_id, galaxy.black_hole_mass,
                      galaxy.position[0], galaxy.position[1], galaxy.position[2],
                      galaxy.velocity[0], galaxy.velocity[1], galaxy.velocity[2])]

    for particle_id in range(start_id + 1, start_id + galaxy.particles):
        radius = 0.1 + max(galaxy.disk_radius - 0.1, 0.0) * math.sqrt(rng.random())
        angle = 2.0 * math.pi * rng.random()
        local_position = (radius * math.sin(angle),
                          radius * math.cos(angle),
                          rng.uniform(-galaxy.disk_thickness, galaxy.disk_thickness))

        gravity_distance = norm(local_position)
        orbital_speed = math.sqrt(G * galaxy.black_hole_mass / gravity_distance)
        local_velocity = ((local_position[1] / radius) * orbital_speed,
                          (-local_position[0] / radius) * orbital_speed,
                          0.0)

        relative_position = transform(local_position, galaxy.basis_u, galaxy.basis_v, galaxy.basis_n)
        relative_velocity = transform(local_velocity, galaxy.basis_u, galaxy.basis_v, galaxy.basis_n)
        position = add(galaxy.position, relative_position)
        velocity = add(galaxy.velocity, relative_velocity)

        particles.append(Body(particle_id, rng.uniform(galaxy.particle_min, galaxy.particle_max),
                              position[0], position[1], position[2],
                              velocity[0], velocity[1], velocity[2]))

    return particles


def generate_particles(config):
    rng = random.Random(config["seed"])
    particles = []

    for galaxy in create_galaxies(config):
        particles.extend(generate_galaxy_particles(galaxy, rng, len(particles)))

    return particles


def write_csv(path, particles):
    with open(path, "w+", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(("id", "mass", "pos_x", "pos_y", "pos_z", "vel_x", "vel_y", "vel_z"))
        for particle in particles:
            writer.writerow((particle.id, particle.mass,
                             particle.pos_x, particle.pos_y, particle.pos_z,
                             particle.vel_x, particle.vel_y, particle.vel_z))


def main():
    args = parse_args()

    try:
        config = load_config(args.config)
        validate_config(config, args.output)
        output = args.output if args.output is not None else config["output"]
        write_csv(output, generate_particles(config))
    except (OSError, json.JSONDecodeError, ValueError) as error:
        raise SystemExit(f"error: {error}")


if __name__ == "__main__":
    main()
