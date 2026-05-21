#!/usr/bin/env python3
# -*- coding: utf-8 -*-

########################################################################################################################
# Authors: Alexander Van Craen                                                                          #
# Copyright (C): 2026 Alexander Van Craen                                            #
# License: This file is released under the MIT license. See the LICENSE file in the project root for full information. #
########################################################################################################################

import argparse
import csv
from dataclasses import dataclass
import json
import math
import random
from typing import Optional


# calculate G using the correct units
gravitational_constant_si = 6.67440e-11  # m^3 / (kg * s^2)
solar_mass_in_kg = 1.988435e30           # kg
parsec_in_m = 3.08567758129e16           # m
year_in_s = 365.25 * 86400.0             # s
G = (solar_mass_in_kg * year_in_s * (gravitational_constant_si / parsec_in_m) *
     (1.0 / parsec_in_m) * (1.0 / parsec_in_m) * year_in_s)


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
class PlummerSystem:
    name: str
    particles: int
    total_mass: float
    scale_radius: float
    cutoff_radius: Optional[float]
    position: tuple
    velocity: tuple
    correct_center_of_mass: bool
    correct_momentum: bool


def parse_args():
    parser = argparse.ArgumentParser(
        prog="generate_plummer_data",
        description="Generate N-body initial conditions for one or more configurable Plummer systems."
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


def get_bool(config, key, default=None):
    if key not in config:
        return default

    value = config[key]
    require(isinstance(value, bool), f"'{key}' must be a boolean")
    return value


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
    require("systems" in config, "missing required 'systems'")
    require(isinstance(config["systems"], list) and len(config["systems"]) > 0,
            "'systems' must be a non-empty list")

    names = set()
    for index, system in enumerate(config["systems"]):
        prefix = f"systems[{index}]"
        require(isinstance(system, dict), f"{prefix} must be an object")
        require("name" in system and isinstance(system["name"], str) and system["name"],
                f"{prefix}.name must be a non-empty string")
        require(system["name"] not in names, f"duplicate system name '{system['name']}'")
        names.add(system["name"])

        require("particles" in system and isinstance(system["particles"], int) and system["particles"] >= 1,
                f"{prefix}.particles must be an integer >= 1")
        require(get_number(system, "total_mass") is not None and get_number(system, "total_mass") > 0,
                f"{prefix}.total_mass must be > 0")
        require(get_number(system, "scale_radius") is not None and get_number(system, "scale_radius") > 0,
                f"{prefix}.scale_radius must be > 0")

        cutoff_radius = get_number(system, "cutoff_radius")
        if cutoff_radius is not None:
            scale_radius = get_number(system, "scale_radius")
            require(cutoff_radius > scale_radius, f"{prefix}.cutoff_radius must be > scale_radius")

        get_bool(system, "correct_center_of_mass", True)
        get_bool(system, "correct_momentum", True)

    for system in config["systems"]:
        name = system["name"]
        orbit = system.get("orbit")
        if orbit is None:
            get_vector(system, "position")
            get_vector(system, "velocity")
        else:
            require(isinstance(orbit, dict), f"system '{name}' orbit must be an object")
            require("position" not in system and "velocity" not in system,
                    f"system '{name}' must not define root position/velocity when orbit is set")
            require("around" in orbit and orbit["around"] in names,
                    f"system '{name}' orbit.around must reference an existing system")
            require(orbit["around"] != name, f"system '{name}' cannot orbit around itself")
            require(orbit.get("mode") in ("stable", "crash", "static", "custom"),
                    f"system '{name}' orbit.mode must be stable, crash, static, or custom")

            if orbit["mode"] == "custom":
                get_vector(orbit, "position")
                get_vector(orbit, "velocity")
            else:
                require(get_number(orbit, "distance") is not None and get_number(orbit, "distance") > 0,
                        f"system '{name}' orbit.distance must be > 0")
                require(get_number(orbit, "angle") is not None,
                        f"system '{name}' orbit.angle must be set")

    resolve_system_order(config["systems"])


def resolve_system_order(system_configs):
    by_name = {system["name"]: system for system in system_configs}
    temporary = set()
    permanent = set()
    ordered = []

    def visit(name):
        require(name not in temporary, f"host orbit cycle detected at system '{name}'")
        if name in permanent:
            return

        temporary.add(name)
        orbit = by_name[name].get("orbit")
        if orbit is not None:
            visit(orbit["around"])
        temporary.remove(name)
        permanent.add(name)
        ordered.append(by_name[name])

    for system in system_configs:
        visit(system["name"])

    return ordered


def add(a, b):
    return a[0] + b[0], a[1] + b[1], a[2] + b[2]


def subtract(a, b):
    return a[0] - b[0], a[1] - b[1], a[2] - b[2]


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


def random_unit_vector(rng):
    z = rng.uniform(-1.0, 1.0)
    angle = 2.0 * math.pi * rng.random()
    xy = math.sqrt(max(0.0, 1.0 - z * z))
    return xy * math.cos(angle), xy * math.sin(angle), z


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
    tangential_direction = cross((0.0, 0.0, 1.0), radial_direction)
    if norm(tangential_direction) == 0.0:
        tangential_direction = (0.0, 1.0, 0.0)
    tangential_direction = normalize(tangential_direction)

    orbital_speed = math.sqrt(G * host.total_mass / norm(relative_position))

    if orbit["mode"] == "stable":
        tangential_factor = float(orbit.get("tangential_velocity_factor", 1.0))
        radial_factor = float(orbit.get("radial_velocity_factor", 0.0))
    else:
        tangential_factor = float(orbit.get("tangential_velocity_factor", 0.5))
        radial_factor = float(orbit.get("radial_velocity_factor", 0.7))

    relative_velocity = add(scale(tangential_direction, orbital_speed * tangential_factor),
                            scale(radial_direction, -orbital_speed * radial_factor))

    return position, add(host.velocity, relative_velocity)


def create_systems(config):
    ordered_configs = resolve_system_order(config["systems"])
    systems = {}
    result = []

    for system_config in ordered_configs:
        orbit = system_config.get("orbit")
        if orbit is None:
            position = get_vector(system_config, "position")
            velocity = get_vector(system_config, "velocity")
        else:
            position, velocity = orbit_state(orbit, systems[orbit["around"]])

        system = PlummerSystem(system_config["name"],
                               system_config["particles"],
                               float(system_config["total_mass"]),
                               float(system_config["scale_radius"]),
                               get_number(system_config, "cutoff_radius"),
                               position,
                               velocity,
                               get_bool(system_config, "correct_center_of_mass", True),
                               get_bool(system_config, "correct_momentum", True))
        systems[system.name] = system
        result.append(system)

    return result


def sample_plummer_radius(scale_radius, cutoff_radius, rng):
    for _ in range(100000):
        u = rng.random()
        if u == 0.0:
            continue

        radius = scale_radius / math.sqrt(u ** (-2.0 / 3.0) - 1.0)
        if cutoff_radius is None or radius <= cutoff_radius:
            return radius

    raise ValueError("failed to sample Plummer radius below cutoff_radius")


def sample_plummer_speed(radius, total_mass, scale_radius, rng):
    while True:
        q = rng.random()
        acceptance = q * q * (1.0 - q * q) ** 3.5
        if 0.1 * rng.random() <= acceptance:
            break

    potential = G * total_mass / math.sqrt(radius * radius + scale_radius * scale_radius)
    escape_speed = math.sqrt(2.0 * potential)
    return q * escape_speed


def correct_system_offsets(local_particles, correct_center_of_mass, correct_momentum):
    if not local_particles:
        return local_particles

    total_mass = sum(particle.mass for particle in local_particles)

    if correct_center_of_mass:
        center = scale(tuple(sum(particle.mass * component
                                 for particle, component in zip(local_particles, components))
                             for components in zip(*(particle.position for particle in local_particles))),
                       1.0 / total_mass)
        for particle in local_particles:
            particle.position = subtract(particle.position, center)

    if correct_momentum:
        mean_velocity = scale(tuple(sum(particle.mass * component
                                        for particle, component in zip(local_particles, components))
                                    for components in zip(*(particle.velocity for particle in local_particles))),
                              1.0 / total_mass)
        for particle in local_particles:
            particle.velocity = subtract(particle.velocity, mean_velocity)

    return local_particles


@dataclass
class LocalParticle:
    mass: float
    position: tuple
    velocity: tuple


def generate_system_particles(system, rng, start_id):
    particle_mass = system.total_mass / system.particles
    local_particles = []

    for _ in range(system.particles):
        radius = sample_plummer_radius(system.scale_radius, system.cutoff_radius, rng)
        position = scale(random_unit_vector(rng), radius)

        speed = sample_plummer_speed(radius, system.total_mass, system.scale_radius, rng)
        velocity = scale(random_unit_vector(rng), speed)

        local_particles.append(LocalParticle(particle_mass, position, velocity))

    correct_system_offsets(local_particles, system.correct_center_of_mass, system.correct_momentum)

    particles = []
    for index, particle in enumerate(local_particles):
        position = add(system.position, particle.position)
        velocity = add(system.velocity, particle.velocity)
        particles.append(Body(start_id + index, particle.mass,
                              position[0], position[1], position[2],
                              velocity[0], velocity[1], velocity[2]))

    return particles


def generate_particles(config):
    rng = random.Random(config["seed"])
    particles = []

    for system in create_systems(config):
        particles.extend(generate_system_particles(system, rng, len(particles)))

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
