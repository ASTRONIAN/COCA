Satellite Conjunction Assessment and Collision Avoidance
This Python script performs satellite conjunction assessment and collision avoidance using Two-Line Element (TLE) data. The script calculates the time of closest approach, distance of closest approach, probability of collision, and executes a collision avoidance maneuver if necessary.
Functions
The script includes the following main functions:
distance(sat1, sat2): Calculates the distance between two satellite positions.
time_of_closest_approach(sat1, sat2, t_span): Calculates the time of closest approach between two satellites within a specified time span.
collision_avoidance_maneuver(sat1, delta_v): Performs a collision avoidance maneuver by applying a delta-V to the satellite's velocity.
probability_of_collision(min_distance, R1, R2, sigma_radial, sigma_tangential): Calculates the probability of collision between two satellites within a specified time span.
calc_epoch(tle_file): Calculates the epoch from the TLE file.
Usage
Obtain orbital elements from TLE files for the satellites of interest.
Create Orbit objects for the satellites using the Poliastro library.
Define the time span for the conjunction assessment.
Calculate the time of closest approach and minimum distance between the satellites.
Calculate the probability of collision.
Define the minimum separation distance for the conjunction assessment.
If the probability of collision is higher than the acceptable limit, perform a collision avoidance maneuver.
Plot the orbits of both satellites before and after the maneuver.

Dependencies
numpy
scipy
astropy
poliastro
datetime
re
Example
The script demonstrates the conjunction assessment and collision avoidance for the International Space Station (ISS) and an Iridium satellite (IRI) using their TLE data. The script calculates the closest approach time, closest approach distance, and probability of collision. If the probability of collision is higher than the acceptable limit or the closest approach distance is less than the minimum separation distance, the script performs a collision avoidance maneuver and updates the closest approach time, distance, and probability of collision.