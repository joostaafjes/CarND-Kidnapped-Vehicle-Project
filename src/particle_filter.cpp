/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h>
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include <cfloat>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  // Set the number of particles. Initialize all particles to first position (based on estimates of
  //   x, y, theta and their uncertainties from GPS) and all weights to 1.
  // Add random Gaussian noise to each particle.
  // NOTE: Consult particle_filter.h for more information about this method (and others in this file).
  default_random_engine gen;

  // This line creates a normal (Gaussian) distribution for x.
  normal_distribution<double> dist_x(x, std[0]);

  // Create normal distributions for y and theta.
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

  num_particles = 10; // how to determine??
  for (int index = 0; index < num_particles; ++index) {
    Particle p;
    p.id = index;
    p.x = dist_x(gen);
    p.y = dist_y(gen);
    p.theta = dist_theta(gen);
    p.weight = 1.0;

    particles.push_back(p);

    weights.push_back(1.0);

    cout << "Particle init" << index << " " << p.x << " " << p.y << " " << p.theta << endl;
  }

  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
  // Add measurements to each particle and add random Gaussian noise.
  // NOTE: When adding no
  //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
  //  http://www.cplusplus.com/reference/random/default_random_engine/
  default_random_engine gen;

  normal_distribution<double> dist_x(0.0, std_pos[0]);
  normal_distribution<double> dist_y(0.0, std_pos[1]);
  normal_distribution<double> dist_theta(0.0, std_pos[2]);

  for (int index = 0; index < num_particles; ++index) {
    double theta = particles[index].theta;
    if (fabs(yaw_rate) == 0) {
      // straight
      particles[index].x += cos(theta) * velocity * delta_t;
      particles[index].y += sin(theta) * velocity * delta_t;
    } else {
      // turn
      particles[index].x += velocity / yaw_rate * (sin(theta + yaw_rate * delta_t) - sin(theta));
      particles[index].y += velocity / yaw_rate * (-cos(theta + yaw_rate * delta_t) + cos(theta));
      particles[index].theta += yaw_rate * delta_t;
    }

    // add noise
    particles[index].x += dist_x(gen);
    particles[index].y += dist_y(gen);
    particles[index].theta += dist_theta(gen);

    std::cout << "x       " << "y       " << "theta" << std::endl;
    std::cout << particles[index].x << " " << particles[index].y << " " << particles[index].theta << std::endl;
  }
}

//void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs> &observations) {
void ParticleFilter::dataAssociation(const Map &map_landmarks, LandmarkObs &observation) {
  // Find the predicted measurement that is closest to each observed measurement and assign the
  //   observed measurement to this particular landmark.
  // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
  //   implement this method and use it as a helper during the updateWeights phase.

  int best_match_id = -1;
  double min_dist = DBL_MAX;
  for (auto landmark: map_landmarks.landmark_list) {
    double distance = sqrt(pow((landmark.x_f - observation.x), 2) + pow((landmark.y_f - observation.y), 2));
    if (distance < min_dist) {
      min_dist = distance;
      best_match_id = landmark.id_i;
    }
  }
  // assign best match to observation
  observation.id = best_match_id;
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
  // Update the weights of each particle using a mult-variate Gaussian distribution. You can read
  //   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
  // NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
  //   according to the MAP'S coordinate system. You will need to transform between the two systems.
  //   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
  //   The following is a good resource for the theory:
  //   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
  //   and the following is a good resource for the actual equation to implement (look at equation
  //   3.33
  //   http://planning.cs.uiuc.edu/node99.html

  for (auto &particle: particles) {
    // match observations relative to this particle with existing landmarks (in Map)
    // transformation needed?? vehicle to map
    // use method dataAssociation

    long double sum_mvgp = 1.0f;

    std::vector<LandmarkObs> observation_map_coord;

    std::vector<int> associations;
    std::vector<double> sense_x, sense_y;

    // map observations from car to map coordinate systems for this particle, new structure
    for (auto observation: observations) {
      // observation in vehicle coordinate system transform to map coordinate system
      double x_map = particle.x + cos(particle.theta) * observation.x - sin(particle.theta) * observation.y;
      double y_map = particle.y + sin(particle.theta) * observation.x + cos(particle.theta) * observation.y;

      LandmarkObs landmarkObs;
      landmarkObs.x = x_map;
      landmarkObs.y = y_map;

      // match each observation to a map landmark (see dataAssociation above) -> set the id in the observation
      dataAssociation(map_landmarks, landmarkObs);
      observation_map_coord.push_back(landmarkObs);

      // preprare association for visualisation
      associations.push_back(landmarkObs.id);
      sense_x.push_back(x_map);
      sense_y.push_back(y_map);
    }

    // associate particle with observations for visualisation
    particle = SetAssociations(particle, associations, sense_x, sense_y);

    for (auto observation: observation_map_coord) {
      // calculate weight using Multivariate Gaussian probability density
      long double mvgp = exp(
          -(
              pow(observation.x - map_landmarks.landmark_list[observation.id - 1].x_f, 2) / (2 * pow(std_landmark[0], 2))
                  + pow((observation.y - map_landmarks.landmark_list[observation.id - 1].y_f), 2)
                      / (2 * pow(std_landmark[1], 2))
          )) /
          (2 * M_PI * std_landmark[0] * std_landmark[1]);

      int dec = abs(log10(mvgp)) + 3;
      mvgp = round(pow(10, dec) * mvgp + 0.5) / pow(10, dec);

      sum_mvgp *= mvgp;
    }

    // update weight for this particle
    particle.weight = sum_mvgp;
    std::cout << particle.weight << std::endl;
  }

  std::cout << particles[0].weight << std::endl;
}

void ParticleFilter::resample() {
  // Resample particles with replacement with probability proportional to their weight.
  // NOTE: You may find std::discrete_distribution helpful here.
  //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

  // get total weight
  long double sum_weights = 0.0;
  for (auto particle: particles) {
    std::cout << particle.weight << std::endl;
    sum_weights += particle.weight;
  }

  // normalize
  std::vector<int> weights;
  for (auto particle: particles) {
    weights.push_back(particle.weight / sum_weights * 1000);
  }

  std::discrete_distribution<int> d(weights.begin(), weights.end());

  std::vector<Particle> updated_particles;

  std::random_device rd;
  std::mt19937 gen(rd());

  // how many times??
  for (int i = 0; i < particles.size(); i++) {
    updated_particles.push_back(particles[d(gen)]);
  }
  particles = updated_particles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int> &associations,
                                         const std::vector<double> &sense_x, const std::vector<double> &sense_y) {
  //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates

  particle.associations = associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;

  return particle;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  stringstream ss;
  copy(v.begin(), v.end(), ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length() - 1);  // get rid of the trailing space
  return s;
}
string ParticleFilter::getSenseX(Particle best) {
  vector<double> v = best.sense_x;
  stringstream ss;
  copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length() - 1);  // get rid of the trailing space
  return s;
}
string ParticleFilter::getSenseY(Particle best) {
  vector<double> v = best.sense_y;
  stringstream ss;
  copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length() - 1);  // get rid of the trailing space
  return s;
}
