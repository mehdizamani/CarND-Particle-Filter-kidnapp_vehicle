/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h> 
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>
#include "map.h"
#include "helper_functions.h"

using std::string;
using std::vector;

using namespace std;


default_random_engine gen;
void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   *  Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  if (initialized() == false)
  {
	  // set the initialization flag to true
	  is_initialized = true;

	  //  Set the number of particles 
	  num_particles = 100; 
	  
	  // Set the yaw thereshold to be used in the prediction method
	  // this parameter might be useful in tweeking the final result to be faster
	  yaw_zero_thereshold = 0.0001;
	  // Create the normal distribution for x , y and theta with 
	  // mean x and y and theta and with standard deviation 
	  // std[0],std[1],std[2] respectively
	  normal_distribution<double> nDist_x(x,std[0]);
	  normal_distribution<double> nDist_y(y,std[1]);
	  normal_distribution<double> nDist_theta(theta,std[2]);
	  // This loop initialize our particels one by one and push them
	  // back to the particles vector.
	  // Using the normal distribution random variable generator, 
	  // we create some randome variable x,y,theta and assign them
	  // to our particles one by one
	  for (int i=0 ;  i <  num_particles ; i++)
	  {
		  Particle temp_particle;
		  temp_particle.id = i;
		  temp_particle.x = nDist_x(gen);
		  temp_particle.y = nDist_y(gen);
		  temp_particle.theta = nDist_theta(gen);
		  temp_particle.weight = 1.0;
		  particles.push_back(temp_particle);
		  weights.push_back(temp_particle.weight);
	  }
  }
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   *  Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
   
    //create 3 set of distribution with mean zero and std_pos[0], std_pos[1], std_pos[2]
	normal_distribution<double> nDist_x(0,std_pos[0]);
	normal_distribution<double> nDist_y(0,std_pos[1]);
	normal_distribution<double> nDist_theta(std_pos[2]);
	
	// going through all particles and update them and at the end, adding noise to them
	for (int i = 0; i < num_particles ; i++)
	{
		Particle p;
		p = particles[i];
		if ( fabs(yaw_rate) > yaw_zero_thereshold ) 
		{
			p.x = p.x + velocity / yaw_rate * (sin(p.theta + yaw_rate * delta_t) - sin(p.theta)); 
			p.y = p.y + velocity / yaw_rate * (cos(p.theta ) - cos(p.theta + yaw_rate * delta_t)); 
			p.theta = p.theta  + yaw_rate * delta_t ; 
		} 
		else
		{
			p.x = p.x + velocity * delta_t * cos(p.theta); 
			p.y = p.y + velocity * delta_t * sin(p.theta); 
		}
		// add noise
		p.x += nDist_x(gen);
		p.y += nDist_y(gen);
		p.theta += nDist_theta(gen);
		// update the ith particles with new theta and x and y values
		particles[i] = p;
	}
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   *  Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */

}



int ParticleFilter::findClosestLandmarkIndex(LandmarkObs obs, Map map_landmark)
{
	// This function goes through all of the map landmark and find the id of closest landmakr to the observaion
	double dClosestDistanceToLandmark = 999999;
	double dTemp = 0.0;
	int iID = -1;
	for ( int i = 0 ; i < map_landmark.landmark_list.size() ; i++)
	{
		Map::single_landmark_s temp = map_landmark.landmark_list[i];
		dTemp = dist(obs.x , obs.y , temp.x_f , temp.y_f);
		if (dTemp < dClosestDistanceToLandmark) 
		{
			dClosestDistanceToLandmark = dTemp;
			iID = i;
		}			
	}
	return iID;
}


void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   *  Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
   
   for (int i = 0 ; i < num_particles ; i ++)
   {
	   Particle p = particles[i];
	   vector<int> associations;
	   vector<double> sense_x;
	   vector<double> sense_y;
	   // vector<LandmarkObs> transformedLandmakObs_list;
	   LandmarkObs landmarkObs;
	   LandmarkObs transformedLandmarkObs;
	   double dTempWeight = 1.0;
	   double dProb = 0.0;
	   int itemp = 0 ;
	   // go through all of the observation and transformed them to vehicle coordinate 
	   for (int j = 0 ; j < observations.size() ; j++ )
	   {
		   landmarkObs = observations[j];
		   transformedLandmarkObs.x = p.x + landmarkObs.x * cos(p.theta) - landmarkObs.y * sin(p.theta);
		   transformedLandmarkObs.y = p.y + landmarkObs.x * sin(p.theta) + landmarkObs.y * cos(p.theta);
		    //transformedLandmakObs_list.push_back(tempLandmarkObs);		
		   itemp = findClosestLandmarkIndex(transformedLandmarkObs,map_landmarks);
		   if (itemp > -1) 
		   {
			    transformedLandmarkObs.id = itemp;
				auto closestLandmark = map_landmarks.landmark_list[transformedLandmarkObs.id];
		   
				dProb = (1/ (2 * M_PI * std_landmark[0] * std_landmark[1])) * 
					exp(-0.5 * ( pow ( (transformedLandmarkObs.x - closestLandmark.x_f) ,2 ) / pow (std_landmark[0],2) 
					           + pow ( (transformedLandmarkObs.y - closestLandmark.y_f) ,2 ) / pow (std_landmark[1],2) )) ;
				dTempWeight = dTempWeight * dProb ;
		   }
	   }
	   
	   SetAssociations(p, associations , sense_x, sense_y , dTempWeight);
	   particles[i] = p ;
	   weights[i] = dTempWeight;
   }
}


void ParticleFilter::resample() {
  /**
   *  Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
    // a discrete distribution of the weights is created
	discrete_distribution<> dDist_weights(weights.begin(), weights.end());
	vector<Particle> temp_particles;
	for (int i = 0 ; i < num_particles ; i++)
	{
		Particle p = particles[dDist_weights(gen)];
		p.id = i;
		temp_particles.push_back(p);
	}
    particles.clear();
	particles = temp_particles;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, 
                                     const std::vector<double>& sense_y, const double &weight) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
  particle.weight = weight;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
