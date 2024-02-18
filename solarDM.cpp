#include <iostream>
#include <math.h>
#include<stdio.h>
#include<stdlib.h>
#include <vector>

bool debug = true;

namespace physics
{
	const double G = 1;
	const double H0 = -0.5;
	const double phi1 = 0;
	const double e = 0.6;
};

/// A structure to store each point of a given orbit
struct Particle{
	double x = 0;
	double y = 0;
	double z = 0;

	double vx = 0;
	double vy = 0;
	double vz = 0;
};

// A structure to store the coordinates and mass of a static body
struct Body{
	double x = 0;
	double y = 0;
	double z = 0;

	double Mass = 0;
};

struct NBodySystem
{
	std::vector <Body> bodies;

	Particle particle;
};

struct Vector{
	double X = 0;
	double Y = 0;
	double Z = 0;
};



inline double VelocitySquared( const Particle &point )
{
	return (point.vx*point.vx + point.vy*point.vy + point.vz*point.vz);
}

inline double DistanceToCenterSquared( const Particle &point )
{
	return (point.x*point.x + point.y*point.y + point.z*point.z);
}

Particle IncreasePosition( const Particle &point, double x, double y, double z )
{
	Particle newPoint = point;
	newPoint.x += x;
	newPoint.y += y;
	newPoint.z += z;

	return newPoint;
}

Particle IncreaseVelocity( const Particle &point, double vx, double vy, double vz )
{
	Particle newPoint = point;
	newPoint.vx += vx;
	newPoint.vy += vy;
	newPoint.vz += vz;

	return newPoint;
}

inline double DistanceToBodySquared( const Particle &point, const Body &body )
{
	double dx = point.x - body.x;
	double dy = point.y - body.y;
	double dz = point.z - body.z;

	return (dx*dx + dy*dy + dz*dz);
}

Vector KeplerForce(const Particle &particle, const std::vector <Body> &bodies )
{ 
	Vector force;
	std::cout << "X: " << force.X << " Y: " << force.Y << " Z: " << force.Z << std::endl;
	for( const auto &body: bodies )
	{
		double r3 = pow( DistanceToBodySquared( particle, body ), 1.5 );
		force.X -= physics::G * body.Mass * ( particle.x - body.x )/r3;
		force.Y -= physics::G * body.Mass * ( particle.y - body.y )/r3;
		force.Z -= physics::G * body.Mass * ( particle.z - body.z )/r3;
	}

	return force;
}

double PotentialEnergy(const Particle& particle, const Body& body) {
	// const double G = 6.67430e-11;  // Gravitational constant
	const double G = 1;
	double r = sqrt(DistanceToBodySquared(particle, body));
	return -G * body.Mass / r;
}

double KineticEnergy(const Particle& particle) {
	return 0.5 * (pow(particle.vx, 2) + pow(particle.vy, 2) + pow(particle.vz, 2));
}

double TotalEnergy( const Particle &particle, const std::vector <Body> &bodies ) { 
	double totalU = 0.0;  // Total potential energy
	for (const auto& body : bodies)
		totalU += PotentialEnergy(particle, body);

	double totalK = KineticEnergy(particle);
	return totalU + totalK;
}

NBodySystem InitAraujo_1( std::vector<double> params )
{
	NBodySystem nBodySystem;

	if( params.size() != 4 )
	{
		std::cout << "Error. Araujo. It requires 4 parameters (M1,M2,d,v)" << std::endl;
		return nBodySystem;
	}

	double M1 = params[0];
	double M2 = params[1];
	double d = params[2];
	double v = params[3];

	double mu1 = M1/(M1+M2);
	double mu2 = M2/(M1+M2);

	Body body;

	body.x = mu1;
	body.y = 0;
	body.z = 0;
	body.Mass = M1;

	std::vector <Body> bodies;
	bodies.push_back( body );

	body.x = -mu2;
	body.y = 0;
	body.z = 0;
	body.Mass = M2;

	bodies.push_back( body );

	// We assign the bodies
	nBodySystem.bodies = bodies;

	// We assign the particle
	Particle particle;
	particle.x = mu1 + d;
	particle.vy = v - d;
	nBodySystem.particle = particle;

}

NBodySystem InitializeNBodySystem(std::string config, std::vector <double> params) {
	
	NBodySystem nBodySystem;

	if( config == "Araujo_1" )
		return InitAraujo_1(params);

	return nBodySystem;
}

void PrintBody( const Body &body ) {
	std::cout << "Mass: " << body.Mass << " X: " << body.x << " Y: " << body.y << " Z: " << body.z << std::endl;
}

void PrintBodies( const std::vector<Body> &bodies ) {
	std::cout << " = Bodies = " << std::endl;
	std::cout << " ------------ " << std::endl;
	for( const auto &b : bodies )
		PrintBody( b );
	std::cout << " ------------ " << std::endl;
}

void PrintParticle( const Particle &particle ) {
	std::cout << " = Particle = " << std::endl;
	std::cout << " ------------ " << std::endl;
	std::cout << "x: " << particle.x << " y: " << particle.y << " z: " << particle.z << std::endl;
	std::cout << "vx: " << particle.vx << " vy: " << particle.vy << " vz: " << particle.vz << std::endl;
	std::cout << " ------------ " << std::endl;
}

void PrintNBodySystem( const NBodySystem &nBodySystem )
{
	PrintParticle( nBodySystem.particle );
	PrintBodies( nBodySystem.bodies );
}

Particle RK4( const Particle &particle, const std::vector <Body> &bodies, double h = 0.01){

	Particle newParticle = particle;
	Vector force = KeplerForce(particle, bodies);

 	/// ----->  k1 
	double k1_vx = force.X;
	double k1_vy = force.Y;
	double k1_vz = force.Z;

	double k1_x = particle.vx;
	double k1_y = particle.vy;
	double k1_z = particle.vz;

	newParticle = IncreasePosition( particle, 0.5*h*k1_x, 0.5*h*k1_y, 0.5*h*k1_z );
	force = KeplerForce( newParticle, bodies );

 	/// ----->  k2
	double k2_x = particle.vx + h * 0.5 * k1_vx;
	double k2_y = particle.vy + h * 0.5 * k1_vy;
	double k2_z = particle.vz + h * 0.5 * k1_vz;

	double k2_vx = force.X;
	double k2_vy = force.Y;
	double k2_vz = force.Z;

	newParticle = IncreasePosition( particle, 0.5*h*k2_x, 0.5*h*k2_y, 0.5*h*k2_z );
	force = KeplerForce( newParticle, bodies );

 	/// ----->  k3
	double k3_vx = force.X;
	double k3_vy = force.Y;
	double k3_vz = force.Z;

	double k3_x = particle.vx + h * 0.5 * k2_vx;
	double k3_y = particle.vy + h * 0.5 * k2_vy;
	double k3_z = particle.vz + h * 0.5 * k2_vz;

	newParticle = IncreasePosition( particle, h*k3_x, h*k3_y, h*k3_z );
	force = KeplerForce( newParticle, bodies );

 	/// ----->  k4
	double k4_vx = force.X;
	double k4_vy = force.Y;
	double k4_vz = force.Z;

	double k4_x = particle.vx + h * k3_vx;
	double k4_y = particle.vx + h * k3_vy;
	double k4_z = particle.vx + h * k3_vz;

	double dX = h * (k1_x + 2.0 * k2_x + 2.0 * k3_x + k4_x)/6.;
	double dY = h * (k1_y + 2.0 * k2_y + 2.0 * k3_y + k4_y)/6.;
	double dZ = h * (k1_z + 2.0 * k2_z + 2.0 * k3_z + k4_z)/6.;

	double dVX = h * (k1_vx + 2.0 * k2_vx + 2.0 * k3_vx + k4_vx)/6.;
	double dVY = h * (k1_vy + 2.0 * k2_vy + 2.0 * k3_vy + k4_vy)/6.;
	double dVZ = h * (k1_vz + 2.0 * k2_vz + 2.0 * k3_vz + k4_vz)/6.;

	newParticle = IncreasePosition( particle, dX, dY, dZ);
	newParticle = IncreaseVelocity( particle, dVX, dVY, dVZ);

	return newParticle;
}

int main() {
	/// M1, M2, d , v 
	std::vector <double> araujaParams {1., 1e-7, 0.00287, 0.0050};
	NBodySystem nBodySystem = InitializeNBodySystem( "Araujo_1", araujaParams );

	PrintNBodySystem( nBodySystem );

	/*
	Particle particle;

	particle = RK4 ( particle, bodies );
	*/
}

