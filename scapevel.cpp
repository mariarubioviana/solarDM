#include <iostream>
#include <math.h>
#include<stdio.h>
#include<stdlib.h>
#include <fstream> 
#include <vector>
#include <string>
#include <sstream>

bool debug = false;
bool info = true;
double DT = 0.001; // differential step
std::vector<double> times;
std::vector<double> totalenergy;

namespace physics
{
	 double G = 1;
	 double H0 = -0.5;
	 double phi1 = 0;
	 double e = 0.6;
};

//// Vector3D class definition and operations
class Vector3D {
	public:
	double x = 0,y = 0,z = 0;
	Vector3D(double a = 0, double b = 0, double c = 0)  : x(a), y(b), z(c) { }
	double  X() const { return x; }
	double  Y() const { return y; }
	double  Z() const { return z; }

	void Print( std::string name = "" ) 
	{
		std::cout << "Vector3D=" << name << "("<<X()<<","<<Y()<<","<<Z()<<")" << std::endl;
	}

	/// Operators
	Vector3D& operator += (const Vector3D & b)  {
		x += b.X(); y += b.Y(); z += b.Z();
		return *this;
	}

	Vector3D& operator -= (const Vector3D & b)  {
		x -= b.X(); y -= b.Y(); z -= b.Z();
		return *this;
	}

	Vector3D& operator *= (const double & b)  {
		x = b*X(); y = b*Y(); z = b*Z();
		return *this;
	}
	
	Vector3D& operator /= (const double & b)  {
		x = b/X(); y = b/Y(); z = b/Z();
		return *this;
	}
};

Vector3D operator + (const Vector3D & a, const Vector3D & b)  {
	    return Vector3D(a) += b;
}

Vector3D operator - (const Vector3D & a, const Vector3D & b)  {
	    return Vector3D(a) -= b;
}

Vector3D operator * ( const double& a, const Vector3D & b)  {
	    return Vector3D( a*b.X(), a*b.Y(), a*b.Z());
}

Vector3D operator * (const Vector3D & b, const double& a )  {
	    return Vector3D( a*b.X(), a*b.Y(), a*b.Z());
}

Vector3D operator / ( const Vector3D & b, const double& a ) {
	    return Vector3D( b.X()/a, b.Y()/a, b.Z()/a);
}

// Body structure to store the coordinates and mass
struct Body{
	Vector3D position;
	Vector3D velocity;

	double M = 0;
	double R = 0;

	double Mass( double r3 = 0 ) const
	{
		if( r3 < R )
		{
			return M * r3/R/R/R;
		}
		else
			return M;
	}

	void Print( std::string name = "" ) const {
		std::cout << " = Body (" << name << ") = " << std::endl;
		std::cout << " ------------------- " << std::endl;
		std::cout << "Mass: " << M << " Radius: " << R << " X: " << position.X() << " Y: " << position.Y() << " Z: " << position.Z() << std::endl;
		std::cout << "vX: " << velocity.X() << " vY: " << velocity.Y() << " vZ: " << velocity.Z() << std::endl;
	}
};

/// Particle structure to store each point of a given orbit
struct Particle{
	Vector3D position;
	Vector3D velocity;

	/// Returns the squared velocity of the particle
	inline double VelocitySquared( )
	{
		return (velocity.X()*velocity.X() + velocity.Y()*velocity.Y() + velocity.Z()*velocity.Z());
	}

	 /// Returns the radial velocity of the particle
	double RadialVelocity( ) const
	{
		return (sqrt(velocity.X()*velocity.X() + velocity.Y()*velocity.Y()));
	}

        /// Returns the radial position of the particle
	double RadialPosition( ) const
	{
		return (sqrt(position.X()*position.X() + position.Y()*position.Y()));
	}

	/// Returns the distance squared between the particle and the origin
	inline double DistanceToCenterSquared( )
	{
		return (position.X()*position.X() + position.Y()*position.Y() + position.Z()*position.Z());
	}
	
	/// It returns the distance between the particle and the body given by argument
	inline double DistanceToBodySquared( const Body &b ) const {
		double dx = position.X() - b.position.X();
		double dy = position.Y() - b.position.Y();
		double dz = position.Z() - b.position.Z();

		return (dx*dx + dy*dy + dz*dz);
	}

	/// It returns the potential energy of the particle respect to a given body
	double PotentialEnergy( const Body &body) const {
		double r = sqrt(DistanceToBodySquared( body));
		return -physics::G * body.Mass(r*r*r) / r;
	}

	/// It returns the kinetic energy of the particle respect to a given body
	double KineticEnergy( const Body &body) const {
		return 0.5 * (pow(velocity.X()-body.velocity.X(), 2) + pow(velocity.Y()-body.velocity.Y(), 2) + pow(velocity.Z()-body.velocity.Z(), 2));
	}
	
	/// It returns the total energy, kinetic energy is calculated respect to a given body (given by n)
	double TotalEnergy(  const std::vector <Body> &bodies, int n = 0 ) { 
		if( bodies.size() == 0 ) return 0;

		double totalU = 0.0;  // Total potential energy
		for ( auto& body : bodies)
			totalU += PotentialEnergy( body);

		double totalK = KineticEnergy( bodies[n]);
		return totalU + totalK;
	}

	double BodyEnergy( const Body &body ) const { 
		double totalU = PotentialEnergy( body);  // Total potential energy
		double totalK = KineticEnergy( body);
		return totalU + totalK;
	}

	void Print( std::string name = "" ) {
		std::cout << " = Particle (" << name << ") = " << std::endl;
		std::cout << " ----------- " << std::endl;
		std::cout << "x: " << position.X() << " y: " << position.Y() << " z: " << position.Z() << std::endl;
		std::cout << "vx: " << velocity.X() << " vy: " << velocity.Y() << " vz: " << velocity.Z() << std::endl;
		std::cout << " ----------- " << std::endl;
	}
};

// A structure to store the complete gravitational system
struct NBodySystem {
	std::vector <Body> bodies;

	Particle particle;
};


///// Helper Methods /////
void GetChar() {
	std::cout << "Press a key to continue ... " << std::endl;
	getchar();
}

/// It returns the distance between the particle and the body given by argument
Vector3D KeplerForce( Particle p,  std::vector <Body> bodies ) {
	Vector3D force;
	if( debug ) force.Print("Init KeplerForce");
	int cont = 0;
	for(  auto &b: bodies )
	{
		cont++;
		double r3 = pow( p.DistanceToBodySquared( b ), 1.5 );
		double forceX = - physics::G * b.Mass(r3) * ( p.position.X() - b.position.X() )/r3;
		//std::cout << "Body: " << cont << " r3 : " << r3 << " Mass : " << b.Mass << " diff: " << p.position.X() - b.position.X() << " r3: " << r3 << " ForceX: " << forceX << std::endl;
		if( debug )
		{
			std::cout << "Body: " << cont << " Mass : " << b.Mass() << " diff: " << p.position.X() - b.position.X() << " r3: " << r3 << " ForceX: " << forceX << std::endl;
		}
		force.x -= physics::G * b.Mass(r3) * ( p.position.X() - b.position.X() )/r3;
		force.y -= physics::G * b.Mass(r3) * ( p.position.Y() - b.position.Y() )/r3;
		force.z -= physics::G * b.Mass(r3) * ( p.position.Z() - b.position.Z() )/r3;
		if( debug) force.Print( "Force body " + std::to_string(cont));

	}

	//GetChar();
	return force;
}

///// Initialization Methods /////
//////////////////////////////////
std::pair<Particle,std::vector<Body>> InitInes() {
	Particle particle;
	std::vector <Body> bodies;

	Body b;

	b.position.x = 0;
	b.position.y = 0;
	b.position.z = 0;

	/// Velocity at the synodic system
	b.velocity.x = 0;
	b.velocity.y = 0;
	b.velocity.z = 0;
	b.M = 1;

	bodies.push_back( b );

	// We assign the particle
	particle.position.x = 0.4;
	
	/// Velocity at the synodic system
	particle.velocity.y = 2;

	return {particle,bodies};
}

std::pair<Particle,std::vector<Body>> InitKKAxions ( std::vector<double> params ){
	Particle particle;
	std::vector <Body> bodies;
	
	std::cout << "InitKKAxions: ";
	for (const auto &p:params)
		std::cout << p << " ";
	std::cout << std::endl;

	if( params.size() != 3 )
	{
		std::cout << "Error. KKAxions. It requires 3 parameters (v,x,y)" << std::endl;
		return {particle,bodies};
	}

	double v = params[0];
	double x = params[1];
	double y = params[2];
	
	Body b;

	b.position.x = 0;
	b.position.y = 0;
	b.position.z = 0;

	/// Velocity at the synodic system
	b.velocity.x = 0;
	b.velocity.y = 0;
	b.velocity.z = 0;
	b.M = 1;
	b.R = 1; // Size of the Sun in solar radii

	bodies.push_back( b );

	// We assign the particle (inside the Sun)
	particle.position.x = x; //0.4
	particle.position.y = y; //-0.3
	
	particle.velocity.x = v; //1.2
	particle.velocity.y = 0.0; //0.1

	return {particle,bodies};
}

std::pair<Particle,std::vector<Body>> InitAraujo( std::vector<double> params ) {
	Particle particle;
	std::vector <Body> bodies;

	if( params.size() != 3 )
	{
		std::cout << "Error. Araujo. It requires 3 parameters (M2,v,d)" << std::endl;
		return {particle,bodies};
	}

	double M1 = 1 - params[0];
	double M2 = params[0];
	double v = params[1];
	double d = params[2];

	double mu1 = M1/(M1+M2);
	double mu2 = M2/(M1+M2);

	Body b;

	b.position.x = -mu2;
	b.position.y = 0;
	b.position.z = 0;

	/// Velocity at the synodic system
	b.velocity.x = 0;
	b.velocity.y = -mu2;
	b.velocity.z = 0;
	b.M = M1;

	bodies.push_back( b );

	b.position.x = mu1;
	b.position.y = 0;
	b.position.z = 0;
	
	/// Velocity at the synodic system
	b.velocity.x = 0;
	b.velocity.y = 0;
	b.velocity.z = 0;
	b.M = M2;

	bodies.push_back( b );

	// We assign the particle
	particle.position.x = mu1 + d;
	
	/// Velocity at the synodic system (relative to body 2)
	particle.velocity.y = v - d;

	return {particle, bodies};
}

std::pair<Particle,std::vector<Body>> InitAraujoModified( std::vector<double> params ) {
	Particle particle;
	std::vector <Body> bodies;
	
	if( params.size() != 4 )
	{
		std::cout << "Error. Araujo. It requires 4 parameters (M2,v,d,CentralMass)" << std::endl;
		return {particle,bodies};
	}

	double M1 = 1 - params[0];
	double M2 = params[0];
	double v = params[1];
	double d = params[2];
	double CentralMass = params[3];

	double mu1 = M1/(M1+M2);
	double mu2 = M2/(M1+M2);

	Body b;

	b.position.x = -mu2;
	b.position.y = 0;
	b.position.z = 0;

	/// Velocity at the synodic system
	b.velocity.x = 0;
	b.velocity.y = -mu2;
	b.velocity.z = 0;
	b.M = CentralMass;

	bodies.push_back( b );

	b.position.x = mu1;
	b.position.y = 0;
	b.position.z = 0;
	
	/// Velocity at the synodic system
	b.velocity.x = 0;
	b.velocity.y = 0;
	b.velocity.z = 0;
	b.M = M2;

	bodies.push_back( b );

	// We assign the particle
	particle.position.x = mu1 + d;
	
	/// Velocity at the synodic system (relative to body 2)
	particle.velocity.y = v - d;

	return {particle, bodies};
}

std::pair<Particle,std::vector<Body>> InitializeSystem(std::string config, std::vector <double> params = {}) {
	
	std::pair<Particle, std::vector<Body>> system = InitInes();

	if( config == "Araujo" )
		system = InitAraujo(params);

	if( config == "AraujoModified" )
		system = InitAraujoModified(params);

	if( config == "Ines" )
		system = InitInes();

	if( config == "KKaxions" )
		system = InitKKAxions(params);

	if( info )
	{
		system.first.Print("Init");
		for( const auto &b : system.second )
			b.Print();
	}

	return system;
}

///// Integration Methods /////
///////////////////////////////
void RK4( Particle &particle,  std::vector <Body> &bodies, double h = 0.01){

	// Initial values
	Particle initialParticle = particle;

	// k1, l1
	Vector3D k1 = h * initialParticle.velocity;
	Vector3D l1 = h * KeplerForce(particle, bodies);

	if( debug )
	{
		k1.Print("k1");
		l1.Print("l1");
	}

	// k2
	Particle k2particle;
	if( debug )
	{
		k2particle.Print("k2");
		KeplerForce( particle, bodies).Print("ForceInitial");
		initialParticle.Print("initial");
	}

	k2particle.position = initialParticle.position + k1/2.;
	k2particle.velocity = initialParticle.velocity + l1/2.;
	if( debug )
	{
		k2particle.Print("k2particle");
	}

	Vector3D k2 = h * k2particle.velocity;
	Vector3D l2 = h * KeplerForce( k2particle, bodies);
	if( debug )
	{
		k2.Print("k2");
		l2.Print("l2");
		KeplerForce( k2particle, bodies).Print("ForceK2");
	}

	// k3
	Particle k3particle;
	k3particle.position = initialParticle.position + k2/2.;
	k3particle.velocity = initialParticle.velocity + l2/2.;
	if( debug )
	{
		k3particle.Print("k3");
	}

	Vector3D k3 = h * k3particle.velocity;
	Vector3D l3 = h * KeplerForce( k3particle, bodies);
	if( debug )
	{
		KeplerForce( k3particle, bodies).Print("ForceK3");
		k3.Print("k3");
		l3.Print("l3");
	}

	// k3
	Particle k4particle;
	k4particle.position = initialParticle.position + k3;
	k4particle.velocity = initialParticle.velocity + l3;
	if( debug )
	{
		k4particle.Print("k4");
	}

	Vector3D k4 = h * k4particle.velocity;
	Vector3D l4 = h * KeplerForce( k4particle, bodies);
	if( debug )
	{
		k4.Print("k4");
		l4.Print("l4");
		KeplerForce( k4particle, bodies).Print("ForceK4");
	}

	// Update position and velocity
	particle.position = initialParticle.position + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
	particle.velocity = initialParticle.velocity + (l1 + 2.0 * l2 + 2.0 * l3 + l4) / 6.0;
}

///// Output method /////
/////////////////////////
//void WriteParticle( std::ofstream &ofile, double t, const Particle &p, const std::vector <Body> bodies ) {
	//ofile << t << "\t" <<  p.position.X() << "\t" << p.position.Y() << "\t" << p.position.Z() << "\t" << p.velocity.X() << "\t" << p.velocity.Y() << "\t" << p.velocity.Z();
	//times.push_back(t);
	//for( const auto &b: bodies)
		//ofile << "\t" << p.BodyEnergy(b); //NO SE COMO PONER AQUI LA ENERGÍA TOTAL
	//bodyenergy.push_back(p.BodyEnergy(bodies[1]));

	//ofile << "\n";
//}

///// Display help method /////
///////////////////////////////
void displayHelp(const char* programName) {
	std::cout << "Usage: " << programName << " [options]" << std::endl;
	std::cout << "Options:" << std::endl;
	std::cout << "  --output <output_filename>: Set the output filename (default: out.txt)" << std::endl;
	std::cout << "  --system <system_name>: Set the system name (default: Araujo)" << std::endl;
	std::cout << "  --params <param1,param2,...>: Set the parameters as a comma-separated list of values" << std::endl;
	std::cout << "  --steps <N>: Set the number of iterations (Default 2000)" << std::endl;
	std::cout << "  --help: Display this help message" << std::endl;
}

double GetEscapeVelocity( std::pair<Particle, std::vector<Body>> nBodySystem, std::string outfname, int steps )
{
	Particle particle = nBodySystem.first;
	std::vector <Body> bodies = nBodySystem.second;

	std::ofstream outputFile( outfname, std::ofstream::out ); // Open file for appending

	double deltaT = DT;
	double t = 0;
	//WriteParticle( outputFile, t, particle, bodies );
	
	int cont = 0;
	while ( cont < steps )
	{
		RK4 ( particle, bodies, deltaT );
		/* if( debug ) { PrintNBodySystem( nBodySystem ); } */
		//WriteParticle( outputFile, t, particle, bodies );
		if ( particle.RadialPosition() > 1){
			double en = particle.TotalEnergy(bodies);
			if (en > 0){
				return 1;
			}
		}
		t += deltaT;
		cont++;
	}
	outputFile.close();
	return 0;
}

//// It stars the main program ////
int main(int argc, char* argv[]) {

	std::string outputFilename = "out.txt";
	std::string systemName = "Araujo";
	// M2, v, d
	std::vector<double> params = { 1.e-5, 20, 0.0001};
	int steps = (int) (2./DT);
	double data;
	double vini = 0.8; //no sé que valores dar
	double vfin = 1.8; //no sé que valores dar

	// Parse command-line arguments
	for (int i = 1; i < argc; ++i) {
		std::string arg = argv[i];
		if (arg == "--output") {
			if (i + 1 < argc) {
				outputFilename = argv[i + 1];
				++i;  // skip the next argument
			} else {
				std::cerr << "--output option requires one argument." << std::endl;
				return 1;
			}
		} else if (arg == "--system") {
			if (i + 1 < argc) {
				systemName = argv[i + 1];
				++i;  // skip the next argument
			} else {
				std::cerr << "--system option requires one argument." << std::endl;
				return 1;
			}
		} else if (arg == "--steps") {
			if (i + 1 < argc) {
				steps = std::stoi(argv[i + 1]);
				++i;  // skip the next argument
			} else {
				std::cerr << "--system option requires one argument." << std::endl;
				return 1;
			}
		} else if (arg == "--params") {
			if (i + 1 < argc) {
				std::string paramsString = argv[i + 1];
				std::istringstream iss(paramsString);
				double param;
				char delimiter;
				params.clear();
				while (iss >> param) {
					params.push_back(param);
					if (!(iss >> delimiter && delimiter == ',')) {
						break;  // break if delimiter is not ','
					}
				}
				++i;  // skip the next argument
			} else {
				std::cerr << "--params option requires one argument." << std::endl;
				return 1;
			}
		}
		else if (arg == "--help") {
			displayHelp(argv[0]);
			return 0;
		}
	}

	// Now you can use the parsed arguments
	std::cout << "Output Filename: " << outputFilename << std::endl;
	std::cout << "System Name: " << systemName << std::endl;
	std::cout << "Params: ";

	if (params.empty()) {
		std::cout << "{}";
	} else {
		for (double param : params) {
			std::cout << param << " ";
		}
	}
	std::cout << std::endl;

	if (argc < 2) {
		displayHelp(argv[0]);
		return 1; // indicating an error
	}

	//std::pair<Particle, std::vector<Body>> nBodySystem = InitializeSystem(systemName, params );
	//WriteOrbit( nBodySystem, outputFilename, steps );
	std::ofstream archivo("escapevelocity.txt");

	for (double v = vini; v < vfin; v += 0.05){ //no sé que salto poner
		std::cout << "loop: ";
		for (const auto &p:params)
			std::cout << p << " ";
		std::cout << std::endl;
		std::cout << "v: " << v << std::endl;
		cout << "hola" << "\n";

		std::pair<Particle, std::vector<Body>> nBodySystem = InitializeSystem(systemName, {v, params[1], params[2]} );
		data = GetEscapeVelocity( nBodySystem, outputFilename, steps );
		if (data == 1){
			if (archivo.is_open()) {
        			archivo << v << "\n";
    			} 
			else {
       				std::cout << "No se pudo abrir el archivo." << std::endl;
    			}
			break; //cuando se encuentra la velocidad de escape salgo del bucle
		}
	}
	archivo.close();
}
