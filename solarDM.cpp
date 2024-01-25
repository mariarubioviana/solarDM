#include <iostream>
#include <math.h>
#include<stdio.h>
#include<stdlib.h>
#include <vector>

namespace physics
{
	const double G = 1;
	const double H0 = -0.5;
	const double phi1 = 0;
	const double e = 0.6;
};

struct OrbitPoint{
	double x = 0;
	double y = 0;
	double z = 0;

	double vx = 0;
	double vy = 0;
	double vz = 0;
};

// A static body
struct Body{
	double x = 0;
	double y = 0;
	double z = 0;

	double Mass = 0;
};

struct Vector{
	double X = 0;
	double Y = 0;
	double Z = 0;
};



inline double VelocitySquared( const OrbitPoint &point )
{
	return (point.vx*point.vx + point.vy*point.vy + point.vz*point.vz);
}

inline double DistanceToCenterSquared( const OrbitPoint &point )
{
	return (point.x*point.x + point.y*point.y + point.z*point.z);
}

inline double DistanceToBodySquared( const OrbitPoint &point, const Body &body )
{
	double dx = point.x - body.x;
	double dy = point.y - body.y;
	double dz = point.z - body.z;

	return (dx*dx + dy*dy + dz*dz);
}

Vector KeplerForce(const OrbitPoint &particle, const std::vector <Body> bodies )
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

//double Energy( const OrbitPoint &point ){ return 0.5*(VelocitySquared(point)) - 1/(sqrt(DistanceSquared(point))); }

std::vector <Body> InitializeTwoBodies()
{
	std::vector <Body> bodies;
	
	Body body1;
	body1.x = 0;
	body1.y = 0;
	body1.z = 0;
	body1.Mass = 1;

	bodies.push_back( body1 );

	Body body2;
	body1.x = 100;
	body1.y = 0;
	body1.z = 0;
	body1.Mass = 2;
	bodies.push_back( body2 );
}

void PrintBody( std::vector<Body> bodies, int n = 0 )
{
	if( n >= bodies.size() )
		std::cout << "Index larger than number of bodies!" << std::endl;

	std::cout << "Mass: " << bodies[n].Mass << " X: " << bodies[n].x << " Y: " << bodies[n].y << " Z: " << bodies[n].z << std::endl;
}

int main()
{
	std::vector <Body> bodies = InitializeTwoBodies();

	PrintBody( bodies, 0 );
	PrintBody( bodies, 1 );
}

/*
void EulerExplicito(int pasos, double h);
void ImpMidpoint(int pasos, double h);
void Verlet(int pasos, double h);
void SympcEuler(int pasos, double h);
void printEexplicito(double x[], double y[],double ener[], double err[], int pasos);
void printEsymp(double x[], double y[],double ener[],double err[], int pasos);
void printMidpoint(std::vector <double> x, std::vector<double> y, std::vector<double> en, std::vector<double> err);
void printVerlet(double x[], double y[],double ener[],double err[], int pasos);
double Calcerror(double x, double y); ///El valor de d hay que ponerlo a mano
void CalculoAnalitico(int pasos);

void RK4(double pasos, double h);
// Es mejor omitir, para tener claro los miembros que vienen de std::
//using namespace std;

bool doVerlet = true;
bool doMidPoint = false;
bool doSymplectic = false;
bool doRK4 =true;
int main()
{

    //definición de d y L0
    double d = 1.0 - e*e;
    double L0 = sqrt(d);

    printf("Parametros del sistema: H0 = %f  e = %f L0 = %f d = %f\n",H0,e, L0, d);
    puts("Pulsa Intro para continuar ... ");
    getchar();

    // crear array de los pasos temporales para distintos métodos

    double h[4];
    int steps[4];
    h[0]= 0.0005;
    steps[0]=400000 ;//euler explicito
    h[1]=0.01;
    steps[1]=40000; //symplectic Euler, implicit midpoint, Verlet


    ///SOLUCION ANALITICA (sale)
    //CalculoAnalitico(steps[1]);
    puts("Calculo analitico acabado ... ");

    ///EULER EXPLÍCITO (sale)
    //EulerExplicito(steps[0], h[0]);
    // puts("Metodo Euler explicito acabado ... ");
    //getchar();

    ///VERLET (sale)
    if( doVerlet )
    {
        Verlet(steps[1],h[1]);
        puts("Metodo de Verlet acabado ... ");
        getchar();
    }

    ///IMPLICIT MIDPOINT
    if( doMidPoint )
    {
        ImpMidpoint(steps[1], h[1]);
        puts("Metodo implicit midpoint acabado ... ");
        getchar();
    }

    ///SYMPLECTIC EULER (sale)
    if( doSymplectic )
    {
        SympcEuler(steps[1],h[1]);
        puts("Metodo symplectic Euler acabado ... ");
        getchar();
    }

    ///RK4
    if( doRK4 )
    {
    RK4(steps[1],h[1]);
        puts("Metodo RK4 acabado ... ");
        getchar();
    }


    return 0;
}

void RK4(double pasos, double h){
{

}

void RK4(double pasos, double h){

	std::vector <double> x, y, vx, vy;
	std::vector <double> energy, err;

	x.push_back(0.0);
	y.push_back(0.0);
	vx.push_back(0);
	vy.push_back(2);
	energy.push_back(-0.5);
	err.push_back(0.0);

	for(int i=0; i < pasos-1; i++) {

		/// Getting the last vector element
		double Xo = x.back();
		double Yo = y.back();

		double VXo = vx.back();
		double VYo = vy.back();

		//printf( "xo = %f,  vx0 = %f, y0 = %f, vy0 = %f \n", Xo, VXo, Yo, VYo);

		double k1_vx = f_Kepler(Xo, Yo);
		double k1_vy = f_Kepler(Yo, Xo);

		double k1_x = VXo;
		double k1_y = VYo;


		//printf( "VX k1 = %f, VY k1 = %f \n", k1_vx,k1_vy);


		double k2_vx = f_Kepler(Xo + 0.5 * h * k1_x, Yo + 0.5 * h * k1_y );
		double k2_vy = f_Kepler(Yo + 0.5 * h * k1_y, Xo + 0.5 * h * k1_x );

		double k2_x = VXo + h * 0.5 * k1_vx;
		double k2_y = VYo + h * 0.5 * k1_vy;


		//printf( "VX k2 = %f, VY k2 = %f \n", k2_vx,k2_vy);

		double k3_vx = f_Kepler(Xo + 0.5 * h * k2_x, Yo + 0.5 * h * k2_y );
		double k3_vy = f_Kepler(Yo + 0.5 * h * k2_y, Xo + 0.5 * h * k2_x );

		double k3_x = VXo + h * 0.5 * k2_vx;
		double k3_y = VYo + h * 0.5 * k2_vy;

		// printf( "VX k3 = %f, VY k3 = %f \n", k3_vx,k3_vy);

		double k4_vx = f_Kepler(Xo + h * k3_x, Yo + h * k3_y );
		double k4_vy = f_Kepler(Yo + h * k3_y, Xo + h * k3_x );

		double k4_x = VXo + h * k3_vx;
		double k4_y = VYo + h * k3_vy;

		//printf( "VX k4 = %f, VY k4 = %f \n", k4_vx,k4_vy);


		vx.push_back( VXo + h * (k1_vx + 2.0 * k2_vx + 2.0 * k3_vx + k4_vx) * 0.1666666 );
		vy.push_back( VYo + h * (k1_vy + 2.0 * k2_vy + 2.0 * k3_vy + k4_vy) * 0.1666666 );


		x.push_back( Xo +  h * (k1_x + 2.0 * k2_x + 2.0 * k3_x + k4_x) * 0.1666666 );
		y.push_back( Yo +  h * (k1_y + 2.0 * k2_y + 2.0 * k3_y + k4_y) * 0.1666666 );



		energy.push_back( Energia( x.back(), y.back(), vx.back(), vy.back() ) );
		err.push_back(-(-0.5 - energy[i])/(-0.5));

		//printf("Paso %d. X=%f  Y= %f Energia=%f  VEL VX = %f VY =%f\n", i, x[i], y[i], energy[i], vx[i],vy[i]);
		// getchar();
	}

	printf("X=%f Y =%f VX=%f VY=%f", x[pasos-1], y[pasos-1], vx[pasos-1],vy[pasos-1]);
	FILE *f5 = fopen("RK4_prueba.txt", "w");
	if(!f5) { printf("Error al abrir el archivo de texto."); exit(1); }

	for(unsigned int i=0; i < x.size(); i++)
		fprintf(f5,"%f   %f    %f   %f   %f\n", 0.01 * (double) i, x[i], y[i], energy[i], err[i]);
	fclose(f5);
}







void printEexplicito(double x[], double y[],double ener[],double err[], int pasos){
    double aux;
FILE *f1;
 f1=fopen("Euler_explicito.txt", "w");
	if(f1==NULL){printf("Error al abrir el archivo de texto.");exit(1);	}

    for(int i=0;i<pasos;i++){
            aux=i*0.05;

            fprintf(f1,"%f   %f    %f   %f   %f\n",aux,x[i],y[i],ener[i], err[i]);
    }
    fclose(f1);



}

void printEsymp(double x[], double y[],double ener[],double err[], int pasos){
    double aux;
FILE *f2;
 f2=fopen("Euler_symp_0001.txt", "w");
	if(f2==NULL){printf("Error al abrir el archivo de texto.");exit(1);	}

    for(int i=0;i<pasos;i++){
            aux=i*0.001;
            fprintf(f2,"%f   %f    %f   %f   %f\n",aux,x[i],y[i],ener[i], err[i]);

    }
    fclose(f2);


}




void printVerlet(double x[], double y[],double ener[],double err[], int pasos){
    double aux;
FILE *f4;
 f4=fopen("Verlet_prueba.txt", "w");
	if(f4==NULL){printf("Error al abrir el archivo de texto.");exit(1);	}

    for(int i=0;i<pasos;i++){
             aux=i*0.01;
            fprintf(f4,"%f   %f    %f   %f   %f\n",aux,x[i],y[i],ener[i], err[i]);
    }
    fclose(f4);
}


double Calcerror(double x, double y){
double phi,r,r_teo,err,d;
double pi;
pi=4*atan(1.0);
d=1-e*e;
if(x>0 && y>=0) phi=(double)atan(y/x);
if(x==0 && y>0) phi=0.5*pi;
if(x>0 && y <0) phi=(double)atan(y/x) + 2*pi;
if(x<0) phi=(double)atan(y/x) + pi;
if(x==0 && y<0) phi=1.5*pi;

r=sqrt(x*x+y*y);
r_teo= (d)/(1.0+e*(double)cos(phi));

err=(r_teo-r)*(r_teo - r);
return err;

}

void CalculoAnalitico(int pasos){
    double pi;
pi=4*atan(1.0);
double delta=(2*pi)/(double)pasos;
//printf("delta %f \n", delta);
double d=1-e*e;
double x[pasos];
double y[pasos];

double phi,r;
phi=0;

for(int i=0; i<pasos; i++){
r= (d)/(1.0+e*(double)cos(phi));
x[i]=r*cos(phi);
y[i]=r*sin(phi);
phi=phi+delta;
//printf("Paso %d. phi = %f, r = %f,  x = %f , y= %f \n", i,phi, r, x[i],y[i]);
//getchar();
}

FILE *f5;
 f5=fopen("Analitico.txt", "w");
	if(f5==NULL){printf("Error al abrir el archivo de texto.");exit(1);	}

    for(int i=0;i<pasos;i++){

            fprintf(f5,"%f   %f\n",x[i],y[i]);
    }
    fclose(f5);

}

*/

