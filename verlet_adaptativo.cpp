#include <stdio.h>
#include <stdlib.h>
#include <math.h>

float potencial(float a, float b){
    float f;
    f = -1/(sqrt(a*a+b*b));
    return f;
}

float derivada_potencial(float a, float b){
    float f, den;
    den = pow((a*a+b*b), 1.5);
    f = a/den;
    return f;
}

float energia(float x, float y, float vx, float vy){ //Considero que la masa es 1
    float E;
    E = 0.5*(vx*vx+vy*vy) - 1/(sqrt(x*x+y*y));
    return E;
}

struct datos {
    float x;
    float y;
    float t;
    float vx;
    float vy;
    float p;
    float H;
    float error;
};

void VerletAdapt(int steps, float h, float e, struct datos d){

    //DATOS INICIALES KEPLER (2.11)
    d.x = 1-e;
    d.y = 0;
    d.t = 0;
    d.vx = 0;
    d.vy = sqrt((1.0+e)/(1.0-e));
    d.p = sqrt(d.vx*d.vx+d.vy*d.vy); //Considero que la masa es 1
    d.H = -1/2;
    d.error = 0;

    FILE * fverlet;
    fverlet=fopen("verlet.txt", "wt");

    for(int i=0; i<steps;i++){

        fprintf(fverlet, "%f %f %f %f\n" , d.x, d.y, d.H, d.error);

        d.x += (0.5*h*d.vx)/(0.5*(d.vx*d.vx+d.vy*d.vy)+d.p);
        d.y += (0.5*h*d.vy)/(0.5*(d.vx*d.vx+d.vy*d.vy)+d.p);
        d.t += (0.5*h)/(0.5*(d.vx*d.vx+d.vy*d.vy)+d.p); //NO SE USA PERO LO DABAN EN LAS ECUACIONES

        d.vx += (h/(potencial(d.x, d.y)))*(derivada_potencial(d.x, d.y));
        d.vy += (h/(potencial(d.y, d.x)))*(derivada_potencial(d.y, d.x));

        d.x += (0.5*h*d.vx)/(0.5*(d.vx*d.vx+d.vy*d.vy)+d.p);
        d.y += (0.5*h*d.vy)/(0.5*(d.vx*d.vx+d.vy*d.vy)+d.p);
        d.t += (0.5*h)/(0.5*(d.vx*d.vx+d.vy*d.vy)+d.p); //NO SE USA PERO LO DABAN EN LAS ECUACIONES

        d.H = energia(d.x, d.y, d.vx, d.vy);
        d.error = 1+2*d.H; //SIMPLIFICACIÓN DE LA EXPRESIÓN ORIGINAL... DE DONDE SALE?

    }

    fclose (fverlet);
}

int main(){
    int steps;
    float h;
    float e;

    steps = 4000;
    h = 0.05;
    e = 0.6;

    struct datos data;

    VerletAdapt(steps, h, e, data);

    return 0;
}
