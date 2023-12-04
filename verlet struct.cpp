#include <stdio.h>
#include <stdlib.h>
#include <math.h>

float f_Kepler(float a, float b){  //Ecuacion de la aceleración 2.2
    float f, den;
    den = pow((a*a+b*b), 1.5);
    f = -a/den;
    return f;
}

float Energia(float x, float y, float vx, float vy){ //Ecuacion de la energia 2.3
    float E;
    E = 0.5*(vx*vx+vy*vy) - 1/(sqrt(x*x+y*y));
    return E;
}

struct datos {
    float x;
    float y;
    float vx;
    float vy;
    float H;
    float error;
};

void Verlet(int steps, float h, float e, struct datos d){

    //DATOS INICIALES KEPLER (2.11)
    d.x = 1-e;
    d.y = 0;
    d.vx = 0;
    d.vy = sqrt((1.0+e)/(1.0-e));
    d.H = -1/2;
    d.error = 0;

    FILE * fverlet;
    fverlet=fopen("verlet.txt", "wt");

    for(int i=0; i<steps;i++){

        fprintf(fverlet, "%f %f %f %f\n" , d.x, d.y, d.H, d.error);

        d.vx += h*0.5*f_Kepler(d.x, d.y);
        d.vy += h*0.5*f_Kepler(d.y, d.x);

        d.x += h*d.vx;
        d.y += h*d.vy;

        d.vx += 0.5*h*f_Kepler(d.x, d.y);
        d.vy += 0.5*h*f_Kepler(d.y, d.x);

        d.H = Energia(d.x, d.y, d.vx, d.vy);
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

    Verlet(steps, h, e, data);

    return 0;
}
