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

void Verlet(int steps, float h, float e){
    float x[steps];
    float y[steps];
    float vx[steps];
    float vy[steps];
    float H[steps];
    float error[steps];

    for (int i=0; i<steps; i++){ //inicializamos
        x[i] = 0;
        y[i] = 0;
        H[i] = 0;
        error[i] = 0;
        vx[i] = 0;
        vy[i] = 0;
    }

    //DATOS INICIALES KEPLER (2.11)
    x[0] = 1-e;
    y[0] = 0;
    vx[0] = 0;
    vy[0] = sqrt((1.0+e)/(1.0-e));
    H[0] = -1/2;

    for(int i=0; i<(steps-1);i++){
        vx[i] += h*0.5*f_Kepler(x[i], y[i]);
        vy[i] += h*0.5*f_Kepler(y[i], x[i]);

        x[i+1] = x[i] + h*vx[i];
        y[i+1] = y[i] + h*vy[i];

        vx[i+1] = vx[i] + 0.5*h*f_Kepler(x[i+1], y[i+1]);
        vy[i+1] = vy[i] + 0.5*h*f_Kepler(y[i+1], x[i+1]);

        H[i+1] = Energia(x[i+1], y[i+1], vx[i+1], vy[i+1]);
        error[i+1] = 1+2*H[i]; //SIMPLIFICACIÓN DE LA EXPRESIÓN ORIGINAL... DE DONDE SALE?
    }

    FILE * fverlet;
        if((fverlet=fopen("verlet.txt", "wt"))==NULL)
            puts ("error al abrir el archivo");
        else{
            for (int i=0; i<steps; i++){
                fprintf(fverlet, "%f %f %f %f\n" , x[i], y[i], H[i], error[i]);
            }
            fclose(fverlet);
        }
}

int main(){
    int steps;
    float h;
    float e;

    steps = 4000;
    h = 0.05;
    e = 0.6;

    Verlet(steps, h, e);

    return 0;
}
