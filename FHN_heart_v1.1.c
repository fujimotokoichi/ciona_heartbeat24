#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main(void){
    FILE *fp;
    char filename[21];
    // Variables at time t: v1, w1; Work memory of variables at time t+1 v2, w2.
    // v1 and v2 correspond to the variable v; w1 and w2 correspond to variable w in Eqs. 1--2 in p. 12. 
    // Thiese are used for update with Runge-Kutta method.
    // The index of variables cooresponds to the "i" in p.12.
    double v1[21], v2[21], w1[21], w2[21];
    // Differences in variables used for 4th order Runge-Kutta method.
    // kp and jp(p=1--4) correspnd to v and w, respectively. 
    double k1[21], k2[21], k3[21], k4[21], j1[21], j2[21], j3[21], j4[21];
    // The patameters in Eqs. 1--2.
    double a[21], b[21], c[21];
    // The cell index for loop.
    int l;
    // t is time.
    // dt is time difference.
    // I value for debug.
    // C is defined by C dt^2 = D with D is the diffusion constant except for boundary i = 1 and 20.
    // C1 is defined by C1 dt^2 =D1 with D is the diffusion constact at the boundary i = 1 and 20.
    double t, dt, I, C, C1;
 
    fp=fopen( "FHN_output.txt", "w");  
    if(fp==NULL){
    printf("cannot open\n");
    exit(1);
    }

   
    for(int j=1; j<21; ++j){
        dt=0.0010;
        //dx=0.001;
        C=5.0;  // C*dt<0.3
        C1=5.0;
        
        //The setting for the parameters. 
        //The frequency increases with a and b decreasing or c increasing.
        // For pacemekers at boundary i = 1.
        a[1]=0.37;         
        b[1]=0.80;                            
        c[1]=10; 
        // For W cells i = 2--19.
        for(l=2; l<20; ++l){ 
            a[l]=0.44;         
            b[l]=0.80;
            c[l]=10;
            }
        // For pacemekers at boundary i = 20.
        a[20]=0.44;        
        b[20]=0.80;
        c[20]=10;
        // The initial conditions for v1 and w1.
        // For pacemekers at boundary i = 1.
        v1[1]=-5.0;
        w1[1]=0.0;
        // For W cells i = 2--19.
        for(l=2; l<20; ++l){ 
            v1[l]=-8.0;
            w1[l]=0.0;
        }
        // For pacemekers at boundary i = 20.
        v1[20]=-5.0;
        w1[20]=0.0;
        // Starting point of time iteration.
        int i;
        for(i=0; i<100000; ++i){ //
        
            t=i*dt;
            I=0.0;
            // output for debug.
            if(t>=10 && t<90){
            // b[1]=b[1]+0.0000050;
            // a[20]=a[20]-0.0000050;
            }
            if(t>=40 && t<60){
            // a[1]=a[1]+0.0000050;
            // a[20]=a[20]-0.0000050;
            }
            if(t>=40 && t<60){
            // a[20]=a[20]-0.0000050;
            }
        
            if(i%1000==0){
                printf("t=%lf;a1=%lf;a20=%lf;b1=%lf;b20=%lf\n", t, a[1], a[20],b[1], b[20]);
            }
            // 
            double s=1.0;
            // Changes at the pacemaker at i=1.
            k1[1]=dt*(c[1]*(-v1[1]*v1[1]*v1[1]/3.0+v1[1]-w1[1]+I)+C1*(v1[2]-v1[1]));
            j1[1]=dt*(v1[1]-b[1]*w1[1]+a[1]);

            k2[1]=dt*(c[1]*(-(v1[1]+(k1[1]/2.0))*(v1[1]+(k1[1]/2.0))*(v1[1]+(k1[1]/2.0))/3.0+(v1[1]+(k1[1]/2.0))-(w1[1]+(j1[1]/2.0))+I)+C1*((v1[2]+k1[1]/2.0)-(v1[1]+(k1[1]/2.0))));
            j2[1]=dt*((v1[1]+(k1[1]/2.0))-b[1]*(w1[1]+(j1[1]/2.0))+a[1]);

            k3[1]=dt*(c[1]*(-(v1[1]+(k2[1]/2.0))*(v1[1]+(k2[1]/2.0))*(v1[1]+(k2[1]/2.0))/3.0+(v1[1]+(k2[1]/2.0))-(w1[1]+(j2[1]/2.0))+I)+C1*((v1[2]+k2[1]/2.0)-(v1[1]+(k2[1]/2.0))));
            j3[1]=dt*((v1[1]+(k2[1]/2.0))-b[1]*(w1[1]+(j2[1]/2.0))+a[1]);

            k4[1]=dt*(c[1]*(-(v1[1]+k3[1])*(v1[1]+k3[1])*(v1[1]+k3[1])/3.0+(v1[1]+k3[1])-(w1[1]+j3[1])+I)+C1*((v1[2]+k3[1])-(v1[1]+k3[1])));
            j4[1]=dt*((v1[1]+k3[1])-b[1]*(w1[1]+j3[1])+a[1]);

            v2[1]=v1[1]+s*(k1[1]+2.0*k2[1]+2.0*k3[1]+k4[1])/6.0;
            w2[1]=w1[1]+s*(j1[1]+2.0*j2[1]+2.0*j3[1]+j4[1])/6.0;
            // Changes at the pacemaker at W cells for i=2--19.
            for(l=2; l<20; ++l){  
                k1[l]=dt*(c[l]*(-v1[l]*v1[l]*v1[l]/3.0+v1[l]-w1[l]+I)+C*(v1[l+1]+v1[l-1]-2.0*v1[l]));
                j1[l]=dt*(v1[l]-b[l]*w1[l]+a[l]);

                k2[l]=dt*(c[l]*(-(v1[l]+(k1[l]/2.0))*(v1[l]+(k1[l]/2.0))*(v1[l]+(k1[l]/2.0))/3.0+(v1[l]+(k1[l]/2.0))
                -(w1[l]+(j1[l]/2.0))+I)+C*((v1[l+1]+k1[l]/2.0)+(v1[l-1]+(k1[l]/2.0))-2.0*(v1[l]+(k1[l]/2.0))));
                j2[l]=dt*((v1[l]+(k1[l]/2.0))-b[l]*(w1[l]+(j1[l]/2.0))+a[l]);

                k3[l]=dt*(c[l]*(-(v1[l]+(k2[l]/2.0))*(v1[l]+(k2[l]/2.0))*(v1[l]+(k2[l]/2.0))/3.0+(v1[l]+(k2[l]/2.0))
                -(w1[l]+(j2[l]/2.0))+I)+C*((v1[l+1]+k2[l]/2.0)+(v1[l-1]+(k2[l]/2.0))-2.0*(v1[l]+(k2[l]/2.0))));
                j3[l]=dt*((v1[l]+(k2[l]/2.0))-b[l]*(w1[l]+(j2[l]/2.0))+a[l]);

                k4[l]=dt*(c[l]*(-(v1[l]+k3[l])*(v1[l]+k3[l])*(v1[l]+k3[l])/3.0+(v1[l]+k3[l])
                -(w1[l]+j3[l])+I)+C*((v1[l+1]+k3[l])+(v1[l-1]+k3[l])-2.0*(v1[l]+k3[l])));
                j4[l]=dt*((v1[l]+k3[l])-b[l]*(w1[l]+j3[l])+a[l]);

                v2[l]=v1[l]+(k1[l]+2.0*k2[l]+2.0*k3[l]+k4[l])/6.0;
                w2[l]=w1[l]+(j1[l]+2.0*j2[l]+2.0*j3[l]+j4[l])/6.0;
                    
            }
            // Changes at the pacemaker at i=20.
            double u=1.0;
            k1[20]=dt*(c[20]*(-v1[20]*v1[20]*v1[20]/3.0+v1[20]-w1[20]+I)+C1*(v1[19]-v1[20]));
            j1[20]=dt*(v1[20]-b[20]*w1[20]+a[20]);

            k2[20]=dt*(c[20]*(-(v1[20]+(k1[20]/2.0))*(v1[20]+(k1[20]/2.0))*(v1[20]+(k1[20]/2.0))/3.0+(v1[20]+(k1[20]/2.0))-(w1[20]+(j1[20]/2.0))+I)+C1*((v1[19]+k1[20]/2.0)-(v1[20]+(k1[20]/2.0))));
            j2[20]=dt*((v1[20]+(k1[20]/2.0))-b[20]*(w1[20]+(j1[20]/2.0))+a[20]);

            k3[20]=dt*(c[20]*(-(v1[20]+(k2[20]/2.0))*(v1[20]+(k2[20]/2.0))*(v1[20]+(k2[20]/2.0))/3.0+(v1[20]+(k2[20]/2.0))-(w1[20]+(j2[20]/2.0))+I)+C1*((v1[19]+k2[20]/2.0)-(v1[20]+(k2[20]/2.0))));
            j3[20]=dt*((v1[20]+(k2[20]/2.0))-b[20]*(w1[20]+(j2[20]/2.0))+a[20]);

            k4[20]=dt*(c[20]*(-(v1[20]+k3[20])*(v1[20]+k3[20])*(v1[20]+k3[20])/3.0+(v1[20]+k3[20])-(w1[20]+j3[20])+I)+C1*((v1[19]+k3[20])-(v1[20]+k3[20])));
            j4[20]=dt*((v1[20]+k3[20])-b[20]*(w1[20]+j3[20])+a[20]);

            v2[20]=v1[20]+u*(k1[20]+2.0*k2[20]+2.0*k3[20]+k4[20])/6.0;
            w2[20]=w1[20]+u*(j1[20]+2.0*j2[20]+2.0*j3[20]+j4[20])/6.0;

            if(i % 5 == 0){//To reduce HDD access cost in plotting data. 
                fprintf(fp, "%lf\t%d\t%lf\t%lf\n", t, j, v1[j], w1[j]); 
            }
 // update of v1 and w1 for the next iteration
            for(l=1; l<21; ++l){
                v1[l]=v2[l];
                w1[l]=w2[l];
            }

        }

    }
    fclose(fp);
    return 0;
}