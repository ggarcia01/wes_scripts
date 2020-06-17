#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#define DT 0.001
#define N 10
#define MAX_RNDM 2147483647


/*

                          u1 - returns the potential energy of a particle interacting with the trap walls

*/

double u1(double x, double y){

  //variables
  double mu,r,u1_2d;
  mu = 10000;
  r = sqrt(x*x + y*y);
  u1_2d = mu * r*r / (1.0+r*r);
  return  u1_2d;
}

/*

                          new U2 - more realistic atomic interaction
                                 - repulsive close range; attractive long range
                                 - gives us potential of repulsion/attraction

*/

double u2(int nparticles, int i, double x[],double y[]){
  double m = 1.0, r_naught = 1.0;
  double sum,r,term,ratio;
  int j;
  for( i=0,sum=0 ; i < nparticles ; i+=1 ){
    for( j = i+1 ; j < nparticles ; j+=1){
      r = sqrt( (x[i] - x[j])*(x[i] - x[j]) + (y[i]-y[j])*(y[i]-y[j]) ) ;
      term = 0.5 * m * r_naught * r_naught;
      ratio = ( (r-r_naught) * (r-r_naught) ) / (r*r);
      sum += term * ratio;
    }
  }
  return sum;
}

/*

                          potential - returns the total potential energy of the nparticles

*/

double potential(int nparticles, double x[],double y[]){

  //fxn
  double u1(double x, double y);
  double u2(int nparticles, int i, double x[],double y[]);
  //variables
  int i;
  double tot_potential_energy;

  for(tot_potential_energy=0.0, i=0 ; i < nparticles ; i+=1){
    tot_potential_energy += u1(x[i],y[i]); //  u2(nparticles, i, x,y);
  }
  return tot_potential_energy + u2(nparticles,i,x,y);
}

/*

                          kinetic - returns the total kinetic energy of the nparticles

*/

double Kinetic(int nparticles, double vx[], double vy[]){

  int i;
  double m=1.0,tot_kinetic_energy;

  for(tot_kinetic_energy=0.0,i=0; i < nparticles ; i+=1){
    tot_kinetic_energy += 0.5 * m * (vx[i]*vx[i] + vy[i]*vy[i]); //sqrt root and square cancel
  }
  return tot_kinetic_energy;
}


/*

                          Gradient - use the change in position to calculate the velocity step

*/

void Gradient(int nparticles, double x[], double y[], double gradx[],double grady[]){

  //fxn
  //double u1(double x);
  double potential(int nparticles, double x[],double y[]);
  //variables
  int i;
  //variables for x
  double temp_x,potential1,potential2;
  //variables for y
  double temp_y,y_potential1,y_potential2;
  //h step
  double h=1e-6;

  for(i=0; i< nparticles; i+=1){


    temp_x = x[i];
    x[i] = temp_x+h;
    potential1 = potential(nparticles, x,y );
    x[i] = temp_x-h;
    potential2 = potential(nparticles,x,y);
    gradx[i] = ( potential1 - potential2 ) / (2.0*h);
    x[i] = temp_x;

    temp_y = y[i];
    y[i] = temp_y+h;
    y_potential1 = potential(nparticles,x,y);
    y[i] = temp_y-h;
    y_potential2 = potential(nparticles,x,y);
    grady[i] = (y_potential1 - y_potential2) / (2.0*h);
    y[i] = temp_y;
  }
}


/*

                          EulerStep - takes one integration step

*/



void Euler(int nparticles,double x[], double y[], double vx[], double vy[]){
        int i;
	double gradx[N],grady[N];

        Gradient(nparticles,x,y,gradx,grady);

        for(i=0;i<N;i++){
                x[i]+=vx[i]*DT;
                y[i]+=vy[i]*DT;
                vx[i]-=gradx[i]*DT;
                vy[i]-=grady[i]*DT;
        }
}

/*

                          RK2 - takes one integration step

*/


void RK2(int nparticles, double x[], double y[], double vx[], double vy[]){
        int i;
        double xtemp[N],gradx1[N],gradx2[N],x1[N],x2[N],vx1[N],vx2[N];
        double ytemp[N],grady1[N],grady2[N],y1[N],y2[N],vy1[N],vy2[N];

        Gradient(nparticles,x,y,gradx1,grady1);

        for(i=0;i<N;i++){
                x1[i]=vx[i];
                vx1[i]=-gradx1[i];
                xtemp[i]=x[i]+x1[i]*DT;

                y1[i]=vy[i];
                vy1[i]=-grady1[i];
                ytemp[i]=y[i]+y1[i]*DT;
        }

        Gradient(nparticles,xtemp,ytemp,gradx2,grady2);

        for(i=0;i<N;i++){
                x2[i]=vx[i]+vx1[i]*DT;
                vx2[i]=-gradx2[i];
                x[i]+=(DT/2.)*(x1[i]+x2[i]);
                vx[i]+=(DT/2.)*(vx1[i]+vx2[i]);

                y2[i]=vy[i]+vy1[i]*DT;
                vy2[i]=-grady2[i];
                y[i]+=(DT/2.)*(y1[i]+y2[i]);
                vy[i]+=(DT/2.)*(vy1[i]+vy2[i]);
        }
}

/*

                          RK4 - takes one integration step

*/


void RK4(int nparticles, double x[], double y[], double vx[], double vy[]){
        int i;
        double xtemp[N],gradx1[N],gradx2[N],gradx3[N],gradx4[N],x1[N],x2[N],x3[N],x4[N],vx1[N],vx2[N],vx3[N],vx4[N];
        double ytemp[N],grady1[N],grady2[N],grady3[N],grady4[N],y1[N],y2[N],y3[N],y4[N],vy1[N],vy2[N],vy3[N],vy4[N];

        Gradient(nparticles,x,y,gradx1,grady1);

        for(i=0;i<N;i++){
                x1[i]=vx[i];
                vx1[i]=-gradx1[i];
                xtemp[i]=x[i]+x1[i]*DT/2.;

                y1[i]=vy[i];
                vy1[i]=-grady1[i];
                ytemp[i]=y[i]+y1[i]*DT/2.;
        }
        Gradient(nparticles,xtemp,ytemp,gradx2,grady2);

        for(i=0;i<N;i++){
                x2[i]=vx[i]+vx1[i]*DT/2.;
                vx2[i]=-gradx2[i];
                xtemp[i]=x[i]+x2[i]*DT/2.;

                y2[i]=vy[i]+vy1[i]*DT/2.;
                vy2[i]=-grady2[i];
                ytemp[i]=y[i]+y2[i]*DT/2.;
        }
        Gradient(nparticles,xtemp,ytemp,gradx3,grady3);

        for(i=0;i<N;i++){
                x3[i]=vx[i]+vx2[i]*DT/2.;
                vx3[i]=-gradx3[i];
                xtemp[i]=x[i]+x3[i]*DT;

                y3[i]=vy[i]+vy2[i]*DT/2.;
                vy3[i]=-grady3[i];
                ytemp[i]=y[i]+y3[i]*DT;
        }
        Gradient(nparticles,xtemp,ytemp,gradx4,grady4);

        for(i=0;i<N;i++){
                x4[i]=vx[i]+vx3[i]*DT;
                vx4[i]=-gradx4[i];
                x[i]+=(DT/6.)*(x1[i]+2.*x2[i]+2.*x3[i]+x4[i]);
                vx[i]+=(DT/6.)*(vx1[i]+2.*vx2[i]+2.*vx3[i]+vx4[i]);

                y4[i]=vy[i]+vy3[i]*DT;
                vy4[i]=-grady4[i];
                y[i]+=(DT/6.)*(y1[i]+2.*y2[i]+2.*y3[i]+y4[i]);
                vy[i]+=(DT/6.)*(vy1[i]+2.*vy2[i]+2.*vy3[i]+vy4[i]);
        }
}

/*

                          main

*/

int main(){

/*

        The Set Up

*/

  //opening a file to write into
  FILE *fp, *fp1;
  fp = fopen("molecular_dynamics_2d_test.out","w");
  fp1 = fopen("mo_dy_velocities.out","w");
  // initializing fxns
  double potential(int nparticles, double x[N], double y[N]);
  double Kinetic( int nparticles, double vx[N], double vy[N]);
  // initializing variables
  int nparticles,i,k;
  double t, tot_potential_energy,tot_kinetic_energy,final_total_energy;
  double initial_kinect, initial_potent, initial_tot_energy,percentage;
  // initializing arrays
  double x[N], y[N],vx[N], vy[N],gradx[N],grady[N];
  srandom(time(0));
  //srandom(9);
  for(i=0 ; i<N ; i+=1){
    x[i] = (N/4)*(double) (random()/(0.5*MAX_RNDM))-1;
    y[i] = (N/4)*(double) (random()/(0.5*MAX_RNDM))-1;
    vx[i] = 0.0;
    vy[i] = 0.0;
    }
  // number of nparticles
  nparticles = N;

  initial_potent = potential(nparticles, x,y);
  initial_kinect = Kinetic(nparticles,vx,vy);
  initial_tot_energy = initial_potent + initial_kinect;

  // intigrating the motion
  for(t=0 ; t < 40.0*M_PI ;t+=DT){
        tot_potential_energy = potential(nparticles, x,y);
        tot_kinetic_energy = Kinetic(nparticles,vx,vy);
        final_total_energy = tot_potential_energy + tot_kinetic_energy;
    fprintf(fp,"%0.4e\t %0.5e\t %0.5e\t  %0.5e\t %0.5e\t  %0.5e\t %0.5e\t  %0.5e\t %0.5e\t  %0.5e\t %0.5e\t  %0.5e\t %0.5e\t  %0.5e\t %0.5e\t  %0.5e\t %0.5e\t  %0.5e\t %0.5e\t   %0.5e\t  %0.5e\t %0.5e\t %0.5e\t  %0.5e\n", t,x[0],y[0],x[1],y[1],x[2],y[2],x[3],y[3],x[4],y[4],x[5],y[5],x[6],y[6],x[7],y[7],x[8],y[8],x[9],y[9],final_total_energy,tot_potential_energy,tot_kinetic_energy);
// choose ur fav. integrator: Euler, RK2, or RK4
    RK4(nparticles,x,y,vx,vy);
    }
    percentage = ( final_total_energy- initial_tot_energy)/initial_tot_energy;
    printf("percent energy change %0.5e\n  ", percentage*100   );
    // outputing our velocities to a file
//    for(i=0 ; i < nparticles ; i+=1){
//      fprintf(fp1, "%d\t  %0.5e\n ",i, (sqrt(vx[i]*vx[i] + vy[i]*vy[i])) );
//    }


    /*
      Calculating mean and std dev of position and velocity of particles
    */

    //variables
    int j;
    double pos_sum=0,vel_sum=0,sqrd_pos_sum=0, sqrd_vel_sum=0;
    double pos_mean,vel_mean,x_stdev,v_stdev,sigma_x,sigma_v;
    //
    for(j=0 ; j<nparticles-1 ; j+=1){
      pos_sum += sqrt(x[j]*x[j] + y[j]*y[j]);
      vel_sum += sqrt(vx[j]*vx[j] + vy[j]*vy[j]);
      sqrd_pos_sum += (pos_sum*pos_sum);
      sqrd_vel_sum += (vel_sum*vel_sum);
    }

    //mean
    pos_mean = pos_sum / nparticles;
    vel_mean = vel_sum / nparticles;
    //stdev
    for(i = 0 ; i < nparticles -1; i+=1){
    sigma_x =   (sqrt(x[i]*x[i] + y[i]*y[i]) - pos_mean)*(sqrt(x[i]*x[i] + y[i]*y[i]) - pos_mean) ;
    sigma_v =   (sqrt(vx[i]*vx[i] + vy[i]*vy[i]) - pos_mean)*(sqrt(vx[i]*vx[i] + vy[i]*vy[i]) - pos_mean) ;
    }
  //  x_stdev = sqrt(( sqrd_pos_sum - (nparticles*pos_mean*pos_mean) )/(nparticles));
    x_stdev = sqrt(sigma_x/nparticles);
  //  v_stdev = sqrt(( sqrd_vel_sum - (nparticles*vel_mean*vel_mean))/(nparticles-1));
    v_stdev = sqrt(sigma_v /nparticles);
    printf("avg. position: %0.5e\t stdev. position: %0.5e \n", pos_mean, x_stdev);
    printf("avg. velocity: %0.5e\t  stdev. velocity: %0.5e \n", vel_mean, v_stdev );

    //putting final velocities into bins
    int l,dist_x[21],dist_y[21], v_ex[N],v_why[N],toobig, toosmall; // writing the bins according to the handout
        for(i=0;i<21;i++){ // we have 20 bins
            dist_x[i]=0;
            dist_y[i]=0;
        }
        for(i=0 ; i< nparticles ; i+=1){
        //  v[i] = sqrt( vx[i]*vx[i] + vy[i]*vy[i] );
          v_ex[i] = vx[i];
          v_why[i]=vy[i];
        }
        for(i=0;i<nparticles;i++){
            v_ex[i]+=10.5;
            v_why[i]+=10.5;
            if(v_ex[i]>21) toobig++;
            else if(v_ex[i]<0) toosmall++;
            else dist_x[(int)v_ex[i]]++;
            if(v_why[i]>21) toobig++;
            else if(v_why[i]<0) toosmall++;
            else dist_y[(int)v_why[i]]++;
        }
        for(i=0;i<21;i++){
            l=i-10;
            // uncomment below to print out bins
            printf("%d\t%d\t%d\n",l,dist_x[i],dist_y[i]); // print out the elements in each bins
            fprintf(fp1, "%d\t %d\t  %d\n ",l, dist_x[i],dist_y[i] );

        }
/*

                How many particles escaped the well?

*/

int im_freeee;
double position,ke,trap_wall;
double u1(double x, double y);

trap_wall = 50;

for(im_freeee = 0,i=0 ; i < nparticles ; i+=1){
  position = sqrt(x[i]*x[i] + y[i]*y[i]);
  ke = 0.5*(vx[i]*vx[i] + vy[i]*vy[i]);
  if(position > trap_wall && ke > u1(x[i],y[i])){
    im_freeee +=1;
  }
}
printf("%d particles escaped\n", im_freeee);

fclose(fp);
fclose(fp1);

}
