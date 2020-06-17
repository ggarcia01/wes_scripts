// march 27 2019
// in this file, we calculate the period of the pendulum as a function of angles
#include <stdio.h>
#include<math.h>
#define DT 0.000001 // we found that this step size gives sufficiently accurate results
int main(){
	int n,j; //these variables will serve as counters
	double x1,x11,oldx1,x2,t,oldt,y,y1,period,approx,diff,q; // the extra variables needed
	double Interp(double x_1,double y_1,double x_2,double y_2,double x); //we will want to call on our Interp function
	// thus, we define it here so that it is ready when needed
	FILE* fp; //opening a file
	fp = fopen("damped_pendulum.out","w"); //we will be writing out to this file
	for(x11=0.1;x11<M_PI;x11=x11+0.1 ){ //the bigger for loop will iterate over different angles in step of 0.1
		for(x1=x11,n=0,j=0,x2=0,t=0,q=1.5;t<6.*M_PI;t+=DT,n++){ //for each angle, we will find its period, q is the damping factor
			oldx1 = x1; //save the old x value before updating it
			x1+=x2*DT; //updating x1, one of the 1st order diff eq
			x2-=sin(oldx1)*DT + (1/q)* DT; //the second diff eq that was found
			if(x1<0 && oldx1>0 && j==0){ //when we cross the x axis, enter this if statement
				y = Interp(oldx1,oldt,x1,t,0); //we will calc the x-intercept
				j+=1; //in this iteration, we will not enter this if statement again. thus preserving out y value
			}
			else if(x1<0 && oldx1>0){ //next time we cross the x-axis...
				y1 = Interp(oldx1,oldt,x1,t,0); //calculate the next x-intercept
				break; //we have out values needed to calculate the period, so lets get out of this for loop
			}
			oldt = t; //save the previous t value before updating it
		}
	period=fabs(y1-y); //calc period
	fprintf(fp,"%0.5e\t%0.5e\n",x11,period); //printing to the file
	}
	fclose(fp); //closing file
}
