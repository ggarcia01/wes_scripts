#include <stdio.h>
#include <math.h>
int main(){
	int i;
	float m,k,xi,vi,x,v,dx,dt,dv,t,c,diff;
	float a,omega,phi;
	FILE *fp;
	fp=fopen("eq_of_motion.out","w");
	printf("enter xi, vi,m,k: ");
	scanf("%f%f%f%f",&xi,&vi,&m,&k);
	printf("Enter amp, ang freq, and phase offset: ");
	scanf("%f%f%f",&a,&omega,&phi);
	dt = 0.00001;
	x = xi;
	v = vi;
	for(i=0;i<600000;i+=1){
		dv = (-k*x*dt)/m;
		v = v + dv;
		dx = v*dt;
		x = x + dx;
		t = t + dt;
		c = a*cosf(omega*t+phi);
		diff = x-c ;
		fprintf(fp,"%f\t%f\t%f\t%f\n",t,x,c,diff);
	}
	fclose(fp);
}	
