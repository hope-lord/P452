#include <cstdlib>
#include <iostream>
#include <cmath>

using namespace std;
int a = 572; //multiplier for gen2
int m = 16381; //modulus for gen2

int x0 = 94; //seed 

int my_random(){
	int l = (a*x0)%m;
	x0 = l;
	return l;
}

double rnd_db(double a, double b){
	return (b-a)*my_random()/m + a;
	/* return (b-a)*(double)random()/RAND_MAX + a; */
}

int main(){
/*	Given that the radius of the cylinder is 1. So the direct distance between two points the Steinmetz solid is less than 2. So by taking the cube of side length 2 we can simulate the montecarlo.*/
 	/* Assume that one cylinder as axis along z axis and another one has axis along y axis*/
	int N = 10000000; // Max No of points to be thrown
	for(int n=1000;n<=N;n*=10){
		int sum = 0;
		for(int i=0;i<n;i++){
			double x = rnd_db(-1,1);
			double y = rnd_db(-1,1);
			double z = rnd_db(-1,1);
			bool isINSIDE = (x*x+y*y<=1) and (x*x+z*z<=1);
			if(isINSIDE) sum++;
		}
		double volume = 8*(double)sum/n;
		printf ("For %8d points thrown volume = %lf\n",n,volume);
	}
	return 0;
}
