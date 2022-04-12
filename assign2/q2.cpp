#include <cstdlib>
# include <iostream>
# include <cmath>

int x1 = 13;  //seed for gen1
int x2 = 101;  //seed for gen2

int a1 = 65; //multiplier for gen1
int m1 = 1021; //modulus for gen1

int a2 = 572; //multiplier for gen2
int m2 = 16381; //modulus for gen2


double my_rand1(double min,double max){
	x1 = (a1 * x1) % m1;
	return min + (max - min) * x1 / m1;
}

double my_rand2(double min,double max){
	x2 = (a2 * x2) % m2;
	return min + (max - min) * x2 / m2;
}


// Simulate value of pi by throwing dart method and using my_rand1 as random number generator
double pi_with_throwpoints_gen1(int n){
	double inside = 0;
	for(int i = 0; i < n; i++){
		double x = my_rand1(-0.5,0.5);
		double y = my_rand1(-0.5,0.5);
		if(x*x + y*y <= 0.25){
			inside++;
		}
	}
	//The circle is inscribed in a square with side length 1
	return 4 * inside / n;
}

// Simulate value of pi by throwing dart method and using my_rand2 as random number generator
double pi_with_throwpoints_gen2(int n){
	double inside = 0;
	for(int i = 0; i < n; i++){
		double x = my_rand2(-0.5,0.5);
		double y = my_rand2(-0.5,0.5);
		if(x*x + y*y <= 0.25)
			inside++;
	}
	//The circle is inscribed in a square with side length 1
	return 4 * inside / n;
}

double func(double x){
	return sqrt(1 - x*x);
}

// Calculate the integral of func from 0 to 1 using Monte Carlo method with my_rand1 as random number generator

double pi_with_MC_integration_gen1(int n){
	double inside = 0;
	for(int i = 0; i < n; i++){
		double x = my_rand1(0,1);
		inside += func(x);
	}
	return inside / n * 4;
}

// Calculate the integral of func from 0 to 1 using Monte Carlo method with my_rand1 as random number generator
double pi_with_MC_integration_gen2(int n){
	double inside = 0;
	for(int i = 0; i < n; i++){
		double x = my_rand2(0,1);
		inside += func(x);
	}
	return inside / n * 4;
}


int main(){
	std::cout << "pi with throwpoints with a = 65, m = 1021   : " << pi_with_throwpoints_gen1(1000000) << std::endl;
	std::cout << "pi with throwpoints with a = 572, m = 16381 : " << pi_with_throwpoints_gen2(1000000) << std::endl;
	std::cout << "pi with MC integration a = 65, m = 1021      : " << pi_with_MC_integration_gen1(1000000) << std::endl;
	std::cout << "pi with MC integration a = 572, m = 16981    : " << pi_with_MC_integration_gen2(1000000) << std::endl;
	return 0;
}



/*
pi with throwpoints with a = 65, m = 1021   : 3.21568
pi with throwpoints with a = 572, m = 16381 : 3.12187
pi with MC integration a = 65, m = 1021      : 3.14269
pi with MC integration a = 572, m = 16981    : 3.1417
*/
