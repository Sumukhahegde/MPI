#include <stdio.h>
#include <math.h>


double f(double x){

	double return_val;
	
	return_val = pow((4.0 - x*x),0.5);

	return return_val;	
}

main() {

	double integral;
	double a,b;
	int n, i;
	double h;
	double x;



	printf("Enter a, b, and n\n");
	scanf("%lf %lf %d", &a, &b, &n);


	h = (b - a)/n;

	integral = (f(a) + f(b))/2.0;
	
	x = a;

	for(i=1;i<=n-1;i++){
		x = x+h;	
		integral = integral + f(x);
	}
	
	integral = integral * h;
	
	printf(" With %d trapezoids, our estimate \n", n);

	printf("of the integral from %lf to %lf = %lf\n", a, b, integral);

}


 
	
	

		
	
