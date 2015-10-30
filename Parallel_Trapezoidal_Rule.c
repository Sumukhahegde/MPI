#include <stdio.h>
#include <math.h>
#include "mpi.h"

double f(double x){
	
	double height;
	height = pow((4.0 - x * x), 0.5);

	return height;
}

double Trap(double local_a, double local_b, long long int local_n, double h){

	double integral;
	double x;
	double i;	

	integral = (f(local_a) + f(local_b))/2.0;

	x = local_a;	

	for(i = 1;i <= local_n - 1;i++){
		x = x + h;
		integral = integral + f(x);
	}

	integral = integral * h;

	return integral;
}

main(int argc, char** argv){

	int p;
	int my_rank;
	double a = -2.0;
	double b = 2.0;
	long long int n = pow(2, 30);
	double h;	//base lenght of each trapezoid
	double local_a;
	double local_b;
	long long int local_n;
	double integral;
	double total;
	int source;
	int dest = 0;
	int tag = 0;
	MPI_Status status;

//	double Trap(double local_a, double local_b, double local_n, double h);

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	MPI_Comm_size(MPI_COMM_WORLD, &p);

	h = ( b - a ) / n;

	local_n = n / p;

	local_a = a + my_rank * local_n * h;

	local_b = local_a + local_n * h;

	integral = Trap(local_a, local_b, local_n, h);


	if(my_rank==0){
		total = integral;

		for(source = 1;source<p;source++){
			MPI_Recv(&integral, 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
			total += integral;
		}
	} else{
		MPI_Send(&integral, 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
	}

	if(my_rank == 0){
		printf("With %lld trapezoids, our estimate of the integral from %lf to %lf = %lf\n", n, a, b, total);
	}	  
	
	MPI_Finalize();
}






		
	
		 			
	
		

		 
	

