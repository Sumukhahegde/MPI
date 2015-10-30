#include <stdio.h>
#include <math.h>
#include "mpi.h"



void Get_data(double *a_ptr, double *b_ptr, long long int *n_ptr, int my_rank, int p){

	int source = 0;
	int dest;
	int tag;
	MPI_Status status;	

	if(my_rank == 0){
		printf("Enter the starting point a, ending point b and the number of trapezoids:");
		scanf("%lf %lf %lld",a_ptr, b_ptr, n_ptr);

		for(dest = 1; dest<p; dest++){
			tag = 0;
			MPI_Send(a_ptr, 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
			tag = 1;
			MPI_Send(b_ptr, 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
			tag = 2;
			MPI_Send(n_ptr, 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
		}
	} else{
		tag = 0;
		MPI_Recv(a_ptr, 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
		tag = 1;
		MPI_Recv(b_ptr, 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
		tag = 2;
		MPI_Recv(n_ptr, 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
	}

}

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
	double a;
	double b;
	long long int n;
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
	
	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	MPI_Comm_size(MPI_COMM_WORLD, &p);


	Get_data(&a, &b, &n, my_rank, p); 


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






		
	
		 			
	
		

		 
	

