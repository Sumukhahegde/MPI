#include <stdio.h>
#include <math.h>
#include "mpi.h"

void Get_data1(double *a_ptr, double *b_ptr, long long int *n_ptr, int my_rank){

	if(my_rank == 0){
		printf("Enter a, b and n\n");
		scanf("%lf %lf %lld", a_ptr, b_ptr, n_ptr);
	}

	MPI_Bcast(a_ptr, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
	MPI_Bcast(b_ptr, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
	MPI_Bcast(n_ptr, 1, MPI_LONG, 0, MPI_COMM_WORLD);

} 

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
	printf("I am process number %d", p);
			MPI_Send(b_ptr, 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
			tag = 2;
			MPI_Send(n_ptr, 1, MPI_LONG, dest, tag, MPI_COMM_WORLD);
		}
	} else{
		tag = 0;
		MPI_Recv(a_ptr, 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
		tag = 1;
		MPI_Recv(b_ptr, 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
		tag = 2;
		MPI_Recv(n_ptr, 1, MPI_LONG, source, tag, MPI_COMM_WORLD, &status);
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
	double a = -2.0;
	double b = 2.0;
	long long int n = 1000;
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
//	printf("I am process number %d", p);
	
	MPI_Init(&argc, &argv);
//	printf("I am process number %d", p);

	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	MPI_Comm_size(MPI_COMM_WORLD, &p);

	printf("I am process number %d\n", my_rank);

	Get_data1(&a, &b, &n, my_rank);
//	Get_data(&a, &b, &n, my_rank, p); 


	h = ( b - a ) / n;

	local_n = n / p;

	local_a = a + my_rank * local_n * h;

	local_b = local_a + local_n * h;

	integral = Trap(local_a, local_b, local_n, h);

	MPI_Reduce(&integral, &total, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	

	if(my_rank == 0){
		printf("With %lld trapezoids, our estimate of the integral from %lf to %lf = %lf\n", n, a, b, total);
	}	  
	
	MPI_Finalize();
}






		
	
		 			
	
		

		 
	

