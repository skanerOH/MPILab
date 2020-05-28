#include "pch.h"
#include <mpi.h>
#include <stdio.h>
#include <chrono>
#include <iostream>
#include <cmath>

inline double F(double x)
{
	return cos(sin(x*x) - x * x*x + x * x*x*x / 47 - 54 * cos(pow(x, -1.5)));
}


int main(int argc, char** argv) {


	auto start_time = std::chrono::high_resolution_clock::now();
	int pid, np;
	const double begin = 1;
	const double end = 200;
	const long int values_count = 1000001;


	MPI_Init(NULL, NULL);

	// Get the number of processes
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	// Get the rank of the process
	MPI_Comm_rank(MPI_COMM_WORLD, &pid);
	const long int calcs_per_process = (values_count / np);

	if (pid == 0)
		//main process begins here
	{
		double* web = new double[values_count];
		web[0] = begin;
		for (int i = 1; i < values_count-1; ++i)
			web[i] = begin + ((end - begin) / (values_count-1))*i;
		web[values_count - 1] = end;

		if (np > 1)
		{
			int i=1;
			long int from_ind = 0;
			long int to_ind = calcs_per_process;
			long int values_to_send = calcs_per_process + 1;
			for (i = 1; i < np - 1; ++i)
			{
				//sending values count
				MPI_Send(&values_to_send, 1, MPI_LONG, i, 0, MPI_COMM_WORLD);
				//sending values
				MPI_Send(&(web[i*(calcs_per_process)]), values_to_send, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
				//printf("from %d to %d\n", i*calcs_per_process, i*calcs_per_process - 1 + values_to_send);
			}

			long int values_left_to_send = values_count - i * calcs_per_process;
			MPI_Send(&values_left_to_send, 1, MPI_LONG, i, 0, MPI_COMM_WORLD);
			MPI_Send(&(web[i*calcs_per_process]), values_left_to_send, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
			//printf("from %d to %d\n", i*calcs_per_process, i*calcs_per_process - 1 + values_left_to_send);
		}
		//printf("Sent values");

		//main process calculations
		double res = 0;
		double a, b;
		long int values_last_ind_to_calc_itself = ((np == 1) ? calcs_per_process - 1 : calcs_per_process);
		for (long int i = 0; i < values_last_ind_to_calc_itself; ++i)
		{
			a = web[i];
			b = web[i + 1];
			res += (b - a)*(F((a + b) / 2.0 - (b - a) / (2 * sqrt(3))) + F((a + b) / 2.0 + (b - a) / (2 * sqrt(3)))) / 2;
		}
		//printf("Calculated itself");


		double returned_integral;
		MPI_Status status;
		for (int i = 1; i < np; i++) {
			MPI_Recv(&returned_integral, 1, MPI_DOUBLE,
				MPI_ANY_SOURCE, 0,
				MPI_COMM_WORLD,
				&status);
			int sender = status.MPI_SOURCE;
			res += returned_integral;
		}

		printf("Integral = %f\n", res);

		delete[] web;
	}
	else
		//"child" process way
	{
		//printf("Hello from child");
		long int recieved_values_count;
		MPI_Status status;

		MPI_Recv(&recieved_values_count,
			1, MPI_LONG, 0, 0,
			MPI_COMM_WORLD,
			&status);

		double* recieved_values = new double[recieved_values_count];

		MPI_Recv(recieved_values, recieved_values_count,
			MPI_DOUBLE, 0, 0,
			MPI_COMM_WORLD,
			&status);

		double prores = 0;
		double a, b;
		for (long int i = 0; i < recieved_values_count-1; ++i)
		{
			a = recieved_values[i];
			b = recieved_values[i+1];
			prores += (b - a)*(F((a + b) / 2.0 - (b - a) / (2 * sqrt(3))) + F((a + b) / 2.0 + (b - a) / (2 * sqrt(3)))) / 2;
		}

		MPI_Send(&prores, 1, MPI_DOUBLE,
			0, 0, MPI_COMM_WORLD);

		delete[] recieved_values;
	}

	MPI_Finalize();
	if (pid == 0)
	{
		auto stop_time = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop_time - start_time);

		std::cout << "Time spent: " << duration.count() << std::endl;
	}
	return 0;
}
