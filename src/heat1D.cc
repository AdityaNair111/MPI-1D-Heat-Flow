#include<iostream>
#include<fstream>
#include<iomanip>
#include <limits>
#include"mpi.h"
#include<math.h>
using namespace std;



int main(int argc, char *argv[])
{
	double r=0.25;
	MPI_Init(&argc,&argv);
	int numProcs,proc_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

	int NumGridPoints,numTimeSteps;
	double T1,T2;

	// check if the number of argumnets are correct

	if ((argc != 5) && (proc_rank==0))
	{
		cout<<"There seems to be a problem with number of arguments"<<endl;
		MPI_Finalize();
		exit(EXIT_FAILURE);
	}
	NumGridPoints = atoi(argv[3]);
	numTimeSteps = atoi(argv[4]);
	T1 = atof(argv[1]);
	T2 = atof(argv[2]);

	if ((NumGridPoints < 1)&&(proc_rank==0))
	{
	cout<<"Number of points must atleast be 1! Aborting!"<<endl;
		MPI_Finalize();
		exit(EXIT_FAILURE);
	}

	if ((numProcs > NumGridPoints) && (proc_rank==0))
	{
		cout<<"Processors requested greater than tasks! but its okay"<<endl;
	}

	
	
	int start_points[numProcs];
	int end_points[numProcs];
	int num_points[numProcs];
	
	if (numProcs<=NumGridPoints)
	{
		int approxPointsPerProc = NumGridPoints/numProcs;
		int remPoints= NumGridPoints % numProcs;

		for (int i = 0; i < numProcs; i++)
		{
			int start,end,num;
			if(i<remPoints)
			{
				start = i*ceil( float(NumGridPoints)/float(numProcs));
				end = i + ceil( float(NumGridPoints)/float(numProcs)) - 1;
				num = ceil( float(NumGridPoints)/float(numProcs));
			}
			else
			{
				start = i*approxPointsPerProc + remPoints;
				end = start + approxPointsPerProc - 1;
				num = end-start	+ 1;
			}

			num =num+1;
			end =end+1;
			num =num+1;
			start_points[i]=start;
			end_points[i]=end;
			num_points[i]=num;
		}
	}
	else
	{
	
		for (int i = 0; i < numProcs; i++)
		{
			if (i<NumGridPoints)
			{
				start_points[i]=i;
				num_points[i]=3;
			}
			else
			{
				start_points[i]=i;
				num_points[i]=0;
			}
		}
		numProcs=NumGridPoints;
	}
	NumGridPoints=NumGridPoints+2;
	

	int current_num = num_points[proc_rank];
	int current_start = start_points[proc_rank];
	int current_end = end_points[proc_rank];

	double current_grid[current_num];
	double temp_current_grid[current_num];
	
	if (proc_rank==0)
	{
		
		current_grid[0]=T1;
		for (int i = 1; i < current_num; i++)
		{
			current_grid[i]=0;
		}
		if (numProcs==1)
		{
			current_grid[current_num-1]=T2;
		}
	}
	else
	{
		for (int i = 0; i < current_num; i++)
		{
			current_grid[i]=0;
		}
		if (proc_rank==numProcs-1)
		{
			current_grid[current_num-1]=T2;
		}
	
	}	

	for (int t = 0; t < numTimeSteps; t++)
	{
		for (int x = 1; x < current_num - 1; x++)
		{
			temp_current_grid[x]=(1-2*r)*current_grid[x]+r*current_grid[x-1]+r*current_grid[x+1];
		}
		

		if (((proc_rank % 2)==0)&&(proc_rank<numProcs))
		{
			if (proc_rank<numProcs-1)
			{
				MPI_Send(&temp_current_grid[current_num-2],1,MPI_DOUBLE,proc_rank+1,2,MPI_COMM_WORLD);
			}
			if (proc_rank>0)
			{
				MPI_Send(&temp_current_grid[1],1,MPI_DOUBLE,proc_rank-1,2,MPI_COMM_WORLD);
			}
			if (proc_rank<numProcs-1)
			{
				MPI_Status status;
				MPI_Recv(&temp_current_grid[current_num-1],1,MPI_DOUBLE,proc_rank+1,2,MPI_COMM_WORLD,&status);
			}
			if (proc_rank>0)
			{
				MPI_Status status;
				MPI_Recv(&temp_current_grid[0],1,MPI_DOUBLE,proc_rank-1,2,MPI_COMM_WORLD,&status);
			}
		}
		else
		{
			if ((proc_rank<numProcs))
			{
				
				if (proc_rank<numProcs-1)
				{
					MPI_Status status;
					MPI_Recv(&temp_current_grid[current_num-1],1,MPI_DOUBLE,proc_rank+1,2,MPI_COMM_WORLD,&status);	
				}
				if (proc_rank>0)
				{
					MPI_Status status;
					MPI_Recv(&temp_current_grid[0],1,MPI_DOUBLE,proc_rank-1,2,MPI_COMM_WORLD,&status);
				}
				if (proc_rank<numProcs-1)
				{
					MPI_Send(&temp_current_grid[current_num-2],1,MPI_DOUBLE,proc_rank+1,2,MPI_COMM_WORLD);
				}
				if (proc_rank>0)
				{
					MPI_Send(&temp_current_grid[1],1,MPI_DOUBLE,proc_rank-1,2,MPI_COMM_WORLD);
				}

			}

		}
		if (proc_rank==0)
		{
			temp_current_grid[0]=T1;
			
		}
		if (proc_rank==numProcs-1)
		{
			temp_current_grid[current_num-1]=T2;
		}
			
		
		MPI_Barrier(MPI_COMM_WORLD);
		
		for (int i = 0; i < current_num; i++)
		{
			current_grid[i]=temp_current_grid[i];
		}

	}

	if (proc_rank==0)
	{
		double grid[NumGridPoints];
		
		grid[0] = T1;
		
		grid[NumGridPoints-1] = T2;
		// for 1st processor
		for (int i = 1; i <= current_num - 2 ; i++)
		{
			grid[i]=current_grid[i];
		}
		
		//everything in between
		for (int i = 1; i <= numProcs-1; i++) // starting from 1
		{
			MPI_Status status;
			
			MPI_Recv(&grid[0]+start_points[i]+1,num_points[i]-2,MPI_DOUBLE,i,3,MPI_COMM_WORLD,&status);
				
		}
		for (int i = 0; i < NumGridPoints ; i++)
		{
		 	
		 	float temperature = grid[i];
	 		fprintf(stdout, "%.8f\n",temperature);
		}
	
		ofstream build ("heat1Doutput.csv", std::ofstream::out);
		for (int i = 1; i < NumGridPoints-2; ++i)
		{
			build <<std::fixed<<std::setprecision(std::numeric_limits<long double>::digits10 + 1)<<grid[i]<<", ";
		}
		build<<grid[NumGridPoints-2];
		build.close();
		
	}
	else
	{
		if(proc_rank<numProcs)
		{
		int points=current_num-2;
		
		MPI_Send(&current_grid[1],points,MPI_DOUBLE,0,3,MPI_COMM_WORLD);
		
		}
	}

	MPI_Finalize();
	return EXIT_SUCCESS;

}
