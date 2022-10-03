#include"libA.h"
#include"SmoothParticle.h"
#include<iostream>
using namespace std;
int tmain(int argc, char* argv[]);
int runmain(int argc, char* argv[]);
int main(int argc, char* argv[])
{
	WORKER* wkr;

	MPI_Init(&argc, &argv);
	WorkerStart(&wkr);
	mpiwait(wkr, argc, argv);
	//int code = tmain(argc, argv);
	int code1 = runmain(argc, argv);
	MPI_Finalize();
	return code1;
}

#include"SmoothParticle.h"

int tmain(int argc, char* argv[])
{
	
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int n = sizeof(ParticleProperties);
	ParticleProperties particle;
	if (rank == 0)
	{
		particle.x = 1;
		particle.y = 2;
	}

	vector<ParticleProperties> a(20,particle);
	if (rank == 0)
	{
		MPI_Send(&a[0], n * 0, MPI_CHAR, 1, 0, MPI_COMM_WORLD);
	}
	MPI_Status stats;
	if (rank == 1)
	{
		cout << "RANK " << rank << ' ' << a[1].x << ' ' << a[19].y << endl;
		MPI_Recv(&a[0], n * 0, MPI_CHAR, 0, 0, MPI_COMM_WORLD,&stats);
		cout << "RANK " << rank << ' ' << a[1].x << ' ' << a[19].y << endl;
	}




	return 0;

}

int runmain(int argc, char* argv[])
{
	double exv = 1e0;
	WaterProperty water;
	water.K = pow(5.1e-10,-1.0)/3/1000 / 1e4;
	water.rho0 = 1000;
	water.mu = 1e-3*(exv)*1000;
	//water.K = 200 * water.rho0 * 9.8 * 1;
	WallBoxProperty box;
	box.g = 9.8*1e0*exv;//?????
	box.x0 = 0;
	box.y0 = 0;
	box.x1 = 1;
	box.y1 = 1;
	box.p = 0;
	box.T = 0;
	box.vx = 0;
	box.vy = 0;

	SPWater spw(water,box,0.02,0.04);
	spw.Initialize(cout);
	//spw.Run(100000, "D:\\CAResults\\SPH\\out",cout);
	spw.Run(100000, "./SPH/out", cout);
	return 0;
}