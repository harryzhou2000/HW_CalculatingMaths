#pragma once

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <list>
#include <math.h>
#ifdef _MSC_VER
#include<Windows.h>
#include"mpi.h"
//#pragma comment(lib,"C:\\Program Files (x86)\\MPI\\Lib\\x64\\msmpi.lib")
#endif
#ifdef __GNUC__
#include <unistd.h>
#include<mpich/mpi.h>

#endif
//#include <cblas.h>


typedef double SPH_REAL;
typedef int SPH_INT;

#define vector std::vector
#define string std::string
#define list std::list
#define ostream std::ostream
#define SPH_MPI_INT MPI_INT

//template<typename T1>
//constexpr inline auto IBR(T1 i, T1 j, T1 nc) { return ((i)*(nc)+(j)); }
#define IBR(i,j,nc) ((i)*(nc)+(j))
template<typename T1>
constexpr inline auto sqr(T1 i) { return i * i; }
template<typename T1>
constexpr inline auto cube(T1 i) { return  i * i * i; }

struct ParticleProperties
{
	SPH_REAL x;//global!!
	SPH_REAL y;
	SPH_REAL vx;
	SPH_REAL vy;
	SPH_REAL rho;
	SPH_REAL p;
	SPH_REAL ax;
	SPH_REAL ay;
	size_t gridx;
	size_t gridy;
};

struct RhoPass
{
	SPH_REAL rho;
	size_t gridx;
	size_t gridy;
	inline void operator = (ParticleProperties&R)
	{
		rho = R.rho, gridx = R.gridx, gridy = R.gridy;
	}
};



struct SmoothParticleSet
{
	size_t nParticles;
	list<ParticleProperties> p;
	//vector<ParticleProperties*> ppointers;

	
	//vector<bool> particleValidity;

	size_t nOParticles;//particles determined in other processes
	vector<ParticleProperties> outerones[8];
	list<ParticleProperties> solidwalls;
	

	
	size_t nGridx; //ncols;
	size_t nGridy;
	size_t nGrid;
	size_t iGridxm;
	size_t iGridym;
	
	vector<vector<size_t>> particlesInGrid; // grid[IBR(gy+1,gx+1,ngridx)]
	vector<vector<ParticleProperties*>> particlesInGridpointers;
	vector<unsigned char> wallReflector;
#define SPRF_Upper 0x10
#define SPRF_Lower 0x20
#define SPRF_Right 0x40
#define SPRF_Left  0x80
#define SPRF_UpperRight 0x01
#define SPRF_LowerRight 0x02
#define SPRF_UpperLeft 0x04
#define SPRF_LowerLeft  0x08


	vector<SPH_REAL> gridVecx;
	vector<SPH_REAL> gridVecy;
	//if the grid is not outer: index in ppointers;
	//if the grid is exchangeouter: index in "op"s
	vector<unsigned char> gridType;
#define SPGT_Common 0x00
#define SPGT_Border 0x01
#define SPGT_Solid 0x02
#define SPGT_SolidWall 0x03

#define SPGT_PressureOuter 0x04
#define SPGT_ExchangeOuter 0x0E
#define SPGT_FreeOuter 0x05
#define SPGT_VelocityOuter 0x06


#define SPGTInnerType 0x0F// switch( gridType & SPGT_BoundaryType)

#define SPGT_Upper 0x10
#define SPGT_Lower 0x20
#define SPGT_Right 0x40
#define SPGT_Left  0x80
#define SPGTCommType 0xF0 // switch( gridType & ....)



	SPH_REAL xabs0;
	SPH_REAL yabs0;
	SPH_REAL xorig;
	SPH_REAL yorig;
	SPH_REAL xmax;
	SPH_REAL ymax;//absolute x y
	SPH_REAL h;//particle radius and grid size
	SPH_REAL m;
	

	SmoothParticleSet(SPH_REAL nh,SPH_REAL nm) { 
		nOParticles = nParticles = nGridx = nGridy = nGrid = 0;
		h = nh; m = nm; xabs0 = yabs0 = xmax = ymax = xorig = yorig = 0;
	}
	inline size_t getGridx(SPH_REAL xi) { return size_t(ceil((xi - xorig) / h)); }
	inline size_t getGridy(SPH_REAL yi) { return size_t(ceil((yi - yorig) / h)); }
	inline SPH_REAL localx(SPH_REAL absx) { return absx - xorig; }
	inline SPH_REAL localy(SPH_REAL absy) { return absy - yorig; }
	inline bool isInside(size_t gx, size_t gy) { return gx >= 0 && (gx < nGridx - 1) && gy >= 0 && (gy < nGridy - 1); }
};

struct IdealGasProperty
{
	SPH_REAL Rg;
	SPH_REAL gamma;
	SPH_REAL mu1;
	SPH_REAL mu2;
};

struct WaterProperty
{
	SPH_REAL K;
	SPH_REAL rho0;
	SPH_REAL mu;
};


struct InfiniteFieldBoundaryBoxProperty
{
	SPH_REAL x0;
	SPH_REAL y0;
	SPH_REAL x1;
	SPH_REAL y1;
	SPH_REAL p;//@infinity
	SPH_REAL vx;
	SPH_REAL vy;
	SPH_REAL T;
};

struct WallBoxProperty
{
	SPH_REAL x0;
	SPH_REAL y0;
	SPH_REAL x1;
	SPH_REAL y1;
	SPH_REAL g;

	SPH_REAL p;//@infinity
	SPH_REAL vx;
	SPH_REAL vy;
	SPH_REAL T;
};


class SPIsentropicAero
{
	SmoothParticleSet set;
	IdealGasProperty gas;
	InfiniteFieldBoundaryBoxProperty box;
	SPH_REAL cpressure;//p = cpressure rho^gamma
	SPH_REAL cpRho;//rho = cprho sum(prho)
	SPH_REAL cpAp; //ap = cpAp sum(pAp)
	SPH_REAL cpAvis; //avis = cpAvis sum(avis)
	SPH_REAL cviss;
	SPH_REAL inletDensityX;
	SPH_REAL inletDensityY;
	size_t boxtype;

	SPH_REAL dt;//time step
	SPH_REAL hinit;//for spawn and initialization 


	int nproc;
	int rank;
	int rankc;//cartesian rank
	int rankx;//cartesian x, y rank	
	int ranky;
	int rankupper;
	int ranklower;
	int rankright;
	int rankleft;
	MPI_Comm commc;


	
public:
	SPIsentropicAero(IdealGasProperty ngas, InfiniteFieldBoundaryBoxProperty nbox, SPH_REAL nh, SPH_REAL nm)
		:gas(ngas), box(nbox), set(nh,nm) {}

	//inline DOUBLE pRho(DOUBLE r);//nonconstant part of summation: rho//note that rho need to add rho0
	//inline DOUBLE pAP(DOUBLE r, DOUBLE xi, DOUBLE xj, DOUBLE pi, DOUBLE pj, DOUBLE rhoi, DOUBLE rhoj); //p acceleration
	//inline DOUBLE p

	inline SPH_REAL pRho(SPH_REAL rsqr) { return  cube(sqr(set.h) - rsqr); };
	inline SPH_REAL pAp(SPH_REAL r, SPH_REAL xi, SPH_REAL xj, SPH_REAL pi, SPH_REAL pj,
		SPH_REAL rhoi, SPH_REAL rhoj)
	{
		return set.h*(pi + pj) / 2 / rhoi / rhoj * sqr(set.h - r) * (xi - xj) / r;
	}
	inline SPH_REAL pAvis(SPH_REAL ujmui, SPH_REAL IVrhoiIVrhojXhmr, SPH_REAL udotrXcviss)
	{
		return (ujmui - udotrXcviss) * IVrhoiIVrhojXhmr ;
	}
	inline SPH_REAL pressureIG(SPH_REAL rho)
	{
		return cpressure * pow(rho, gas.gamma);
	}


	void Initialize();
	/*
	init:
	(1) get info of rank and world, set cartesian communicator, get ranks of neibours
	(2) process 0: allocate ngrids, offset for each process and broadcast,
	(3) set border grids' type, outer grids'type,  wall grids' type
	(4) get constants, box condition
	*/


	void Run(size_t maxiteration);
	/*
	run:
	(1) initial particles: not in walls, density, velocity

	(2)for iteration:
		(a).give all border particles(build a vector buffer, send, recieve with "op"s)
		(b).rebuild vector of ppointers from p, us ppointers for particle finding
		(c).cauculate r->rho->p->a >>> (vector indexing)
				bothsides,left+leftup+up+rightup, border with additions:
					adding all the near outers
			if outer is freeout, none
			if outer is pressure out, use const normal pressure
			if outer is velocity out, use const normal pressure + viscosity

		(d).renew x,v, gridx with a >>> (list iterating)
				if x is in freeout/pressureout/velocityout, delete
				if x is in exchange out, add to sendlist to corresponding process
				if x is in wall, return to previous point, vnorm = 0
		(e).send & receive exchange borders and add to list of p
		(f).spawn new particles at velocityinlet
			<in a squared particle initailizaion, hinit should be sqrt(m/rho0), dt = hinit/v0

					







	*/
};

class SPWater
{
	SmoothParticleSet set;
	WaterProperty water;
	WallBoxProperty box;
	SPH_REAL cpressure;//p = cpressure rho^gamma
	SPH_REAL cpRho;//rho = cprho sum(prho)
	SPH_REAL cpAp; //ap = cpAp sum(pAp)
	SPH_REAL cpAvis; //avis = cpAvis sum(avis)
	SPH_REAL hsqr;
	SPH_REAL cubesqrh;

	SPH_REAL cviss;
	
	SPH_REAL inletDensityX;
	SPH_REAL inletDensityY;
	size_t boxtype;
	size_t ssize;
	size_t prsize;

	SPH_REAL dt;//time step
	SPH_REAL hinit;//for spawn and initialization 

	vector<size_t> xgrids;
	vector<size_t> ygrids;
	vector<SPH_REAL> xorigs;
	vector<SPH_REAL> yorigs;
	vector<SPH_REAL> xspan;
	vector<SPH_REAL> yspan;
	int totalgridx;
	int totalgridy;


	int nproc;
	int rank;
	int rankc;//cartesian rank
	int rankx;//cartesian x, y rank	
	int ranky;

	int rankupper;
	int ranklower;
	int rankright;
	int rankleft;
	int rankupperright;
	int rankupperleft;
	int ranklowerright;
	int ranklowerleft;
	int ranks[8];
	int dims[2];
	MPI_Comm commc;
	int checker;
	int outputer;


	MPI_Request senddreqs[8];
	MPI_Request sendsreqs[8];
	MPI_Request recvdreqs[8];
	MPI_Request recvsreqs[8];


	void InitializeParticles();
	void SetGridTypes();
	void CartRankDetermine();
public:
	SPWater(WaterProperty nw, WallBoxProperty nbox, SPH_REAL nh, SPH_REAL nm)
		:water(nw), box(nbox), set(nh, nm) {
		nproc = rank = rankc = rankx = ranky = rankupper = ranklower
			= rankright = rankleft = totalgridx = totalgridy = -1;
		dt = hinit = cpAp = cpAvis = cpressure = cpRho = cviss = hsqr = 0;
			rankupperright = ranklowerleft = rankupperleft = ranklowerright= -1;
		ssize = sizeof(ParticleProperties);
		prsize = sizeof(RhoPass);
		checker = 5; outputer = 20;

	}

	//inline DOUBLE pRho(DOUBLE r);//nonconstant part of summation: rho//note that rho need to add rho0
	//inline DOUBLE pAP(DOUBLE r, DOUBLE xi, DOUBLE xj, DOUBLE pi, DOUBLE pj, DOUBLE rhoi, DOUBLE rhoj); //p acceleration
	//inline DOUBLE p

	inline SPH_REAL pRho(SPH_REAL rsqr) { return  cube(hsqr - rsqr); };
	inline SPH_REAL pAp(SPH_REAL r, SPH_REAL pi, SPH_REAL pj,
		SPH_REAL rhoi, SPH_REAL rhoj)
	{
		return (pi + pj) / 2 / rhoi / rhoj * sqr(set.h - r)   / r;//(*xi-xj)
	}
	inline SPH_REAL pAvis(SPH_REAL rhoi, SPH_REAL rhoj, SPH_REAL r)
	{
		return 1/rhoi/rhoj * (set.h-r);//(ui-uj)
	}
	inline SPH_REAL pressureWater(SPH_REAL rho)
	{
		return (rho-water.rho0)*water.K;
		//return (pow((rho/water.rho0),7.0)-1.0) * water.K;
	}
	

	void Initialize(ostream& logout);
	/*
	init:
	(1) get info of rank and world, set cartesian communicator, get ranks of neibours
	(2) process 0: allocate ngrids, offset for each process and broadcast,
	(3) set border grids' type, outer grids'type,  wall grids' type
	(4) get constants, box condition
	(5) initial particles: not in walls, density, velocity
	*/


	void Run(int maxiteration, string outpath, ostream& logout);
	/*
	run:
	

	(2)for iteration:
		(a).give all border particles(build a vector buffer, send, recieve with "op"s)
		(b).rebuild vector of ppointers from p, us ppointers for particle finding
		(c).cauculate r->rho->p->a >>> (vector indexing)
				bothsides,left+leftup+up+rightup, border with additions:
					adding all the near outers
			///if outer is freeout, none
			///if outer is pressure out, use const normal pressure
			///if outer is velocity out, use const normal pressure + viscosity

		(d).renew x,v, gridx with a >>> (list iterating)
				if x is in freeout/pressureout/velocityout, delete
				if x is in exchange out, add to sendlist to corresponding process
				if x is in wall, return to previous point, v = 0 *****

		(e).send & receive exchange borders and add to list of p
		///(f).spawn new particles at velocityinlet
			<in a squared particle initailizaion, hinit should be sqrt(m/rho0), dt = hinit/v0









	*/

	void interactRho(size_t gxa, size_t gya, size_t gxb, size_t gyb); 
	void interactPA(size_t gxa, size_t gya, size_t gxb, size_t gyb);
	void interactWall(size_t gxa, size_t gya);
	void renewP();
	void renewVXG();
	void renewVXGA();
	void reflectionWall();
	void sendreceiveOuters();
	void sendreceiveOuterRho();
	void sendreceivePassers();
	void innerInGridPointer();
	void balanceAllocate();
	void output(string path);
	
};


#undef vector
#undef string
#undef list
#undef ostream