#include "SmoothParticle.h"
#include <fstream>

#ifdef _MSC_VER
#include<easyx.h>
#include<conio.h>
#include<graphics.h>
#include<cstdio>
#endif

#ifdef __GUNC__

#endif


using namespace std;
const double pi = acos(-1);
#define SP_NOTOUCHWALL


void SPWater::Initialize(ostream& logout)
{
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	hinit = sqrt(set.m / water.rho0) * 1;
	dt = 0.5f * hinit / sqrt(water.K) * 1e0f;//k = c^2


	if (rank == 0)
	{
		logout << "h = " << set.h << ", hinit = " << hinit << endl;
	}
	dims[0] = nproc;
	dims[1] = 1;
	int periods[2] = { 0,0 };
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &commc);
	MPI_Comm_rank(commc, &rankc);
	MPI_Cart_coords(commc, rankc, 2, periods);
	rankx = periods[0], ranky = periods[1];

	totalgridx = int(ceil((box.x1 - box.x0) / set.h));
	totalgridy = int(ceil((box.y1 - box.y0) / set.h));
	xgrids.resize(dims[0]), ygrids.resize(dims[1]);
	xorigs.resize(dims[0]), yorigs.resize(dims[1]);
	xspan.resize(dims[0]), yspan.resize(dims[1]);

	CartRankDetermine();

	/// distribution 0
	for (int i = 0; i < dims[0]; i++)
	{
		if (i < totalgridx % dims[0])
			xgrids[i] = totalgridx / dims[0] + 1;
		else
			xgrids[i] = totalgridx / dims[0];
		if (i == 0)
			xorigs[i] = box.x0;
		else
			xorigs[i] = xorigs[i - 1] + xgrids[i - 1] * set.h;
		xspan[i] = xgrids[i] * set.h;
	}
	for (int i = 0; i < dims[1]; i++)
	{
		if (i < totalgridy % dims[1])
			ygrids[i] = totalgridy / dims[1] + 1;
		else
			ygrids[i] = totalgridy / dims[1];
		if (i == 0)
			yorigs[i] = box.y0;
		else
			yorigs[i] = yorigs[i - 1] + ygrids[i - 1] * set.h;
		yspan[i] = ygrids[i] * set.h;
	}
	

	set.nGridx = xgrids[rankx]+2;
	set.nGridy = ygrids[ranky]+2;
	set.iGridxm = xgrids[rankx];
	set.iGridym = ygrids[ranky];
	set.nGrid = set.nGridx * set.nGridy;
	set.xorig = xorigs[rankx];
	set.yorig = yorigs[ranky];
	set.xmax = set.xorig + (set.nGridx - 2) * set.h;
	set.ymax = set.yorig + (set.nGridy - 2) * set.h;

	set.particlesInGrid.clear();
	set.particlesInGrid.resize(set.nGrid, vector<size_t>(0));
	set.gridType.clear();
	set.gridType.resize(set.nGrid, SPGT_Common);
	set.particlesInGridpointers.clear();
	set.particlesInGridpointers.resize(set.nGrid);
	set.gridVecx.clear();
	set.gridVecy.clear();
	set.gridVecx.resize(set.nGrid, 0), set.gridVecy.resize(set.nGrid, 0);
	///



	SetGridTypes();
	InitializeParticles();

	////calculation constants
	cpRho = set.m * 4.f / pi / set.h / set.h / set.h / set.h / set.h / set.h / set.h / set.h;
	cpAp = set.m * 30.f / pi / set.h / set.h / set.h / set.h / set.h;
	cpAvis = set.m * water.mu * 40.f / pi / set.h / set.h / set.h / set.h / set.h;
	hsqr = set.h * set.h;
	cubesqrh = cube(hsqr);
	//cpressure not needed
	////_________

}

void SPWater::InitializeParticles()
{
	/// initialize particles
	size_t gx, gy;
	unsigned char gridtype;
	SPH_REAL absx = 0.5f * hinit + box.x0;
	SPH_REAL absy = 0.5f * hinit + box.y0;
	SPH_REAL xlocstart, ylocstart;
	box.x1 = totalgridx * set.h + box.x0;
	box.y1 = totalgridy * set.h + box.y0;
	SPH_REAL x, y;
	size_t particleindex = 0;
	set.p.clear();

	while (absx <= set.xorig)
		absx += hinit;
	while (absy <= set.yorig)
		absy += hinit;
	xlocstart = absx, ylocstart = absy;

	list<ParticleProperties>::iterator iter;
	SPH_REAL hinitgo;
	while (absx < set.xmax)
	{
		absy = ylocstart;
		while (absy < set.ymax)
		{
			hinitgo = hinit / sqrt(((water.rho0 * box.g * (box.y1 - absy) / water.K + water.rho0) / water.rho0));
			//hinitgo = hinit;
			if  //(sqr(absx-0.5) + sqr(absy - 0.5) > sqr(0.25))
				(absx>0.2||absy>0.4)  
				//(absx / (box.x1 - box.x0) + absy / (box.y1 - box.y0) >= 1)
			{
				absy += hinitgo;//!!!
				continue;
			}
			//x = set.localx(absx);
			//y = set.localy(absy);//we get global here
			x = absx;
			y = absy;
			gx = set.getGridx(x);
			gy = set.getGridy(y);
			gridtype = set.gridType[IBR(gy, gx, set.nGridx)];
			if ((gridtype & SPGTInnerType) == SPGT_SolidWall)
			{
				absy += hinit;
				continue;
			}
			set.p.push_back(ParticleProperties());
			iter = set.p.end();
			iter--;
			iter->x = x;
			iter->y = y;
			iter->vx = 0;
			iter->vy = 0;
			iter->gridx = gx;
			iter->gridy = gy;
			iter->rho = 0;//initializaion: for adding after
			iter->p = 0;
			set.particlesInGrid[IBR(gy, gx, set.nGridx)].push_back(particleindex);
			particleindex++;
			absy += hinitgo;
		}
		absx += hinitgo;
	}
	//set.ppointers.resize(set.p.size());

	///initialize pointer
	innerInGridPointer();

	////_________

	//// solidwall particles
	set.solidwalls.clear();
	absx = 0.5 * hinit + box.x0 - set.h;
	absy = 0.5 * hinit + box.y0 - set.h;
	while (absx <= (set.xorig - set.h))
		absx += hinit;
	while (absy <= (set.yorig - set.h))
		absy += hinit;
	xlocstart = absx, ylocstart = absy;

	while (absx < (set.xmax + set.h))
	{
		absy = ylocstart;
		while (absy < (set.ymax + set.h))
		{
			
			x = absx;
			y = absy;
			gx = set.getGridx(x);
			gy = set.getGridy(y);
			gridtype = set.gridType[IBR(gy, gx, set.nGridx)];
			if ((gridtype & SPGTInnerType) != SPGT_SolidWall)
			{
				absy += hinit;
				continue;
			}
			set.solidwalls.push_back(ParticleProperties());
			iter = set.solidwalls.end();
			iter--;
			iter->x = x;
			iter->y = y;
			iter->vx = 0;
			iter->vy = 0;
			iter->gridx = gx;
			iter->gridy = gy;
			iter->rho = water.rho0;//solidwaller
			iter->p = 0;//solidwaller
			//set.particlesInGrid[IBR(gy, gx, set.nGridx)].push_back(particleindex);
			//particleindex++;
			absy += hinit;
		}
		absx += hinit;
	}
	for (iter = set.solidwalls.begin(); iter != set.solidwalls.end(); iter++)
	{
		gx = iter->gridx, gy = iter->gridy;
		set.particlesInGridpointers[IBR(gy, gx, set.nGridx)].push_back(&(*iter));
	}

	////__________
}

void SPWater::SetGridTypes()
{
	////grid types (border wall here)
	size_t gx, gy;
	for (gx = 0; gx < set.nGridx; gx++)
		for (gy = 0; gy < set.nGridy; gy++)
		{
			if (gx == 0 || (gx == set.nGridx - 1) || gy == 0 || (gy == set.nGridy - 1))
			{
				if ((rankx == 0 && gx == 0) || (ranky == 0 && gy == 0) ||
					(rankx == dims[0] - 1 && (gx == set.nGridx - 1)) ||
					(ranky == dims[1] - 1 && (gy == set.nGridy - 1)))
				{
					set.gridType[IBR(gy, gx, set.nGridx)] = SPGT_SolidWall;
					if (gx == 0)
					{
						set.gridType[IBR(gy, gx, set.nGridx)] += SPGT_Left;
						set.gridVecx[IBR(gy, gx, set.nGridx)] = -0.99;
					}
					if (gy == 0)
					{
						set.gridType[IBR(gy, gx, set.nGridx)] += SPGT_Lower;
						set.gridVecy[IBR(gy, gx, set.nGridx)] = -0.99;
					}
					if (gx == set.nGridx - 1)
					{
						set.gridType[IBR(gy, gx, set.nGridx)] += SPGT_Right;
						set.gridVecx[IBR(gy, gx, set.nGridx)] = -0.99;
					}
					if (gy == set.nGridy - 1)
					{
						set.gridType[IBR(gy, gx, set.nGridx)] += SPGT_Upper;
						set.gridVecy[IBR(gy, gx, set.nGridx)] = -0.99;
					}
				}
				else
				{
					set.gridType[IBR(gy, gx, set.nGridx)] = SPGT_ExchangeOuter;
					if (gx == 0)
						set.gridType[IBR(gy, gx, set.nGridx)] += SPGT_Left;
					if (gy == 0)
						set.gridType[IBR(gy, gx, set.nGridx)] += SPGT_Lower;
					if (gx == set.nGridx - 1)
						set.gridType[IBR(gy, gx, set.nGridx)] += SPGT_Right;
					if (gy == set.nGridy - 1)
						set.gridType[IBR(gy, gx, set.nGridx)] += SPGT_Upper;
				}
			}
			else if (gx == 1 || (gx == set.nGridx - 2) || gy == 1 || (gy == set.nGridy - 2))
			{
				set.gridType[IBR(gy, gx, set.nGridx)] = SPGT_Border;
				if (gx == 1)
					set.gridType[IBR(gy, gx, set.nGridx)] += SPGT_Left;
				if (gy == 1)
					set.gridType[IBR(gy, gx, set.nGridx)] += SPGT_Lower;
				if (gx == set.nGridx - 2)
					set.gridType[IBR(gy, gx, set.nGridx)] += SPGT_Right;
				if (gy == set.nGridy - 2)
					set.gridType[IBR(gy, gx, set.nGridx)] += SPGT_Upper;
			}
		}
	////_________
}

void SPWater::CartRankDetermine()
{
	int query[2];
	if (rankx > 0)
	{
		query[0] = rankx - 1, query[1] = ranky;
		MPI_Cart_rank(commc, query, &rankleft);
	}
	else rankleft = -1;
	if (rankx < dims[0] - 1)
	{
		query[0] = rankx + 1, query[1] = ranky;
		MPI_Cart_rank(commc, query, &rankright);
	}
	else rankright = -1;
	if (ranky > 0)
	{
		query[0] = rankx, query[1] = ranky - 1;
		MPI_Cart_rank(commc, query, &ranklower);
	}
	else ranklower = -1;
	if (ranky < dims[1] - 1)
	{
		query[0] = rankx, query[1] = ranky + 1;
		MPI_Cart_rank(commc, query, &rankupper);
	}
	else rankupper = -1;

	if (rankx > 0 && ranky > 0)
	{
		query[0] = rankx - 1, query[1] = ranky - 1;
		MPI_Cart_rank(commc, query, &ranklowerleft);
	}
	else ranklowerleft = -1;
	if ((rankx < dims[0] - 1) && ranky > 0)
	{
		query[0] = rankx + 1, query[1] = ranky - 1;
		MPI_Cart_rank(commc, query, &ranklowerright);
	}
	else ranklowerright = -1;
	if (rankx > 0 && (ranky < dims[1] - 1))
	{
		query[0] = rankx - 1, query[1] = ranky + 1;
		MPI_Cart_rank(commc, query, &rankupperleft);
	}
	else rankupperleft = -1;
	if ((rankx < dims[0] - 1) && (ranky < dims[1] - 1))
	{
		query[0] = rankx + 1, query[1] = ranky + 1;
		MPI_Cart_rank(commc, query, &rankupperright);
	}
	else rankupperright = -1;
	ranks[0] = rankupper;
	ranks[1] = ranklower;
	ranks[2] = rankleft;
	ranks[3] = rankright;
	ranks[4] = rankupperleft;
	ranks[5] = rankupperright;
	ranks[6] = ranklowerleft;
	ranks[7] = ranklowerright;
}

void SPWater::Run(int maxiteration, string outpath, ostream& logout)
{


	size_t gx, gy;
	unsigned char gridtype;
	char namebuf[32];
	
	

	
	

#ifdef _MSC_VER
	int winy = 800;
	int winx = int(ceil(winy / (set.ymax - set.yorig) * (set.xmax - set.xorig)));
	int dx, dy;
	double dwinx, dwiny;
	double vmax = 0.01;
	int r = int(set.h / (set.ymax - set.yorig) * winy);
	r = 1;
	TCHAR timetext[100];
	initgraph(winx, winy, EW_SHOWCONSOLE);
	setfillcolor(0xCCCCCC);
	//system("pause");
	RECT textrect = { 0,0,200,100 };
#endif
	for (int it = 1; it <= maxiteration; it++)
	{
		///Send border and receive eouter and renew pointers
		sendreceiveOuters();
#ifndef SP_NOTOUCHWALL
		reflectionWall();
#endif

		///interactrho /////+ rwallversion£¿
		for(gy = 1; gy<=set.iGridym; gy++)
			for (gx = 1; gx <= set.iGridxm; gx++)
			{
				gridtype = set.gridType[IBR(gy, gx, set.nGridx)];
				if ((gridtype & SPGTInnerType) == SPGT_Border)
				{
					switch (gridtype & SPGTCommType)
					{
					case SPGT_Lower | SPGT_Right:
						interactRho(gx, gy, gx + 1, gy); 
					case SPGT_Lower:
					case SPGT_Lower | SPGT_Left:
						interactRho(gx, gy, gx, gy - 1);
						interactRho(gx, gy, gx + 1, gy - 1);
					case SPGT_Upper | SPGT_Left:
					case SPGT_Left:
						interactRho(gx, gy, gx - 1, gy - 1);
					case SPGT_Upper:
						interactRho(gx, gy, gx - 1, gy);
						interactRho(gx, gy, gx - 1, gy + 1);
						interactRho(gx, gy, gx, gy + 1);
						interactRho(gx, gy, gx + 1, gy + 1);
						break;

					case SPGT_Right:
					case SPGT_Upper | SPGT_Right:
						interactRho(gx, gy, gx + 1, gy - 1);
						interactRho(gx, gy, gx + 1, gy);

						interactRho(gx, gy, gx - 1, gy);
						interactRho(gx, gy, gx - 1, gy + 1);
						interactRho(gx, gy, gx, gy + 1);
						interactRho(gx, gy, gx + 1, gy + 1);
						break;
					}
				}
				else
				{
					interactRho(gx, gy, gx - 1, gy);
					interactRho(gx, gy, gx - 1, gy + 1);
					interactRho(gx, gy, gx, gy + 1);
					interactRho(gx, gy, gx + 1, gy + 1);
				}
				interactRho(gx, gy, gx, gy);
				
			}
		
		///renewp 
		renewP();

		///Send receive eouter rho!!!
		sendreceiveOuterRho();

		///interactPA /////+ rwallversion ?
		for (gy = 1; gy <= set.iGridym; gy++)
		{
			for (gx = 1; gx <= set.iGridxm; gx++)
			{
				gridtype = set.gridType[IBR(gy, gx, set.nGridx)];
				if ((gridtype & SPGTInnerType) == SPGT_Border)
				{
					
					switch (gridtype & SPGTCommType)
					{
					case SPGT_Lower | SPGT_Right:
						interactPA(gx, gy, gx + 1, gy);
					case SPGT_Lower:
					case SPGT_Lower | SPGT_Left:
						interactPA(gx, gy, gx, gy - 1);
						interactPA(gx, gy, gx + 1, gy - 1);
					case SPGT_Upper | SPGT_Left:
					case SPGT_Left:
						interactPA(gx, gy, gx - 1, gy - 1);
					case SPGT_Upper:
						interactPA(gx, gy, gx - 1, gy);
						interactPA(gx, gy, gx - 1, gy + 1);
						interactPA(gx, gy, gx, gy + 1);
						interactPA(gx, gy, gx + 1, gy + 1);
						break;

					case SPGT_Right:
					case SPGT_Upper | SPGT_Right:
						interactPA(gx, gy, gx + 1, gy - 1);
						interactPA(gx, gy, gx + 1, gy);

						interactPA(gx, gy, gx - 1, gy);
						interactPA(gx, gy, gx - 1, gy + 1);
						interactPA(gx, gy, gx, gy + 1);
						interactPA(gx, gy, gx + 1, gy + 1);
						break;
					}
					
					//interactWall(gx, gy);
					interactPA(gx, gy, gx, gy);
				}
				else
				{
					interactPA(gx, gy, gx - 1, gy);
					interactPA(gx, gy, gx - 1, gy + 1);
					interactPA(gx, gy, gx, gy + 1);
					interactPA(gx, gy, gx + 1, gy + 1);
					interactPA(gx, gy, gx, gy);
				}
				
			}
		}

		///add gravity
		///renew v x grid return walls, 
		if (it > 1)
			renewVXG();
		else
			renewVXGA();

		///send receive inpoints and renew pointers
		sendreceivePassers();
		innerInGridPointer();
		
		if (it % checker == 0)
		{
			///check numbers
			///reallocate
			balanceAllocate();
			logout << "Iteration " << it << endl;
		}
		if (it % outputer == 0)
		{
#ifdef _MSC_VER
			///output
			sprintf_s(namebuf, "%.3f.txt", dt* it);
			output(outpath + namebuf);
			cleardevice();
			BeginBatchDraw();
			list<ParticleProperties>::iterator iter;
			setfillcolor(0xCCCCCC);
			setlinecolor(0xFFFFFF);
			sprintf_s(timetext, "t = %f",double(dt*it));
			drawtext(timetext,&textrect , DT_LEFT | DT_TOP);
			int i = 0;
			for (iter = set.p.begin(); iter != set.p.end(); iter++)
			{
				dwinx = iter->x - set.xorig;
				dwiny = iter->y - set.yorig;
				dx = int(dwinx /( set.xmax - set.xorig) * winx);
				dy = winy - int(dwiny / (set.ymax - set.yorig) * winy);
				//if (i % 100 == 0)
				//	setfillcolor(0xAA9900);
				setfillcolor(HSVtoRGB(sqrt(sqr(iter->vx) + sqr(iter->vy))/vmax,1,0.5));
				solidcircle(dx, dy, r);
				//setfillcolor(0xCCCCCC);
				i++;
			}
			EndBatchDraw();
#endif
		}
	}
}

void SPWater::interactRho(size_t gxa, size_t gya, size_t gxb, size_t gyb)
{
	size_t gia = IBR(gya, gxa, set.nGridx), gib = IBR(gyb, gxb, set.nGridx);
	ParticleProperties* pa, *pb;
	SPH_REAL rsqr, rhoadd;
	if ((set.gridType[gib] & SPGTInnerType) == SPGT_SolidWall)
		for (size_t i = 0; i < set.particlesInGridpointers[gia].size(); i++)
			for (size_t j = 0; j < set.particlesInGridpointers[gib].size(); j++)
			{
				pa = set.particlesInGridpointers[gia][i];
				pb = set.particlesInGridpointers[gib][j];
				rsqr = sqr(pa->x - pb->x) + sqr(pa->y - pb->y);
				if (rsqr < hsqr)
				{
					rhoadd = pRho(rsqr);
					pa->rho += rhoadd;
#ifndef SP_NOTOUCHWALL
					pb->rho += rhoadd;
#endif
				}
			}
	else if(gia - gib)
		for(size_t i = 0; i < set.particlesInGridpointers[gia].size(); i++)
			for (size_t j = 0; j < set.particlesInGridpointers[gib].size(); j++)
			{
				pa = set.particlesInGridpointers[gia][i];
				pb = set.particlesInGridpointers[gib][j];
				rsqr = sqr(pa->x - pb->x) + sqr(pa->y - pb->y);
				if (rsqr < hsqr)
				{
					rhoadd = pRho(rsqr);
					pa->rho += rhoadd;
					pb->rho += rhoadd;
				}
			}	
	else
		for(size_t i = 0; i<set.particlesInGridpointers[gia].size(); i++)
			for (size_t j = i+1; j < set.particlesInGridpointers[gib].size(); j++)
			{
				pa = set.particlesInGridpointers[gia][i];
				pb = set.particlesInGridpointers[gib][j];
				rsqr = sqr(pa->x - pb->x) + sqr(pa->y - pb->y);
				//if (rsqr < hsqr)
				//{
				rhoadd = pRho(rsqr);
				pa->rho += rhoadd;
				pb->rho += rhoadd;
				//}
			}
	/////////// cpRho is put in in renewP();
}

void SPWater::interactPA(size_t gxa, size_t gya, size_t gxb, size_t gyb)
{
	size_t gia = IBR(gya, gxa, set.nGridx), gib = IBR(gyb, gxb, set.nGridx);
	ParticleProperties* pa, * pb;
	SPH_REAL r,aaddp,aaddvis, axadd, ayadd;
	if((set.gridType[gib] & SPGTInnerType) == SPGT_SolidWall)
		for (size_t i = 0; i < set.particlesInGridpointers[gia].size(); i++)
			for (size_t j = 0; j < set.particlesInGridpointers[gib].size(); j++)
			{
				pa = set.particlesInGridpointers[gia][i];
				pb = set.particlesInGridpointers[gib][j];
				r = sqrt(sqr(pa->x - pb->x) + sqr(pa->y - pb->y));
				if (r < set.h)
				{
					aaddp = cpAp * pAp(r, pa->p, pb->p, pa->rho, pb->rho);
					aaddvis = cpAvis * pAvis(pa->rho, pb->rho, r);
					axadd = aaddp * (pa->x - pb->x) + aaddvis * (pb->vx - pa->vx);
					ayadd = aaddp * (pa->y - pb->y) + aaddvis * (pb->vy - pa->vy);
					pa->ax += axadd;
					pa->ay += ayadd;
#ifndef SP_NOTOUCHWALL
					pb->ax -= axadd;
					pb->ay -= ayadd;
#endif
				}
			}
	else if (gia - gib)
		for (size_t i = 0; i < set.particlesInGridpointers[gia].size(); i++)
			for (size_t j = 0; j < set.particlesInGridpointers[gib].size(); j++)
			{
				pa = set.particlesInGridpointers[gia][i];
				pb = set.particlesInGridpointers[gib][j];
				r = sqrt(sqr(pa->x - pb->x) + sqr(pa->y - pb->y));
				if (r < set.h)
				{
					aaddp = cpAp * pAp(r, pa->p, pb->p, pa->rho, pb->rho);
					aaddvis =  cpAvis * pAvis(pa->rho, pb->rho, r);
					axadd = aaddp * (pa->x - pb->x) + aaddvis * (pb->vx - pa->vx);
					ayadd = aaddp * (pa->y - pb->y) + aaddvis * (pb->vy - pa->vy);
					pa->ax += axadd;
					pa->ay += ayadd;
					pb->ax -= axadd;
					pb->ay -= ayadd;
				}
			}
	else
		for (size_t i = 0; i < set.particlesInGridpointers[gia].size(); i++)
			for (size_t j = i + 1; j < set.particlesInGridpointers[gib].size(); j++)
			{
				pa = set.particlesInGridpointers[gia][i];
				pb = set.particlesInGridpointers[gib][j];
				r = sqrt(sqr(pa->x - pb->x) + sqr(pa->y - pb->y));
				
				aaddp = cpAp * pAp(r, pa->p, pb->p, pa->rho, pb->rho);
				aaddvis = cpAvis * pAvis(pa->rho, pb->rho, r);
				axadd = aaddp * (pa->x - pb->x) + aaddvis * (pb->vx - pa->vx);
				ayadd = aaddp * (pa->y - pb->y) + aaddvis * (pb->vy - pa->vy);
				pa->ax += axadd;
				pa->ay += ayadd;
				pb->ax -= axadd;
				pb->ay -= ayadd;
			}
	

}

void SPWater::interactWall(size_t gxa, size_t gya)
{
	size_t gia = IBR(gya, gxa, set.nGridx);
	ParticleProperties* pa;
	for (size_t i = 0; i < set.particlesInGridpointers[gia].size(); i++)
	{
		pa = set.particlesInGridpointers[gia][i];
		pa->ax = pa->ay = pa->vx = pa->vy = 0;
	}
			



}

void SPWater::renewP()
{
	/*
	for (INT i = 0; i < set.ppointers.size(); i++)
	{
		set.ppointers[i]->rho = (set.ppointers[i]->rho + cubesqrh) * cpRho;
		set.ppointers[i]->p = pressureWater(set.ppointers[i]->rho);
	}
	*/
	list<ParticleProperties>::iterator iter;
	for(iter = set.p.begin(); iter!=set.p.end(); iter++)
	{
		iter->rho = (iter->rho + cubesqrh) * cpRho;
		iter->p = pressureWater(iter->rho);
	}
#ifndef SP_NOTOUCHWALL
	for(iter = set.solidwalls.begin(); iter!=set.solidwalls.end(); iter++)
		{
			iter->rho = (iter->rho + cubesqrh) * cpRho;
			iter->p = pressureWater(iter->rho);
		}
#endif
}



void SPWater::renewVXG()
{
	list<ParticleProperties>::iterator iter;
	unsigned char gridtype;
	
	for (iter = set.p.begin(); iter != set.p.end(); )
	{
		iter->vx += iter->ax * dt;
		iter->vy += (iter->ay - box.g) * dt;
		iter->x += iter->vx * dt;
		iter->y += iter->vy * dt;
		iter->rho = 0;
		iter->ax = iter->ay = 0;
		iter->gridx = set.getGridx(iter->x);
		iter->gridy = set.getGridy(iter->y);
		/*
		if (iter->gridx < 0)
			iter->x -= iter->vx*dt, iter->gridx = 0, iter->vx = 0;
		if (iter->gridx >= set.nGridx)
			iter->x -= iter->vx * dt , iter->gridx = set.nGridx - 1, iter->vx = 0;
		if (iter->gridy < 0)
			iter->y -= iter->vy *dt, iter->gridy = 0, iter->vy = 0;
		if (iter->gridy >= set.nGridy)
			iter->y -= iter->vy * dt, iter->gridy = set.nGridy - 1, iter->vy = 0;
			*/
		if ((iter->gridx < 0) || (iter->gridx >= set.nGridx) || (iter->gridy < 0) || (iter->gridy >= set.nGridy))
		{
			iter = set.p.erase(iter);
			continue;
		}
		gridtype = set.gridType[IBR(iter->gridy, iter->gridx, set.nGridx)];
		if ((SPGTInnerType & gridtype) == SPGT_SolidWall)//solidwall PAGAGAPA
		{
			iter->x -= iter->vx * dt;
			iter->y -= iter->vy * dt;
			if (gridtype & SPGT_Left || gridtype & SPGT_Right)
			{
				//iter->vx *= -1;
				iter->vx *= set.gridVecx[IBR(iter->gridy, iter->gridx, set.nGridx)];
			}
			if (gridtype & SPGT_Upper || gridtype & SPGT_Lower)
			{
				//iter->vy *= -1;
				iter->vy *= set.gridVecy[IBR(iter->gridy, iter->gridx, set.nGridx)];
			}
			iter->gridx = set.getGridx(iter->x);
			iter->gridy = set.getGridy(iter->y);
			
			//iter->vx = 0, iter->vy = 0;
		}
		iter++;
	}
}// also deals with walls

void SPWater::renewVXGA()
{
	list<ParticleProperties>::iterator iter;
	unsigned char gridtype;

	for (iter = set.p.begin(); iter != set.p.end(); )
	{
		iter->vx += iter->ax * dt;
		iter->vy += (iter->ay - box.g) * dt;
		iter->x += iter->vx * dt/2;
		iter->y += iter->vy * dt/2;
		iter->rho = 0;
		iter->ax = iter->ay = 0;
		iter->gridx = set.getGridx(iter->x);
		iter->gridy = set.getGridy(iter->y);
		/*
		if (iter->gridx < 0)
			iter->x -= iter->vx*dt, iter->gridx = 0, iter->vx = 0;
		if (iter->gridx >= set.nGridx)
			iter->x -= iter->vx * dt , iter->gridx = set.nGridx - 1, iter->vx = 0;
		if (iter->gridy < 0)
			iter->y -= iter->vy *dt, iter->gridy = 0, iter->vy = 0;
		if (iter->gridy >= set.nGridy)
			iter->y -= iter->vy * dt, iter->gridy = set.nGridy - 1, iter->vy = 0;
			*/
		if ((iter->gridx < 0) || (iter->gridx >= set.nGridx) || (iter->gridy < 0) || (iter->gridy >= set.nGridy))
		{
			iter = set.p.erase(iter);
			continue;
		}
		gridtype = set.gridType[IBR(iter->gridy, iter->gridx, set.nGridx)];
		if ((SPGTInnerType & gridtype) == SPGT_SolidWall)//solidwall PAGAGAPA
		{
			iter->x -= iter->vx * dt;
			iter->y -= iter->vy * dt;
			if (gridtype & SPGT_Left || gridtype & SPGT_Right)
				iter->vx *= -1;//set.gridVecx[IBR(iter->gridy, iter->gridx, set.nGridx)];
			if (gridtype & SPGT_Upper || gridtype & SPGT_Lower)
				iter->vy *= -1;//set.gridVecy[IBR(iter->gridy, iter->gridx, set.nGridx)];
			iter->gridx = set.getGridx(iter->x);
			iter->gridy = set.getGridy(iter->y);

			//iter->vx = 0, iter->vy = 0;
		}
		iter++;
	}
}// also deals with walls

void SPWater::reflectionWall()
{
	unsigned char gridtype;
	size_t index, getindex;
	set.solidwalls.clear();
	list<ParticleProperties>::iterator iter;
	for (size_t gy = 0; gy < set.nGridy; gy++)
	{
		for (size_t gx = 0; gx < set.nGridx; gx++)
		{
			index = IBR(gy, gx, set.nGridx);
			gridtype = set.gridType[index];
			if (((gridtype & SPGTInnerType) == SPGT_SolidWall))
			{
				set.particlesInGridpointers[index].clear();
				switch (gridtype & SPGTCommType)
				{
				case SPGT_Right:
					getindex = IBR(gy, gx - 1, set.nGridx);
					for (size_t i = 0; i < set.particlesInGridpointers[getindex].size();i++)
					{
						set.solidwalls.push_back(*set.particlesInGridpointers[getindex][i]);
						iter = set.solidwalls.end();
						iter--;
						set.particlesInGridpointers[index].push_back(&(*iter));
						iter->gridx = gx;
						iter->x = 2*set.xmax -iter->x;
						iter->vx = -iter->vx;
					}
					break;
				case SPGT_Left:
					getindex = IBR(gy, gx + 1, set.nGridx);
					for (size_t i = 0; i < set.particlesInGridpointers[getindex].size(); i++)
					{
						set.solidwalls.push_back(*set.particlesInGridpointers[getindex][i]);
						iter = set.solidwalls.end();
						iter--;
						set.particlesInGridpointers[index].push_back(&(*iter));
						iter->gridx = gx;
						iter->x = 2 * set.xorig -iter->x;
						iter->vx = -iter->vx;
					}
					break;
				case SPGT_Upper:
					getindex = IBR(gy - 1, gx, set.nGridx);
					for (size_t i = 0; i < set.particlesInGridpointers[getindex].size(); i++)
					{
						set.solidwalls.push_back(*set.particlesInGridpointers[getindex][i]);
						iter = set.solidwalls.end();
						iter--;
						set.particlesInGridpointers[index].push_back(&(*iter));
						iter->gridy = gy;
						iter->y = 2 * set.ymax -iter->y;
						iter->vy = -iter->vy;
					}
					break;
				case SPGT_Lower:
					getindex = IBR(gy + 1, gx, set.nGridx);
					for (size_t i = 0; i < set.particlesInGridpointers[getindex].size(); i++)
					{
						set.solidwalls.push_back(*set.particlesInGridpointers[getindex][i]);
						iter = set.solidwalls.end();
						iter--;
						set.particlesInGridpointers[index].push_back(&(*iter));
						iter->gridy = gy;
						iter->y = 2 * set.yorig -iter->y;
						iter->vy = -iter->vy;
					}
					break;
				case SPGT_Upper | SPGT_Right:
					getindex = IBR(gy - 1, gx - 1, set.nGridx);
					for (size_t i = 0; i < set.particlesInGridpointers[getindex].size(); i++)
					{
						set.solidwalls.push_back(*set.particlesInGridpointers[getindex][i]);
						iter = set.solidwalls.end();
						iter--;
						set.particlesInGridpointers[index].push_back(&(*iter));
						iter->gridx = gx;
						iter->x = 2*set.xmax -iter->x;
						iter->vx = -iter->vx;
						iter->gridy = gy;
						iter->y = 2*set.ymax -iter->y;
						iter->vy = -iter->vy;
					}
					break;
				case SPGT_Lower | SPGT_Right:
					getindex = IBR(gy + 1, gx - 1, set.nGridx);
					for (size_t i = 0; i < set.particlesInGridpointers[getindex].size(); i++)
					{
						set.solidwalls.push_back(*set.particlesInGridpointers[getindex][i]);
						iter = set.solidwalls.end();
						iter--;
						set.particlesInGridpointers[index].push_back(&(*iter));
						iter->gridx = gx;
						iter->x = 2 * set.xmax - iter->x;
						iter->vx = -iter->vx;
						iter->gridy = gy;
						iter->y = 2 * set.yorig - iter->y;
						iter->vy = -iter->vy;
					}
					break;
				case SPGT_Upper | SPGT_Left:
					getindex = IBR(gy - 1, gx + 1, set.nGridx);
					for (size_t i = 0; i < set.particlesInGridpointers[getindex].size(); i++)
					{
						set.solidwalls.push_back(*set.particlesInGridpointers[getindex][i]);
						iter = set.solidwalls.end();
						iter--;
						set.particlesInGridpointers[index].push_back(&(*iter));
						iter->gridx = gx;
						iter->x = 2 * set.xorig - iter->x;
						iter->vx = -iter->vx;
						iter->gridy = gy;
						iter->y = 2 * set.ymax - iter->y;
						iter->vy = -iter->vy;
					}
					break;
				case SPGT_Lower | SPGT_Left:
					getindex = IBR(gy + 1, gx + 1, set.nGridx);
					for (size_t i = 0; i < set.particlesInGridpointers[getindex].size(); i++)
					{
						set.solidwalls.push_back(*set.particlesInGridpointers[getindex][i]);
						iter = set.solidwalls.end();
						iter--;
						set.particlesInGridpointers[index].push_back(&(*iter));
						iter->gridx = gx;
						iter->x = 2 * set.xorig - iter->x;
						iter->vx = -iter->vx;
						iter->gridy = gy;
						iter->y = 2 * set.yorig - iter->y;
						iter->vy = -iter->vy;
					}
					break;

				default:
					break;
				}
			}
		}
	}
}

void SPWater::sendreceiveOuters()
{
	vector<vector<ParticleProperties>> sendto(8);
	int sendsize[8];
	int recvsize[8];
	//vector<ParticleProperties> recvfrom;
	list<ParticleProperties> sendlist;
	MPI_Status stats[8];
	//unsigned char gridtype;
	
	list<ParticleProperties>::iterator iter;

	for (int r = 0; r < 8; r++)
	{
		if (ranks[r] >= 0)
		{
			sendlist.clear();
			switch (r)
			{
			case 0: //upper
				for (size_t i = 1; i <= set.iGridxm; i++)
				{
					for (size_t j = 0; j < set.particlesInGridpointers[IBR(set.iGridym,i,set.nGridx)].size(); j++)
					{
						sendlist.push_back(*set.particlesInGridpointers[IBR(set.iGridym, i, set.nGridx)][j]);
						iter = sendlist.end(), iter--;
						iter->gridy = 0;
						//iter->y -= yspan[ranky];
					}
				}
				break; 
			case 1: //lower
				for (size_t i = 1; i <= set.iGridxm; i++)
				{
					for (size_t j = 0; j < set.particlesInGridpointers[IBR(1, i, set.nGridx)].size(); j++)
					{
						sendlist.push_back(*set.particlesInGridpointers[IBR(1, i, set.nGridx)][j]);
						iter = sendlist.end(), iter--;
						iter->gridy = ygrids[ranky-1ll] + 1;
						//iter->y += yspan[ranky - 1];
					}
				}
				break;
			case 2: //left
				for (size_t i = 1; i <= set.iGridym; i++)
				{
					for (size_t j = 0; j < set.particlesInGridpointers[IBR(i, 1, set.nGridx)].size(); j++)
					{
						sendlist.push_back(*set.particlesInGridpointers[IBR(i, 1, set.nGridx)][j]);
						iter = sendlist.end(), iter--;
						iter->gridx = xgrids[rankx-1ll] + 1;
					}
				}

				break;
			case 3: //right
				for (size_t i = 1; i <= set.iGridym; i++)
				{
					for (size_t j = 0; j < set.particlesInGridpointers[IBR(i, set.iGridxm, set.nGridx)].size(); j++)
					{
						sendlist.push_back(*set.particlesInGridpointers[IBR(i, set.iGridxm, set.nGridx)][j]);
						iter = sendlist.end(), iter--;
						iter->gridx = 0;
					}
				}
				break;
			case 4: //ul
				//for (INT i = 1; i < set.nGridy - 1; i++)
				//{
				for (size_t j = 0; j < set.particlesInGridpointers[IBR(set.iGridym, 1, set.nGridx)].size(); j++)
				{
					sendlist.push_back(*set.particlesInGridpointers[IBR(set.iGridym, 1, set.nGridx)][j]);
					iter = sendlist.end(), iter--;
					iter->gridx = xgrids[rankx-1ll]+1, iter->gridy = 0;
				}
				//}
				break;
			case 5: //ur
				//for (INT i = 1; i < set.nGridy - 1; i++)
				//{
				for (size_t j = 0; j < set.particlesInGridpointers[IBR(set.iGridym, set.iGridxm, set.nGridx)].size(); j++)
				{
					sendlist.push_back(*set.particlesInGridpointers[IBR(set.iGridym, set.iGridxm, set.nGridx)][j]);
					iter = sendlist.end(), iter--;
					iter->gridx = 0, iter->gridy = 0;
				}
				//}
				break;
			case 6: //ll
				//for (INT i = 1; i < set.nGridy - 1; i++)
				//{
				for (size_t j = 0; j < set.particlesInGridpointers[IBR(1, 1, set.nGridx)].size(); j++)
				{
					sendlist.push_back(*set.particlesInGridpointers[IBR(1, 1, set.nGridx)][j]);
					iter = sendlist.end(), iter--;
					iter->gridx = xgrids[rankx - 1ll] + 1, iter->gridy = ygrids[ranky - 1ll] + 1;
				}
				//}

				break;
			case 7: //lr
				//for (INT i = 1; i < set.nGridy - 1; i++)
				//{
				for (size_t j = 0; j < set.particlesInGridpointers[IBR(1, set.iGridxm, set.nGridx)].size(); j++)
				{
					sendlist.push_back(*set.particlesInGridpointers[IBR(1, set.iGridxm, set.nGridx)][j]);
					iter = sendlist.end(), iter--;
					iter->gridx = 0, iter->gridy = ygrids[ranky - 1ll] + 1;
				}
				//}
				break;
			default:
				break;
			}
			
			sendto[r].resize(sendlist.size());
			sendsize[r] = sendlist.size();
			size_t i = 0;
			for (iter = sendlist.begin(); iter != sendlist.end(); iter++, i++)
			{
				sendto[r][i] = *iter;
			}
			MPI_Isend(sendsize + r, 1, SPH_MPI_INT, ranks[r], 101, commc, sendsreqs + r);
			if (sendlist.size() > 0)
				MPI_Isend(&sendto[r][0], sendlist.size() * ssize, MPI_CHAR, ranks[r], 100, commc, senddreqs + r);
			else
				senddreqs[r] = MPI_REQUEST_NULL;
		}
		else
		{
			senddreqs[r] = sendsreqs[r] = MPI_REQUEST_NULL;	
		}
	}


	for (int r = 0; r < 8; r++)
	{
		if (ranks[r] >= 0)
		{
			set.outerones[r].clear();
			MPI_Irecv(recvsize + r, 1, SPH_MPI_INT, ranks[r], 101, commc, recvsreqs + r);
			MPI_Wait(recvsreqs + r, stats + r);
			set.outerones[r].resize(recvsize[r]);  
			if (recvsize[r] > 0)
				MPI_Irecv(&set.outerones[r][0], recvsize[r] * ssize, MPI_CHAR, ranks[r], 100, commc, recvdreqs + r);
			else
				recvdreqs[r] = MPI_REQUEST_NULL;
		}
		else
		{
			set.outerones[r].clear();
			recvdreqs[r] = MPI_REQUEST_NULL;
		}
	}
	MPI_Waitall(8, sendsreqs, stats);
	MPI_Waitall(8, senddreqs, stats);
	MPI_Waitall(8, recvdreqs, stats);
	

	for (size_t i = 0; i < set.nGridx; i++)
	{
		if((set.gridType[IBR(0, i, set.nGridx)] & SPGTInnerType) == SPGT_ExchangeOuter)
			set.particlesInGridpointers[IBR(0, i, set.nGridx)].clear();
		if ((set.gridType[IBR(set.nGridy - 1, i, set.nGridx)] & SPGTInnerType) == SPGT_ExchangeOuter)
			set.particlesInGridpointers[IBR(set.nGridy-1, i, set.nGridx)].clear();
	}
	for (size_t i = 0; i < set.nGridy; i++)
	{
		if ((set.gridType[IBR(i, 0, set.nGridx)] & SPGTInnerType) == SPGT_ExchangeOuter)
			set.particlesInGridpointers[IBR(i, 0, set.nGridx)].clear();
		if ((set.gridType[IBR(i, set.nGridx - 1, set.nGridx)] & SPGTInnerType) == SPGT_ExchangeOuter)
			set.particlesInGridpointers[IBR(i, set.nGridx - 1, set.nGridx)].clear();
	}

	size_t gx, gy;
	for (int r = 0; r < 8; r++)
	{
		if (ranks[r] >= 0)
		{
			for (size_t i = 0; i < set.outerones[r].size(); i++)
			{
				gx = set.outerones[r][i].gridx;
				gy = set.outerones[r][i].gridy;
				set.particlesInGridpointers[IBR(gy, gx, set.nGridx)].push_back(&set.outerones[r][i]);
			}
		}
	}
	//MPI_Barrier(commc);
}// send border and receive outer ones, and reorganize those pointers 

/*
void SPWater::sendreceiveOuterRho()
{
	vector<vector<RhoPass>> sendto(8);
	vector<vector<RhoPass>> recvfrom(8);
	int sendsize[8];
	int recvsize[8];
	//vector<ParticleProperties> recvfrom;
	list<RhoPass> sendlist;
	MPI_Status stats[8];
	MPI_Request senddreqs[8];
	MPI_Request sendsreqs[8];
	MPI_Request recvdreqs[8];
	MPI_Request recvsreqs[8];
	list<RhoPass>::iterator iter;
	RhoPass passbuf;

	for (int r = 0; r < 8; r++)
	{
		if (ranks[r] >= 0)
		{
			sendlist.clear();
			switch (r)
			{
			case 0: //upper
				for (INT i = 1; i <= set.iGridxm; i++)
				{
					for (INT j = 0; j < set.particlesInGridpointers[IBR(set.iGridym, i, set.nGridx)].size(); j++)
					{
						passbuf = *set.particlesInGridpointers[IBR(set.iGridym, i, set.nGridx)][j];
						sendlist.push_back(passbuf);
						iter = sendlist.end(), iter--;
						iter->gridy = 0;
					}
				}
				break;
			case 1: //lower
				for (INT i = 1; i <= set.iGridxm; i++)
				{
					for (INT j = 0; j < set.particlesInGridpointers[IBR(1, i, set.nGridx)].size(); j++)
					{
						passbuf = *set.particlesInGridpointers[IBR(1, i, set.nGridx)][j];
						sendlist.push_back(passbuf);
						iter = sendlist.end(), iter--;
						iter->gridy = ygrids[ranky - 1] + 1;
					}
				}
				break;
			case 2: //left
				for (INT i = 1; i <= set.iGridym; i++)
				{
					for (INT j = 0; j < set.particlesInGridpointers[IBR(i, 1, set.nGridx)].size(); j++)
					{
						passbuf = *set.particlesInGridpointers[IBR(i, 1, set.nGridx)][j];
						sendlist.push_back(passbuf);
						iter = sendlist.end(), iter--;
						iter->gridx = xgrids[rankx - 1] + 1;
					}
				}

				break;
			case 3: //right
				for (INT i = 1; i <= set.iGridym; i++)
				{
					for (INT j = 0; j < set.particlesInGridpointers[IBR(i, set.iGridxm, set.nGridx)].size(); j++)
					{
						passbuf = *set.particlesInGridpointers[IBR(i, set.iGridxm, set.nGridx)][j];
						sendlist.push_back(passbuf);
						iter = sendlist.end(), iter--;
						iter->gridx = 0;
					}
				}
				break;
			case 4: //ul
				//for (INT i = 1; i < set.nGridy - 1; i++)
				//{
				for (INT j = 0; j < set.particlesInGridpointers[IBR(set.iGridym, 1, set.nGridx)].size(); j++)
				{
					passbuf = *set.particlesInGridpointers[IBR(set.iGridym, 1, set.nGridx)][j];
					sendlist.push_back(passbuf);
					iter = sendlist.end(), iter--;
					iter->gridx = xgrids[rankx - 1] + 1, iter->gridy = 0;
				}
				//}
				break;
			case 5: //ur
				//for (INT i = 1; i < set.nGridy - 1; i++)
				//{
				for (INT j = 0; j < set.particlesInGridpointers[IBR(set.iGridym, set.iGridxm, set.nGridx)].size(); j++)
				{
					passbuf = *set.particlesInGridpointers[IBR(set.iGridym, set.iGridxm, set.nGridx)][j];
					sendlist.push_back(passbuf);
					iter = sendlist.end(), iter--;
					iter->gridx = 0, iter->gridy = 0;
				}
				//}
				break;
			case 6: //ll
				//for (INT i = 1; i < set.nGridy - 1; i++)
				//{
				for (INT j = 0; j < set.particlesInGridpointers[IBR(1, 1, set.nGridx)].size(); j++)
				{
					passbuf = *set.particlesInGridpointers[IBR(1, 1, set.nGridx)][j];
					sendlist.push_back(passbuf);
					iter = sendlist.end(), iter--;
					iter->gridx = xgrids[rankx - 1] + 1, iter->gridy = ygrids[ranky - 1] + 1;
				}
				//}

				break;
			case 7: //lr
				//for (INT i = 1; i < set.nGridy - 1; i++)
				//{
				for (INT j = 0; j < set.particlesInGridpointers[IBR(1, set.iGridxm, set.nGridx)].size(); j++)
				{
					passbuf = *set.particlesInGridpointers[IBR(1, set.iGridxm, set.nGridx)][j];
					sendlist.push_back(passbuf);
					iter = sendlist.end(), iter--;
					iter->gridx = 0, iter->gridy = ygrids[ranky - 1] + 1;
				}
				//}
				break;
			default:
				break;
			}

			sendto[r].resize(sendlist.size());
			sendsize[r] = sendlist.size();
			INT i = 0;
			for (iter = sendlist.begin(); iter != sendlist.end(); iter++, i++)
			{
				sendto[r][i] = *iter;
			}
			MPI_Isend(sendsize + r, 1, MPI_INT, ranks[r], 101, commc, sendsreqs + r);
			MPI_Isend(&sendto[r][0], sendlist.size() * prsize, MPI_CHAR, ranks[r], 100, commc, senddreqs + r);
		}
		else
		{
			senddreqs[r] = sendsreqs[r] = MPI_REQUEST_NULL;
		}
	}


	for (int r = 0; r < 8; r++)
	{
		if (ranks[r] >= 0)
		{
			set.outerones[r].clear();
			MPI_Irecv(recvsize + r, 1, MPI_INT, ranks[r], 101, commc, recvsreqs + r);
			MPI_Wait(recvsreqs + r, stats + r);
			recvfrom[r].resize(recvsize[r]);
			MPI_Irecv(&recvfrom[r][0], recvsize[r] * prsize, MPI_CHAR, ranks[r], 100, commc, recvdreqs + r);
			for (int i = 0; i < recvfrom[r].size(); i++) //*****here we depend on the previous correctness of sendreceiveouter
			{
				set.outerones[r][i].rho = recvfrom[r][i].rho;
			}
		}
		else
		{
			set.outerones[r].clear();
			recvdreqs[r] = MPI_REQUEST_NULL;
		}
	}
	MPI_Waitall(8, sendsreqs, stats);
	MPI_Waitall(8, senddreqs, stats);
	MPI_Waitall(8, recvdreqs, stats);

	
}
*/

void SPWater::sendreceiveOuterRho()
{
	vector<vector<SPH_REAL>> recvfrom(8);
	vector<vector<SPH_REAL>> sendto(8);
	int sendsize[8];
	int recvsize[8];

	list<SPH_REAL> sendlist;
	MPI_Status stats[8];
	
	list<SPH_REAL>::iterator iter;

	for (int r = 0; r < 8; r++)
	{
		if (ranks[r] >= 0)
		{
			sendlist.clear();
			switch (r)
			{
			case 0: //upper
				for (size_t i = 1; i <= set.iGridxm; i++)
				{
					for (size_t j = 0; j < set.particlesInGridpointers[IBR(set.iGridym, i, set.nGridx)].size(); j++)
					{
						sendlist.push_back(set.particlesInGridpointers[IBR(set.iGridym, i, set.nGridx)][j]->rho);
						iter = sendlist.end(), iter--;
					}
				}
				break;
			case 1: //lower
				for (size_t i = 1; i <= set.iGridxm; i++)
				{
					for (size_t j = 0; j < set.particlesInGridpointers[IBR(1, i, set.nGridx)].size(); j++)
					{
						sendlist.push_back(set.particlesInGridpointers[IBR(1, i, set.nGridx)][j]->rho);
						iter = sendlist.end(), iter--;
					}
				}
				break;
			case 2: //left
				for (size_t i = 1; i <= set.iGridym; i++)
				{
					for (size_t j = 0; j < set.particlesInGridpointers[IBR(i, 1, set.nGridx)].size(); j++)
					{
						sendlist.push_back(set.particlesInGridpointers[IBR(i, 1, set.nGridx)][j]->rho);
						iter = sendlist.end(), iter--;
					}
				}

				break;
			case 3: //right
				for (size_t i = 1; i <= set.iGridym; i++)
				{
					for (size_t j = 0; j < set.particlesInGridpointers[IBR(i, set.iGridxm, set.nGridx)].size(); j++)
					{
						sendlist.push_back(set.particlesInGridpointers[IBR(i, set.iGridxm, set.nGridx)][j]->rho);
						iter = sendlist.end(), iter--;
					}
				}
				break;
			case 4: //ul
				//for (INT i = 1; i < set.nGridy - 1; i++)
				//{
				for (size_t j = 0; j < set.particlesInGridpointers[IBR(set.iGridym, 1, set.nGridx)].size(); j++)
				{
					sendlist.push_back(set.particlesInGridpointers[IBR(set.iGridym, 1, set.nGridx)][j]->rho);
					iter = sendlist.end(), iter--;
				}
				//}
				break;
			case 5: //ur
				//for (INT i = 1; i < set.nGridy - 1; i++)
				//{
				for (size_t j = 0; j < set.particlesInGridpointers[IBR(set.iGridym, set.iGridxm, set.nGridx)].size(); j++)
				{
					sendlist.push_back(set.particlesInGridpointers[IBR(set.iGridym, set.iGridxm, set.nGridx)][j]->rho);
					iter = sendlist.end(), iter--;
				}
				//}
				break;
			case 6: //ll
				//for (INT i = 1; i < set.nGridy - 1; i++)
				//{
				for (size_t j = 0; j < set.particlesInGridpointers[IBR(1, 1, set.nGridx)].size(); j++)
				{
					sendlist.push_back(set.particlesInGridpointers[IBR(1, 1, set.nGridx)][j]->rho);
					iter = sendlist.end(), iter--;
				}
				//}

				break;
			case 7: //lr
				//for (INT i = 1; i < set.nGridy - 1; i++)
				//{
				for (size_t j = 0; j < set.particlesInGridpointers[IBR(1, set.iGridxm, set.nGridx)].size(); j++)
				{
					sendlist.push_back(set.particlesInGridpointers[IBR(1, set.iGridxm, set.nGridx)][j]->rho);
					iter = sendlist.end(), iter--;
				}
				//}
				break;
			default:
				break;
			}

			sendto[r].resize(sendlist.size());
			sendsize[r] = sendlist.size();
			size_t i = 0;
			for (iter = sendlist.begin(); iter != sendlist.end(); iter++, i++)
			{
				sendto[r][i] = *iter;
			}
			MPI_Isend(sendsize + r, 1, SPH_MPI_INT, ranks[r], 101, commc, sendsreqs + r);
			if (sendlist.size() > 0)
				MPI_Isend(&sendto[r][0], sendlist.size() * sizeof(SPH_REAL), MPI_CHAR, ranks[r], 100, commc, senddreqs + r);
			else
				senddreqs[r] = MPI_REQUEST_NULL;
		}
		else
		{
			senddreqs[r] = sendsreqs[r] = MPI_REQUEST_NULL;
		}
	}


	for (int r = 0; r < 8; r++)
	{
		if (ranks[r] >= 0)
		{
			//set.outerones[r].clear(); this can't happen for this is addition
			MPI_Irecv(recvsize + r, 1, SPH_MPI_INT, ranks[r], 101, commc, recvsreqs + r);
			MPI_Wait(recvsreqs + r, stats + r);
 			recvfrom[r].resize(recvsize[r]);
			if (recvsize[r] > 0)
				MPI_Irecv(&recvfrom[r][0], recvsize[r] * sizeof(SPH_REAL),MPI_CHAR, ranks[r], 100, commc, recvdreqs + r);
			else
				recvdreqs[r] = MPI_REQUEST_NULL;
		}
		else
		{
			recvdreqs[r] = MPI_REQUEST_NULL;
		}
	}	
	//cout << rank << " Start waiting" << endl;
	MPI_Waitall(8, sendsreqs, stats);
	//cout << rank << " Start waiting" << endl;
	MPI_Waitall(8, senddreqs, stats);
	//cout << rank << " Start waiting" << endl;
	MPI_Waitall(8, recvdreqs, stats);
	
	//cout << rank << " Start waiting" << endl;
	
	for (int r = 0; r < 8; r++)
	{
		if (ranks[r] >= 0)
			for (int i = 0; i < recvfrom[r].size(); i++)
			{
				set.outerones[r][i].rho = recvfrom[r][i];
			}
	}

	//MPI_Barrier(commc);
}// send border and receive outer ones' rho, upon those of sendreiveOuters

void SPWater::sendreceivePassers()
{
	unsigned char gridtype;
	unsigned char commtype;
	unsigned char innetype;
	list<ParticleProperties>::iterator iter;
	list<ParticleProperties> sendlist[8];

	vector<ParticleProperties> sendbuf[8];
	vector<ParticleProperties> recvbuf[8];
	int sendsize[8];
	int recvsize[8];

	/*
	MPI_Request sendsreqs[8];
	MPI_Request senddreqs[8];
	MPI_Request recvsreqs[8];
	MPI_Request recvdreqs[8];
	*/
	MPI_Status stats[8];

	int sendr;
	for (iter = set.p.begin(); iter != set.p.end(); )
	{
		gridtype = set.gridType[IBR(iter->gridy, iter->gridx, set.nGridx)];
		innetype = gridtype & SPGTInnerType;
		if ( innetype == SPGT_ExchangeOuter)
		{
			commtype = SPGTCommType & gridtype;
			switch (commtype)
			{
			case SPGT_Upper:
				sendr = 0; break;
			case SPGT_Lower:
				sendr = 1; break;
			case SPGT_Left:
				sendr = 2; break;
			case SPGT_Right:
				sendr = 3; break;
			case SPGT_Upper | SPGT_Left:
				sendr = 4; break;
			case SPGT_Upper | SPGT_Right:
				sendr = 5; break;
			case SPGT_Lower | SPGT_Left:
				sendr = 6; break;
			case SPGT_Lower | SPGT_Right:
				sendr = 7; break;
			default:
				break;
			}
			sendlist[sendr].push_back(*iter);
			iter = set.p.erase(iter);
		}
		else
		{
			iter++;
		}
	}
	int i;
	for (sendr = 0; sendr < 8; sendr++)
	{
		if (ranks[sendr] >= 0)
		{
			sendsize[sendr] = sendlist[sendr].size();
			sendbuf[sendr].resize(sendsize[sendr]);
			i = 0;
			for (iter = sendlist[sendr].begin(); iter != sendlist[sendr].end(); iter++, i++)
				sendbuf[sendr][i] = *iter;
			MPI_Isend(&sendsize[sendr], 1, SPH_MPI_INT, ranks[sendr], 101, commc, sendsreqs + sendr);
			if (sendsize[sendr] > 0)
				MPI_Isend(&sendbuf[sendr][0], sendsize[sendr] * ssize, MPI_CHAR, ranks[sendr], 100, commc, senddreqs + sendr);
			else
				senddreqs[sendr] = MPI_REQUEST_NULL;
		}
		else
		{
			senddreqs[sendr] = sendsreqs[sendr] = MPI_REQUEST_NULL;
		}
	}
	for (int recvr = 0; recvr < 8; recvr++)
	{
		if (ranks[recvr] >= 0)
		{
			MPI_Irecv(&recvsize[recvr], 1, SPH_MPI_INT, ranks[recvr], 101, commc, recvsreqs + recvr);
			MPI_Wait(recvsreqs + recvr, stats + recvr);
			recvbuf[recvr].resize(recvsize[recvr]);
			if (recvsize[recvr] > 0)
				MPI_Irecv(&recvbuf[recvr][0], recvsize[recvr] * ssize, MPI_CHAR, ranks[recvr], 100, commc, recvdreqs + recvr);
			else
				recvdreqs[recvr] = MPI_REQUEST_NULL;
		}
		else
		{
			recvdreqs[recvr] = MPI_REQUEST_NULL;
		}
	}
	//MPI_Barrier(commc);
	MPI_Waitall(8, sendsreqs, stats);
	MPI_Waitall(8, senddreqs, stats);
	MPI_Waitall(8, recvdreqs, stats);

	for (int recvr = 0; recvr < 8; recvr++)
	{
		if (ranks[recvr] >= 0)
		{
			for (i = 0; i < recvbuf[recvr].size(); i++)
			{
				recvbuf[recvr][i].gridx = set.getGridx(recvbuf[recvr][i].x);
				recvbuf[recvr][i].gridy = set.getGridy(recvbuf[recvr][i].y);
				set.p.push_back(recvbuf[recvr][i]);
			}
		}
	}
	//MPI_Barrier(commc);
}

void SPWater::balanceAllocate()
{
	
}

void SPWater::output(string path)
{
	ofstream out(path);
	list<ParticleProperties>::iterator iter;
	for (iter = set.p.begin(); iter != set.p.end(); iter++)
	{
		out << iter->x << '\t' << iter->y << '\t' << iter->vx << '\t' << iter->vy << endl;
	}
}

void SPWater::innerInGridPointer()
{
	list<ParticleProperties>::iterator iter;
	size_t gx, gy;
	unsigned char gridtype;
	for (gy = 1; gy <= set.iGridym; gy++)
		for (gx = 1; gx <= set.iGridxm; gx++)
		{
			gridtype = set.gridType[IBR(gy, gx, set.nGridx)];
			if ((gridtype & SPGTInnerType) == SPGT_SolidWall)
				continue;
			set.particlesInGridpointers[IBR(gy,gx,set.nGridx)].clear();
		}
	for (iter = set.p.begin(); iter != set.p.end(); iter++)
	{
		set.particlesInGridpointers[IBR(iter->gridy, iter->gridx, set.nGridx)].push_back(&(*iter));
	}
}


