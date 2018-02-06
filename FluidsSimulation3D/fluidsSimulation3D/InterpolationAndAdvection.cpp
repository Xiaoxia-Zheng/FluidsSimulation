#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include "slVector.H"
#include "slArray.H"
#include "json.h"
#ifdef __APPLE__
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif
#include <cassert>

#ifndef M_PI
#define	M_PI		3.14159265358979323846	/* pi */
#endif

struct SimulationParameters{
	double dt, total_time, density;
};

class StaggeredGrid{
public:
	SlArray3D<double> u, v, w;
	SlVector3 lc, uc;
	SlArray3D<unsigned char> cellLatbels;
	int nx, ny, nz;
	double h, halfh;
	
	unsigned int nFlipParticles;
	double flipRatio;
	std::vector<SlVector3>flipVel, flipPos;
	
};




bool readParticles(const char *fname, unsigned int &nFlipParticles, std::vector<SlVector3> &flipPos, std::vector<SlVector3> &flipVel) {
	std::ifstream in(fname, std::ios::in);
//	std::ifstream in("/Users/zhengli/Desktop/CMSC691/project_5/proj5/particles.in", std::ios::in);
	
	in>>nFlipParticles;
	flipPos.resize(nFlipParticles);
	flipVel.resize(nFlipParticles);
	for (unsigned int i=0; i<nFlipParticles; i++) {
		in>>flipPos[i][0]>>flipPos[i][1]>>flipPos[i][2]>>flipVel[i][0]>>flipVel[i][1]>>flipVel[i][2];
	}
	in.close();
	std::cout<<"inputfile " << fname <<" read"<<std::endl;
	return true;
}


bool writeParticles(const char *fname, unsigned int &nFlipParticles, std::vector<SlVector3> &flipPos, std::vector<SlVector3> &flipVel) {
	std::ofstream out(fname, std::ios::out);
	out<<nFlipParticles<<std::endl;
	for (unsigned int i=0; i<nFlipParticles; i++) {
		out<<flipPos[i][0]<<" "<<flipPos[i][1]<<" "<<flipPos[i][2]<<" "<<flipVel[i][0]<<" "<<flipVel[i][1]<<" "<<flipVel[i][2]<<std::endl;
	}
	out.close();
	std::cout<<"outputfile " << fname <<" written"<<std::endl;
	return true;
}


bool readInputFile(const char *fname,
				   SimulationParameters &params,
				   StaggeredGrid &grid) {
	std::ifstream in(fname, std::ios::in);
	Json::Value root;
	Json::Reader jReader;
	if(!jReader.parse(in, root)){
		std::cout << "couldn't read input file: " << fname << '\n'
		<< jReader.getFormattedErrorMessages() << std::endl;
		exit(1);
	}
	params.dt = root.get("dt", 1.0/300.0).asDouble();
	params.total_time = root.get("total_time", 1.0).asDouble();
	params.density = root.get("density", 1.0/300.0).asDouble();
	grid.flipRatio = root.get("flipRatio", 0.95).asDouble();
	grid.nx = root["res"][0].asInt();
	grid.ny = root["res"][1].asInt();
	grid.nz = root["res"][2].asInt();
	grid.h = root["h"].asDouble();
	SlVector3 lc, uc;
	lc.set(root["lc"][0].asDouble(), root["lc"][1].asDouble(), root["lc"][2].asDouble());
	uc.set(root["uc"][0].asDouble(), root["uc"][1].asDouble(), root["uc"][2].asDouble());
//	uc = lc + SlVector3(grid.h * grid.nx, grid.h * grid.ny, grid.h * grid.nz);
	uc = lc + SlVector3(grid.nx, grid.ny, grid.nz);
	grid.lc = lc;
	grid.uc = uc;
//	grid.allocate(nx, ny, nz, h, lc, uc);
//	printf("test :%s\n", root["particles"].asString().c_str());
	readParticles(root["particles"].asString().c_str(), grid.nFlipParticles, grid.flipPos, grid.flipVel);
	return true;
}


//update particle's next position
bool advection(std::vector<SlVector3> &flipPos, std::vector<SlVector3> &flipVel, double dt){
	for (int particleIndex = 0; particleIndex < flipPos.size(); particleIndex++) {
		flipPos[particleIndex][0] += flipVel[particleIndex][0] * dt;
		flipPos[particleIndex][1] += flipVel[particleIndex][1] * dt;
		flipPos[particleIndex][2] += flipVel[particleIndex][2] * dt;
		
		//boundary
		flipPos[particleIndex][0] = std::min(0.125, flipPos[particleIndex][0]);
		flipPos[particleIndex][0] = std::max(-0.125, flipPos[particleIndex][0]);

		flipPos[particleIndex][1] = std::min(0.125, flipPos[particleIndex][0]);
		flipPos[particleIndex][1] = std::max(-0.125, flipPos[particleIndex][0]);

		flipPos[particleIndex][2] = std::min(0.125, flipPos[particleIndex][0]);
		flipPos[particleIndex][2] = std::max(-0.125, flipPos[particleIndex][0]);

		
//		printf("test: %lf\n", flipPos[0][0]);
	}
	return true;
}


//calculate weights base on the position and where did it locate on the grid.
bool particleToGrid(std::vector< std::vector< std::vector<SlVector3> > > &gridVelocity, std::vector<SlVector3> &flipVel ,std::vector<SlVector3> &flipPos, double h, SlVector3 lc){
	
	//     P6     P7
	//  P2     P3
	//	   P4     P5
	//  P0     P1
	
	std::vector< std::vector<SlVector3> > weights;
	std::vector<SlVector3> temp_weight(8, SlVector3(0, 0, 0));
	
//	std::vector< std::vector< std::vector<SlVector3> > > gridVelocity(25, std::vector< std::vector<SlVector3> >(25, std::vector<SlVector3>(25, SlVector3(0, 0, 0))));
	
	for (int i=0; i<flipPos.size(); i++) {
		//8 vertex, 8 weights
		//w0x = (P1x - Px) / (P1x - P0x)
		//w0y = (P2y - Py) / (P2y - P0y)
		//w0z = (P4z - Pz) / (P4z - P0z)
		//grid Index:
		int gridIndex_x = (flipPos[i][0] - lc[0]) / h;
		int gridIndex_y = (flipPos[i][1] - lc[1]) / h;
		int gridIndex_z = (flipPos[i][2] - lc[2]) / h;
		
		//calculate weights:
		temp_weight[0][0] = fabs((lc[0] + gridIndex_x * h) - flipPos[i][0]) / h;
		temp_weight[0][1] = fabs((lc[1] + gridIndex_y * h) - flipPos[i][1]) / h;
		temp_weight[0][2] = fabs((lc[2] + gridIndex_z * h) - flipPos[i][2]) / h;
		
		temp_weight[1][0] = 1 - temp_weight[0][0];
		temp_weight[1][1] = temp_weight[0][1];
		temp_weight[1][2] = temp_weight[0][2];
		
		temp_weight[2][0] = temp_weight[0][0];
		temp_weight[2][1] = 1 - temp_weight[0][1];
		temp_weight[2][2] = temp_weight[0][2];
		
		temp_weight[3][0] = 1 - temp_weight[0][0];
		temp_weight[3][1] = 1 - temp_weight[0][1];
		temp_weight[3][2] = temp_weight[0][2];
		
		temp_weight[4][0] = temp_weight[0][0];
		temp_weight[4][1] = temp_weight[0][1];
		temp_weight[4][2] = 1 - temp_weight[0][2];
		
		temp_weight[5][0] = 1 - temp_weight[0][0];
		temp_weight[5][1] = temp_weight[0][1];
		temp_weight[5][2] = 1 - temp_weight[0][2];
		
		temp_weight[6][0] = temp_weight[0][0];
		temp_weight[6][1] = 1 - temp_weight[0][1];
		temp_weight[6][2] = 1 - temp_weight[0][2];
		
		temp_weight[7][0] = 1 - temp_weight[0][0];
		temp_weight[7][1] = 1 - temp_weight[0][1];
		temp_weight[7][2] = 1 - temp_weight[0][2];
		
		weights.push_back(temp_weight);
		
		
		//update grid velocity
		//P0
		gridVelocity[gridIndex_x - 1][gridIndex_y - 1][gridIndex_z][0] += flipVel[i][0] * weights[i][0][0];
		gridVelocity[gridIndex_x - 1][gridIndex_y - 1][gridIndex_z][1] += flipVel[i][1] * weights[i][0][1];
		gridVelocity[gridIndex_x - 1][gridIndex_y - 1][gridIndex_z][2] += flipVel[i][2] * weights[i][0][2];
		
		//P1
		gridVelocity[gridIndex_x][gridIndex_y - 1][gridIndex_z][0] += flipVel[i][0] * weights[i][1][0];
		gridVelocity[gridIndex_x][gridIndex_y - 1][gridIndex_z][1] += flipVel[i][1] * weights[i][1][1];
		gridVelocity[gridIndex_x][gridIndex_y - 1][gridIndex_z][2] += flipVel[i][2] * weights[i][1][2];
		
		//P2
		gridVelocity[gridIndex_x - 1][gridIndex_y][gridIndex_z][0] += flipVel[i][0] * weights[i][2][0];
		gridVelocity[gridIndex_x - 1][gridIndex_y][gridIndex_z][1] += flipVel[i][1] * weights[i][2][1];
		gridVelocity[gridIndex_x - 1][gridIndex_y][gridIndex_z][2] += flipVel[i][2] * weights[i][2][2];
		
		//P3
		gridVelocity[gridIndex_x][gridIndex_y][gridIndex_z][0] += flipVel[i][0] * weights[i][3][0];
		gridVelocity[gridIndex_x][gridIndex_y][gridIndex_z][1] += flipVel[i][1] * weights[i][3][1];
		gridVelocity[gridIndex_x][gridIndex_y][gridIndex_z][2] += flipVel[i][2] * weights[i][3][2];
		
		//P4
		gridVelocity[gridIndex_x - 1][gridIndex_y - 1][gridIndex_z - 1][0] += flipVel[i][0] * weights[i][4][0];
		gridVelocity[gridIndex_x - 1][gridIndex_y - 1][gridIndex_z - 1][1] += flipVel[i][1] * weights[i][4][1];
		gridVelocity[gridIndex_x - 1][gridIndex_y - 1][gridIndex_z - 1][2] += flipVel[i][2] * weights[i][4][2];

		//P5
		gridVelocity[gridIndex_x][gridIndex_y - 1][gridIndex_z - 1][0] += flipVel[i][0] * weights[i][5][0];
		gridVelocity[gridIndex_x][gridIndex_y - 1][gridIndex_z - 1][1] += flipVel[i][1] * weights[i][5][1];
		gridVelocity[gridIndex_x][gridIndex_y - 1][gridIndex_z - 1][2] += flipVel[i][2] * weights[i][5][2];
		
		//P6
		gridVelocity[gridIndex_x - 1][gridIndex_y][gridIndex_z - 1][0] += flipVel[i][0] * weights[i][6][0];
		gridVelocity[gridIndex_x - 1][gridIndex_y][gridIndex_z - 1][1] += flipVel[i][1] * weights[i][6][1];
		gridVelocity[gridIndex_x - 1][gridIndex_y][gridIndex_z - 1][2] += flipVel[i][2] * weights[i][6][2];
		
		//P7
		gridVelocity[gridIndex_x][gridIndex_y][gridIndex_z - 1][0] += flipVel[i][0] * weights[i][7][0];
		gridVelocity[gridIndex_x][gridIndex_y][gridIndex_z - 1][1] += flipVel[i][1] * weights[i][7][1];
		gridVelocity[gridIndex_x][gridIndex_y][gridIndex_z - 1][2] += flipVel[i][2] * weights[i][7][2];
		
	}
	return true;
}

bool gridToparticle(std::vector<SlVector3> &flipPos, std::vector<SlVector3> &flipVel, std::vector< std::vector< std::vector<SlVector3> > > &gridVelocity, SlVector3 lc, double h){

	std::vector< std::vector<SlVector3> > weights;
	std::vector<SlVector3> temp_weight(8, SlVector3(0, 0, 0));
	
	//	std::vector< std::vector< std::vector<SlVector3> > > gridVelocity(25, std::vector< std::vector<SlVector3> >(25, std::vector<SlVector3>(25, SlVector3(0, 0, 0))));
	
	for (int i=0; i<flipPos.size(); i++) {
		//8 vertex, 8 weights
		//w0x = (P1x - Px) / (P1x - P0x)
		//w0y = (P2y - Py) / (P2y - P0y)
		//w0z = (P4z - Pz) / (P4z - P0z)
		//grid Index:
		int gridIndex_x = (flipPos[i][0] - lc[0]) / h;
		int gridIndex_y = (flipPos[i][1] - lc[1]) / h;
		int gridIndex_z = (flipPos[i][2] - lc[2]) / h;
		
		//calculate weights:
		temp_weight[0][0] = fabs((lc[0] + gridIndex_x * h) - flipPos[i][0]) / h;
		temp_weight[0][1] = fabs((lc[1] + gridIndex_y * h) - flipPos[i][1]) / h;
		temp_weight[0][2] = fabs((lc[2] + gridIndex_z * h) - flipPos[i][2]) / h;
		
		temp_weight[1][0] = 1 - temp_weight[0][0];
		temp_weight[1][1] = temp_weight[0][1];
		temp_weight[1][2] = temp_weight[0][2];
		
		temp_weight[2][0] = temp_weight[0][0];
		temp_weight[2][1] = 1 - temp_weight[0][1];
		temp_weight[2][2] = temp_weight[0][2];
		
		temp_weight[3][0] = 1 - temp_weight[0][0];
		temp_weight[3][1] = 1 - temp_weight[0][1];
		temp_weight[3][2] = temp_weight[0][2];
		
		temp_weight[4][0] = temp_weight[0][0];
		temp_weight[4][1] = temp_weight[0][1];
		temp_weight[4][2] = 1 - temp_weight[0][2];
		
		temp_weight[5][0] = 1 - temp_weight[0][0];
		temp_weight[5][1] = temp_weight[0][1];
		temp_weight[5][2] = 1 - temp_weight[0][2];
		
		temp_weight[6][0] = temp_weight[0][0];
		temp_weight[6][1] = 1 - temp_weight[0][1];
		temp_weight[6][2] = 1 - temp_weight[0][2];
		
		temp_weight[7][0] = 1 - temp_weight[0][0];
		temp_weight[7][1] = 1 - temp_weight[0][1];
		temp_weight[7][2] = 1 - temp_weight[0][2];
		
		weights.push_back(temp_weight);
		
		
		//update particle velocity
		//X
		flipVel[i][0] =
		gridVelocity[gridIndex_x - 1][gridIndex_y - 1][gridIndex_z][0] * weights[i][0][0] +
		gridVelocity[gridIndex_x][gridIndex_y - 1][gridIndex_z][0] * weights[i][1][0] +
		gridVelocity[gridIndex_x - 1][gridIndex_y][gridIndex_z][0] * weights[i][2][0] +
		gridVelocity[gridIndex_x][gridIndex_y][gridIndex_z][0] * weights[i][3][0] +
		gridVelocity[gridIndex_x - 1][gridIndex_y - 1][gridIndex_z - 1][0] * weights[i][4][0] +
		gridVelocity[gridIndex_x][gridIndex_y - 1][gridIndex_z - 1][0] * weights[i][5][0] +
		gridVelocity[gridIndex_x - 1][gridIndex_y][gridIndex_z - 1][0] * weights[i][6][0] +
		gridVelocity[gridIndex_x][gridIndex_y][gridIndex_z - 1][0] * weights[i][7][0];

		//Y
		flipVel[i][1] =
		gridVelocity[gridIndex_x - 1][gridIndex_y - 1][gridIndex_z][1] * weights[i][0][1] +
		gridVelocity[gridIndex_x][gridIndex_y - 1][gridIndex_z][1] * weights[i][1][1] +
		gridVelocity[gridIndex_x - 1][gridIndex_y][gridIndex_z][1] * weights[i][2][1] +
		gridVelocity[gridIndex_x][gridIndex_y][gridIndex_z][1] * weights[i][3][1]+
		gridVelocity[gridIndex_x - 1][gridIndex_y - 1][gridIndex_z - 1][1] * weights[i][4][1] +
		gridVelocity[gridIndex_x][gridIndex_y - 1][gridIndex_z - 1][1] * weights[i][5][1] +
		gridVelocity[gridIndex_x - 1][gridIndex_y][gridIndex_z - 1][1] * weights[i][6][1] +
		gridVelocity[gridIndex_x][gridIndex_y][gridIndex_z - 1][1] * weights[i][7][1];
		
		//Z
		flipVel[i][2] =
		gridVelocity[gridIndex_x - 1][gridIndex_y - 1][gridIndex_z][2] * weights[i][0][2] +
		gridVelocity[gridIndex_x][gridIndex_y - 1][gridIndex_z][2] * weights[i][1][2] +
		gridVelocity[gridIndex_x - 1][gridIndex_y][gridIndex_z][2] * weights[i][2][2] +
		gridVelocity[gridIndex_x][gridIndex_y][gridIndex_z][2] * weights[i][3][2] +
		gridVelocity[gridIndex_x - 1][gridIndex_y - 1][gridIndex_z - 1][2] * weights[i][4][2] +
		gridVelocity[gridIndex_x][gridIndex_y - 1][gridIndex_z - 1][2] * weights[i][5][2] +
		gridVelocity[gridIndex_x - 1][gridIndex_y][gridIndex_z - 1][2] * weights[i][6][2] +
		gridVelocity[gridIndex_x][gridIndex_y][gridIndex_z - 1][2] * weights[i][7][2];
		
	}
	return true;
}



int main(int argc, char *argv[]){
	char fname[100];
	SimulationParameters params;
	StaggeredGrid grid;
	double time = 0;
	int frame = 0;
	double frameTime = -1.0;
	
	std::vector< std::vector< std::vector<SlVector3> > > gridVelocity(25, std::vector< std::vector<SlVector3> >(25, std::vector<SlVector3>(25, SlVector3(0, 0, 0))));
	
	readInputFile(argv[1], params, grid);
	
	while(time < params.total_time){
		//output files every 30 frames.
		if (frame % 30 == 0) {
			sprintf(fname, argv[2], frame/30);
			writeParticles(fname, grid.nFlipParticles, grid.flipPos, grid.flipVel);
		}

		advection(grid.flipPos, grid.flipVel, params.dt);
		
		particleToGrid(gridVelocity, grid.flipVel, grid.flipPos, grid.h, grid.lc);

		gridToparticle(grid.flipPos, grid.flipVel, gridVelocity, grid.lc, grid.h);
		
		time += params.dt;
		frame++;
	}
}
