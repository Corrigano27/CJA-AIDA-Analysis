#ifndef PARTICLECUTSSN100_H
#define PARTICLECUTSSN100_H

#include "TFile.h"
#include "TTree.h"
#include <vector>
#include "TBranch.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCut.h"
#include "TCutG.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

int numElements = 4;

//string declarations;
std::string particleCut; //PID cuts from file

std::string elements[4]; //isotope elements

//vectors of isotope mass numbers
int isotopeStart[4];
int isotopeEnd[4];

std::vector <int> isotopeDSSDStart[4];
std::vector <int> isotopeDSSDEnd[4];

//graphical cuts for isotopes
TCutG * particleCuts[4][8]; //indices denote elements and isotopes

//histogram vectors
std::vector <TH1D*> implantBeta[4];
std::vector <TH1D*> implantBeta1p[4];
std::vector <TH1D*> implant1p[4];
std::vector <TH1D*> decayEnergy[4];
std::vector <TH1D*> delayed1pEnergy[4];
std::vector <TH1D*> implantZ[4];

//template histogram
TH1D * implantBetaHis;

//pid
TH2D * PID;

//range of isotopes for each element
void SetParticles(){

	elements[0] = "Ag";
	isotopeStart[0] = 93;
	isotopeEnd[0] = 98;

	elements[1] = "Rh";
	isotopeStart[1] = 94;
	isotopeEnd[1] = 99;

	elements[2] = "In";
	isotopeStart[2] = 96;
	isotopeEnd[2] = 100;

	elements[3] = "Sn";
	isotopeStart[3] = 99;
	isotopeEnd[3] = 101;

	
}

//Set DSSD layers of expected implantation for each isotope
void SetImplantDSSD(){

	//Ag isotopes 93-98
	isotopeDSSDStart[0].push_back(2);
	isotopeDSSDEnd[0].push_back(2);

	isotopeDSSDStart[0].push_back(2);
	isotopeDSSDEnd[0].push_back(2);

	isotopeDSSDStart[0].push_back(1);
	isotopeDSSDEnd[0].push_back(2);

	isotopeDSSDStart[0].push_back(1);
	isotopeDSSDEnd[0].push_back(2);

	isotopeDSSDStart[0].push_back(1);
	isotopeDSSDEnd[0].push_back(2);

	isotopeDSSDStart[0].push_back(1);
	isotopeDSSDEnd[0].push_back(2);

	//Rh isotopes 94-99
	isotopeDSSDStart[0].push_back(2);
	isotopeDSSDEnd[0].push_back(2);

	isotopeDSSDStart[0].push_back(2);
	isotopeDSSDEnd[0].push_back(2);

	isotopeDSSDStart[0].push_back(1);
	isotopeDSSDEnd[0].push_back(2);

	isotopeDSSDStart[0].push_back(1);
	isotopeDSSDEnd[0].push_back(2);

	isotopeDSSDStart[0].push_back(1);
	isotopeDSSDEnd[0].push_back(2);

	isotopeDSSDStart[0].push_back(1);
	isotopeDSSDEnd[0].push_back(2);

	//In isotopes 94-99
	isotopeDSSDStart[0].push_back(2);
	isotopeDSSDEnd[0].push_back(2);

	isotopeDSSDStart[0].push_back(1);
	isotopeDSSDEnd[0].push_back(2);

	isotopeDSSDStart[0].push_back(1);
	isotopeDSSDEnd[0].push_back(2);

	isotopeDSSDStart[0].push_back(1);
	isotopeDSSDEnd[0].push_back(2);

	isotopeDSSDStart[0].push_back(1);
	isotopeDSSDEnd[0].push_back(2);

	//Rh isotopes 94-99
	isotopeDSSDStart[0].push_back(1);
	isotopeDSSDEnd[0].push_back(2);

	isotopeDSSDStart[0].push_back(1);
	isotopeDSSDEnd[0].push_back(2);

	isotopeDSSDStart[0].push_back(1);
	isotopeDSSDEnd[0].push_back(2);

}

//reads i cut objects from file
void ReadParticleCuts(std::string cutFile){

	TFile * particleCutFile = TFile::Open(cutFile.c_str(),"READ");

	for (int i = 0; i<numElements; i++){
		for (int j = 0; j <= isotopeEnd[i]-isotopeStart[i]; j++){

			particleCut = "CUTG" + std::to_string(isotopeStart[i]+j) + elements[i];
			particleCuts[i][j] = (TCutG*) particleCutFile->Get(particleCut.c_str());
			//TCutG* test = (TCutG *)particleCutFile->Get(particleCut.c_str()); //seg faults here, to do with particleCut.c_str() string!
			//std::cout << "error drawing" << endl;
			std::cout << particleCut << endl;
			//test->Draw();

		}
	}
}

void DefineHistograms(){

	std::string hisName;

	PID = new TH2D("PID","",1e3,1.95,2.35,1e3,36,53);

	for (int i = 0; i<numElements; i++){
		for (int j = 0; j <= isotopeEnd[i]-isotopeStart[i]; j++){

			hisName = elements[i] + std::to_string(isotopeStart[i]+j) + "ImplantBeta";
			implantBetaHis = new TH1D(hisName.c_str(), "", 2e4, -10, 10);
			implantBeta[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i]+j) + "ImplantBeta1p";
			implantBetaHis = new TH1D(hisName.c_str(), "", 2e4, -10, 10);
			implantBeta1p[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i]+j) + "Implant1p";
			implantBetaHis = new TH1D(hisName.c_str(), "", 2e4, -10, 10);
			implant1p[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i]+j) + "DecayEnergy";
			implantBetaHis = new TH1D(hisName.c_str(),"",5e2,0,5000);
			decayEnergy[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i]+j) + "Delayed1pEnergy";
			implantBetaHis = new TH1D(hisName.c_str(), "", 150, 1500, 5000);
			delayed1pEnergy[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i]+j) + "ImplantZ";
			implantBetaHis = new TH1D(hisName.c_str(), "", 6, 0, 6);
			implantZ[i].push_back(implantBetaHis);



		}
	}

}

#endif