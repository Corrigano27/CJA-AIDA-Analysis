#ifndef 100SNPARTICLECUTS_H
#define 100SNPARTICLECUTS_H

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

//graphical cuts for isotopes
TCutG * particleCuts[4][8]; //indices denote elements and isotopes

//vectors of isotope mass numbers
std::vector <int> isotopeStart[4];
std::vector <int> isotopeEnd[4];

//histogram vectors
std::vector <TH1D*> implantBeta[4];
std::vector <TH1D*> implantBeta1p[4];
std::vector <TH1D*> implant1p[4];
std::vector <TH1D*> decayEnergy[4];
std::vector <TH1D*> delayed1pEnergy[4];
std::vector <TH1D*> implantZ[4];

//template histogram
TH1D * implantBetaHis;

//range of isotopes for each element
void SetParticles(){
	elements[0] = "Sn";
	isotopeStart[0] = 99;
	isotopeEnd[0] = 101;

	elements[1] = "In";
	isotopeStart[1] = 96;
	isotopeEnd[1] = 100;

	elements[2] = "Rh";
	isotopeStart[2] = 94;
	isotopeEnd[2] = 99;

	elements[3] = "Ag";
	isotopeStart[3] = 93;
	isotopeEnd[3] = 98;

}

//Set DSSD layers of expected implantation
void SetImplantDSSD(){

}

//reads i cut objects from file
void ReadParticleCuts(std::string cutFile){

	TFile * particleCutFile = TFile::Open(cutFile.c_str(),"READ");

	for (int i = 0; i<numElements; i++){
		for (int j = 0; j <= isotopeEnd[i]-isotopeStart[i]; j++){

			particleCut = "CUTG" + elements[i] + std::to_string(isotopeStart[i]+j);
			particleCuts[i][j] = (TCutG*) particleFile->Get(particleCut.c_str());

		}
	}
}

void DefineHistograms(){

	std::string hisName;

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