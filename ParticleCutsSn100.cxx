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
std::vector <TH2D*> edT[4];
std::vector <TH2D*> EDiff[4];
std::vector <TH2D*> edTLong[4];
std::vector <TH2D*> EDiffLong[4];
std::vector <TH2D*> edTMid[4];
std::vector <TH2D*> edTMid2[4];
std::vector <TH1D*> dTMidUnder520[4];
std::vector <TH1D*> dTMidOver520[4];
std::vector <TH1D*> dTMid2Under520[4];
std::vector <TH1D*> dTMid2Over520[4];

//template histograms

TH1D * implantBetaHis;
TH2D * implantBetaHis2D;

//pid
TH2D * PID;

TH2D * edT_All;

TH2D * edT_All_beforeVeto;

//range of isotopes for each element
void SetParticles(){

	elements[0] = "Ag";
	isotopeStart[0] = 93;
	isotopeEnd[0] = 98;

	elements[1] = "Cd";
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

	//Cd isotopes 94-99
	isotopeDSSDStart[1].push_back(2);
	isotopeDSSDEnd[1].push_back(2);

	isotopeDSSDStart[1].push_back(2);
	isotopeDSSDEnd[1].push_back(2);

	isotopeDSSDStart[1].push_back(1);
	isotopeDSSDEnd[1].push_back(2);

	isotopeDSSDStart[1].push_back(1);
	isotopeDSSDEnd[1].push_back(2);

	isotopeDSSDStart[1].push_back(1);
	isotopeDSSDEnd[1].push_back(2);

	isotopeDSSDStart[1].push_back(1);
	isotopeDSSDEnd[1].push_back(2);

	//In isotopes 94-99
	isotopeDSSDStart[2].push_back(2);
	isotopeDSSDEnd[2].push_back(2);

	isotopeDSSDStart[2].push_back(1);
	isotopeDSSDEnd[2].push_back(2);

	isotopeDSSDStart[2].push_back(1);
	isotopeDSSDEnd[2].push_back(2);

	isotopeDSSDStart[2].push_back(1);
	isotopeDSSDEnd[2].push_back(2);

	isotopeDSSDStart[2].push_back(1);
	isotopeDSSDEnd[2].push_back(2);

	//Sn isotopes 99-101
	isotopeDSSDStart[3].push_back(1);
	isotopeDSSDEnd[3].push_back(2);

	isotopeDSSDStart[3].push_back(1);
	isotopeDSSDEnd[3].push_back(2);

	isotopeDSSDStart[3].push_back(1);
	isotopeDSSDEnd[3].push_back(2);

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
			std::cout << particleCut << std::endl;
			//test->Draw();

		}
	}
}

void DefineHistograms(){

	std::string hisName;

	PID = new TH2D("PID","",1e3,1.95,2.35,1e3,36,53);

	edT_All = new TH2D("edT_All","",100,0,2000,100,0,10000);

	edT_All_beforeVeto = new TH2D("edT_All_beforeVeto","",100,0,2000,100,0,10000);

	for (int i = 0; i<numElements; i++){
		for (int j = 0; j <= isotopeEnd[i]-isotopeStart[i]; j++){

			hisName = elements[i] + std::to_string(isotopeStart[i]+j) + "ImplantBeta";
			implantBetaHis = new TH1D(hisName.c_str(), "", 2e2, -10, 10);
			implantBeta[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i]+j) + "ImplantBeta1p";
			implantBetaHis = new TH1D(hisName.c_str(), "", 2e2, -10, 10);
			implantBeta1p[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i]+j) + "Implant1p";
			implantBetaHis = new TH1D(hisName.c_str(), "", 2e2, -10, 10);
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

			hisName = elements[i] + std::to_string(isotopeStart[i]+j) + "EdT";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 200, 0, 100, 200, 0, 3000);
			edT[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i]+j) + "EDiff";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 200,0,200,500,-2000,2000);
			EDiff[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i]+j) + "EdTLong";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 200, 0, 10000, 200, 0, 3000);
			edTLong[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i]+j) + "EDiffLong";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 100,0,10,500,-2000,2000);
			EDiffLong[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i]+j) + "EdTMid";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 200, 0, 1000, 200, 0, 10000);
			edTMid[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i]+j) + "EdTMid2";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 200, 1000, 2000, 200, 0, 10000);
			edTMid2[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i]+j) + "dTMidUnder520";
			implantBetaHis = new TH1D(hisName.c_str(),"",200,0,1000);
			dTMidUnder520[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i]+j) + "dTMid2Under520";
			implantBetaHis = new TH1D(hisName.c_str(),"",200,1000,2000);
			dTMid2Under520[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i]+j) + "dTMidOver520";
			implantBetaHis = new TH1D(hisName.c_str(),"",200,0,1000);
			dTMidOver520[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i]+j) + "dTMid2Over520";
			implantBetaHis = new TH1D(hisName.c_str(),"",200,1000,2000);
			dTMid2Over520[i].push_back(implantBetaHis);



		}
	}

}

#endif