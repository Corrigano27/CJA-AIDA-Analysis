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
std::vector <TH2D*> EdT[4];
std::vector <TH2D*> ImplantEdT;
std::vector <TH2D*> NoiseImplantEDecayE;
std::vector <TH2D*> ImplantEDecayE;


//template histograms

TH1D * implantBetaHis;
TH2D * implantBetaHis2D;

//pid
TH2D * PID;

TH2D * edT_All;

TH1D * dT_All;

TH1D * NoiseEnergy;

TH1D * NoiseDT;

TH1D * Energy;

TH1D * FirstPeakEnergy;

TH1D * SecondPeakEnergy;

TH1D * dT_thres_700;

TH1D * dT_thres_750;

TH2D * edT_All_beforeVeto;

TH2D * TePIDfast;

TH2D * IPIDfast;

TH2D * TeEdT;

TH2D * IEdT;
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

	PID = new TH2D("PID","",1e3,1.95,2.35,1e3,39,56);

	TePIDfast = new TH2D("TePIDfast","",1e3,1.95,2.35,1e3,39,56);

	IPIDfast = new TH2D("IPIDfast","",1e3,1.95,2.35,1e3,39,56);

	edT_All = new TH2D("edT_All","",400,-2000,2000,100,0,10000);

	dT_All = new TH1D("dT_All","",400,-2000,-2000);

	edT_All_beforeVeto = new TH2D("edT_All_beforeVeto","",400,-2000,2000,100,0,10000);

	TeEdT = new TH2D("TeEdT","",500,-10,10,100,0,10000);

	IEdT = new TH2D("IEdT","",500,-10,10,100,0,10000);

	NoiseEnergy = new TH1D("NoiseEnergy","",500,0,2000);

	NoiseDT = new TH1D("NoiseDT","",200,0,200);

	Energy = new TH1D("Energy","",1250,0,5000);

	FirstPeakEnergy = new TH1D("FirstPeakEnergy","",500,0,2000);

	SecondPeakEnergy = new TH1D("SecondPeakEnergy","",500,0,2000);

	dT_thres_700 = new TH1D("dT_thres_700","",500,-5000, 5000);

    dT_thres_750 = new TH1D("dT_thres_750","",500,-5000, 5000);

	for (int z = 0; z<6; z++){
		hisName = "ImplantEdT_DSSD" + std::to_string(z);
		implantBetaHis2D = new TH2D(hisName.c_str(), "", 400, -2000, 2000, 100, 0, 7000);
		ImplantEdT.push_back(implantBetaHis2D);

		hisName = "NoiseImplantEDecayE_DSSD" + std::to_string(z);
		implantBetaHis2D = new TH2D(hisName.c_str(), "", 100, 0, 7, 100, 0, 7000);
		NoiseImplantEDecayE.push_back(implantBetaHis2D);

		hisName = "ImplantEDecayE_DSSD" + std::to_string(z);
		implantBetaHis2D = new TH2D(hisName.c_str(), "", 100, 0, 7, 100, 0, 7000);
		ImplantEDecayE.push_back(implantBetaHis2D);

	}

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
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 200, -2000, 2000, 200, 0, 10000);
			EdT[i].push_back(implantBetaHis2D);

		}
	}

}

#endif