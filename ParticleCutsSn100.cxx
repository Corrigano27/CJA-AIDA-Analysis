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

const int numElements = 7;

//string declarations;
std::string particleCut; //PID cuts from file

std::string elements[numElements]; //isotope elements


//vectors of isotope mass numbers
int isotopeStart[numElements];
int isotopeEnd[numElements];

std::vector<int> isotopeDSSDStart[numElements];
std::vector<int> isotopeDSSDEnd[numElements];

//graphical cuts for isotopes
TCutG *particleCuts[numElements][8]; //indices denote elements and isotopes

//directories
std::vector<TDirectory *> ElementDir;
std::vector<TDirectory *> IsotopeDir;
std::vector<TDirectory *> DSSD_Dir;

TDirectory *ElDir;
TDirectory *IsoDir;
TDirectory *ZDir;


//histogram vectors
//decays
std::vector<TH1D *> implantBeta[numElements][6];
std::vector<TH1D *> implantBeta1p[numElements][6];
std::vector<TH1D *> implant1p[numElements][6];
std::vector<TH1D *> decayEnergy[numElements][6];
std::vector<TH1D *> delayed1pEnergy[numElements][6];
std::vector<TH1D *> delayed1pEnergyX[numElements][6];

//std::vector<TH1D *> delayed1pEnergyX_0_1_ms[numElements][6];
//std::vector<TH1D *> delayed1pEnergyX_0_100_ms[numElements][6];
//std::vector<TH1D *> delayed1pEnergyX_0_1_s[numElements][6];

//std::vector<TH1D *> delayed1pEnergyY[numElements][6];
std::vector<TH1D *> delayed1pEnergyRandom[numElements][6];
std::vector<TH1D *> delayed1pEnergyAll[numElements][6];
std::vector<TH2D *> EdT[numElements][6];
std::vector<TH2D *> EdT_ms[numElements][6];
std::vector<TH2D *> EdT_us[numElements][6];

std::vector<TH2D *> implantVelocityimplantZ[numElements];
std::vector<TH2D *> implantVelocityimplantE[numElements][6];

//std::vector<TH2D *> EnergyXChannel[numElements][6];
//std::vector<TH2D *> EnergyYChannel[numElements][6];
//std::vector<TH2D *> ExEy[numElements][6];
//std::vector<TH1D *> ExEyDiff[numElements][6];

std::vector<TH1D *> implantBetaAll[numElements];
std::vector<TH1D *> implantBeta1pAll[numElements];
std::vector<TH1D *> implant1pAll[numElements];
std::vector<TH1D *> decayEnergyAll[numElements];
std::vector<TH1D *> delayed1pEnergy_AllDSSD[numElements];

std::vector<TH1D *> delayed1pEnergyRandom_AllDSSD[numElements];
std::vector<TH1D *> delayed1pEnergyAll_AllDSSD[numElements];
std::vector<TH2D *> EdTAll_NoMultiGate_[numElements];
//std::vector<TH2D *> EdTAll_NoMultiGate_corr[numElements][2];
std::vector<TH2D *> EdTAll_ms[numElements];
//std::vector<TH2D *> EdTAll_ms_corr[numElements][2];
std::vector<TH2D *> EdTAll_us[numElements];
//std::vector<TH2D *> EdTAll_us_corr[numElements][2];
std::vector<TH2D *> EdTAll11[numElements];
std::vector<TH2D *> EdTAll12[numElements];
std::vector<TH2D *> EdTAll21[numElements];
std::vector<TH2D *> EdTAll22[numElements];
std::vector<TH2D *> ExEy11[numElements];
std::vector<TH2D *> ExEy12[numElements];
std::vector<TH2D *> ExEy21[numElements];
std::vector<TH2D *> ExEy22[numElements];
//std::vector<TH1D *> implantEAll[numElements];
std::vector<TH2D *> implantVelocityimplantEAll[numElements];

//all DSSD, gated on 511 summed spectra


TH1D *globalEnergy;

TH1D *isotopeSumEnergy;


//end of gammas

//std::vector<TH2D *> EnergyXChannelAll[numElements];
//std::vector<TH2D *> EnergyYChannelAll[numElements];
//std::vector<TH2D *> ExEyAll[numElements];
//std::vector<TH1D *> ExEyDiffAll[numElements];

std::vector<TH2D *> implantVelocityAOQ[numElements][6];
std::vector<TH2D *> implantVelocityAOQ_AllDSSD[numElements];

std::vector<TH2D *> implantEnergyAOQ[numElements][6];
std::vector<TH2D *> implantEnergyAOQ_AllDSSD[numElements];

//implants
std::vector<TH1D *> implantZ[numElements];
//std::vector<TH1D *> implantE[numElements][6];


//template histograms

TH1D *implantBetaHis;
TH2D *implantBetaHis2D;

//pid
TH2D *PID;

//TH2D *PID_noise;

TH2D *PID_implant;

TH2D *Indium97_gammaveto_EdT;

//range of isotopes for each element
void SetParticles()
{

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

	elements[4] = "Pd";
	isotopeStart[4] = 93;
	isotopeEnd[4] = 99;

	elements[5] = "Rh";
	isotopeStart[5] = 92;
	isotopeEnd[5] = 98;

	elements[6] = "Ru";
	isotopeStart[6] = 91;
	isotopeEnd[6] = 97;

}

//Set DSSD layers of expected implantation for each isotope
void SetImplantDSSD()
{

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

	//In isotopes 96-100
	isotopeDSSDStart[2].push_back(2);
	isotopeDSSDEnd[2].push_back(2);

	isotopeDSSDStart[2].push_back(1);
	isotopeDSSDEnd[2].push_back(3);

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

	//Pd isotopes 93-99
	isotopeDSSDStart[4].push_back(2);
	isotopeDSSDEnd[4].push_back(3);

	isotopeDSSDStart[4].push_back(1);
	isotopeDSSDEnd[4].push_back(3);

	isotopeDSSDStart[4].push_back(1);
	isotopeDSSDEnd[4].push_back(2);

	isotopeDSSDStart[4].push_back(1);
	isotopeDSSDEnd[4].push_back(2);

	isotopeDSSDStart[4].push_back(1);
	isotopeDSSDEnd[4].push_back(2);

	isotopeDSSDStart[4].push_back(1);
	isotopeDSSDEnd[4].push_back(2);

	isotopeDSSDStart[4].push_back(1);
	isotopeDSSDEnd[4].push_back(2);

	//Rh isotopes 92-98
	isotopeDSSDStart[5].push_back(2);
	isotopeDSSDEnd[5].push_back(3);

	isotopeDSSDStart[5].push_back(2);
	isotopeDSSDEnd[5].push_back(3);

	isotopeDSSDStart[5].push_back(1);
	isotopeDSSDEnd[5].push_back(2);

	isotopeDSSDStart[5].push_back(1);
	isotopeDSSDEnd[5].push_back(2);

	isotopeDSSDStart[5].push_back(1);
	isotopeDSSDEnd[5].push_back(2);

	isotopeDSSDStart[5].push_back(1);
	isotopeDSSDEnd[5].push_back(2);

	isotopeDSSDStart[5].push_back(1);
	isotopeDSSDEnd[5].push_back(2);

	//Ru isotopes 91-97
	isotopeDSSDStart[6].push_back(2);
	isotopeDSSDEnd[6].push_back(3);

	isotopeDSSDStart[6].push_back(1);
	isotopeDSSDEnd[6].push_back(3);

	isotopeDSSDStart[6].push_back(1);
	isotopeDSSDEnd[6].push_back(2);

	isotopeDSSDStart[6].push_back(1);
	isotopeDSSDEnd[6].push_back(2);

	isotopeDSSDStart[6].push_back(1);
	isotopeDSSDEnd[6].push_back(2);

	isotopeDSSDStart[6].push_back(1);
	isotopeDSSDEnd[6].push_back(2);

	isotopeDSSDStart[6].push_back(1);
	isotopeDSSDEnd[6].push_back(2);
}

double DTAS_SingleCalib(double E){
	double Ecorr;
	Ecorr = 0.986*E - 20.2;
	return Ecorr;
}

//reads i cut objects from file
void ReadParticleCuts(std::string cutFile)
{

	TFile *particleCutFile = TFile::Open(cutFile.c_str(), "READ");

	for (int i = 0; i < numElements; i++)
	{
		for (int j = 0; j <= isotopeEnd[i] - isotopeStart[i]; j++)
		{

			particleCut = "CUTG" + std::to_string(isotopeStart[i] + j) + elements[i];
			particleCuts[i][j] = (TCutG *)particleCutFile->Get(particleCut.c_str());
			//TCutG* test = (TCutG *)particleCutFile->Get(particleCut.c_str()); //seg faults here, to do with particleCut.c_str() string!
			//std::cout << "error drawing" << endl;
			std::cout << particleCut << std::endl;
			//test->Draw();
		}
	}
}

void DefineHistograms()
{

	std::string hisName;

	//std::string Sn101histName;

	//std::string Ag95histName;

	PID_implant = new TH2D("PID_implant", "", 1e3, 1.95, 2.35, 1e3, 39, 56);
	PID = new TH2D("PID_correlated", "", 1e3, 1.95, 2.35, 1e3, 39, 56);

	globalEnergy = new TH1D("globalEnergy","",700,0,7000);
	isotopeSumEnergy = new TH1D("isotopeSumEnergy","",700,0,7000);

	Indium97_gammaveto_EdT = new TH2D("Indium97_gammaveto_EdT_ms","", 2000, -80, 80,700,0,7000);


	for (int i = 0; i < numElements; i++)
	{
		
		for (int j = 0; j <= isotopeEnd[i] - isotopeStart[i]; j++)
		{
		

			//use dssd arrays here
			for (int z = 0; z < 6; z++)
			{

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "ImplantBeta_DSSD" + std::to_string(z);
				implantBetaHis = new TH1D(hisName.c_str(), "", 2e2, -10, 10);
				implantBeta[i][z].push_back(implantBetaHis);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "ImplantBeta1p_DSSD" + std::to_string(z);
				implantBetaHis = new TH1D(hisName.c_str(), "", 2e2, -10, 10);
				implantBeta1p[i][z].push_back(implantBetaHis);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "Implant1p_DSSD" + std::to_string(z);
				implantBetaHis = new TH1D(hisName.c_str(), "", 2e2, -10, 10);
				implant1p[i][z].push_back(implantBetaHis);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "DecayEnergy_DSSD" + std::to_string(z);
				implantBetaHis = new TH1D(hisName.c_str(), "", 350, 0, 7000);
				decayEnergy[i][z].push_back(implantBetaHis);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "Delayed1pEnergy_DSSD" + std::to_string(z);
				implantBetaHis = new TH1D(hisName.c_str(), "", 350, 0, 7000);
				delayed1pEnergy[i][z].push_back(implantBetaHis);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "Delayed1pEnergyX_DSSD" + std::to_string(z);
				implantBetaHis = new TH1D(hisName.c_str(), "", 350, 0, 7000);
				delayed1pEnergyX[i][z].push_back(implantBetaHis);

				//hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "Delayed1pEnergyY_DSSD" + std::to_string(z);
				//implantBetaHis = new TH1D(hisName.c_str(), "", 350, 0, 7000);
				//delayed1pEnergyY[i][z].push_back(implantBetaHis);

				//hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "ExEyDiff_DSSD" + std::to_string(z);
				//implantBetaHis = new TH1D(hisName.c_str(), "", 1000, -500, 500);
				//ExEyDiff[i][z].push_back(implantBetaHis);

				//hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "Delayed1pEnergyTotal_DSSD" + std::to_string(z);
				//implantBetaHis = new TH1D(hisName.c_str(), "", 350, 0, 7000);
				//delayed1pEnergyAll[i][z].push_back(implantBetaHis);

				//hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "Delayed1pEnergyRandom_DSSD" + std::to_string(z);
				//implantBetaHis = new TH1D(hisName.c_str(), "", 350, 0, 7000);
				//delayed1pEnergyRandom[i][z].push_back(implantBetaHis);

				//hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "ImplantE_DSSD" + std::to_string(z);
				//implantBetaHis = new TH1D(hisName.c_str(), "", 350, 0, 7000);
				//implantE[i][z].push_back(implantBetaHis);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "EdT_DSSD" + std::to_string(z);
				implantBetaHis2D = new TH2D(hisName.c_str(), "", 2e2, -10, 10, 700, 0, 7000);
				EdT[i][z].push_back(implantBetaHis2D);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "EdT_ms_DSSD" + std::to_string(z);
				implantBetaHis2D = new TH2D(hisName.c_str(), "", 4000, -1000, 1000, 700, 0, 7000);
				EdT_ms[i][z].push_back(implantBetaHis2D);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "EdT_us_DSSD" + std::to_string(z);
				implantBetaHis2D = new TH2D(hisName.c_str(), "", 4000, -20000, 20000, 700, 0, 7000);
				EdT_us[i][z].push_back(implantBetaHis2D);

				//hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "ExEy_DSSD" + std::to_string(z);
				//implantBetaHis2D = new TH2D(hisName.c_str(), "", 280, 0, 7000, 280, 0, 7000);
				//ExEy[i][z].push_back(implantBetaHis2D);

				//hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "EvsXStrip_DSSD" + std::to_string(z);
				//implantBetaHis2D = new TH2D(hisName.c_str(), "", 128, 0, 127, 280, 0, 7000);
				//EnergyXChannel[i][z].push_back(implantBetaHis2D);

				//hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "EvsYStrip_DSSD" + std::to_string(z);
				//implantBetaHis2D = new TH2D(hisName.c_str(), "", 128, 0, 127, 280, 0, 7000);
				//EnergyYChannel[i][z].push_back(implantBetaHis2D);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "ImplantVelocityimplantE_DSSD" + std::to_string(z);
				implantBetaHis2D = new TH2D(hisName.c_str(), "", 100, 0.6, 0.65, 280, 0, 7000);
				implantVelocityimplantE[i][z].push_back(implantBetaHis2D);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "ImplantVelocityAOQ_DSSD" + std::to_string(z);
				implantBetaHis2D = new TH2D(hisName.c_str(), "", 100, 1.95, 2.15, 100, 0.6, 0.7);
				implantVelocityAOQ[i][z].push_back(implantBetaHis2D);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "ImplantEnergyAOQ_DSSD" + std::to_string(z);
				implantBetaHis2D = new TH2D(hisName.c_str(), "", 100, 1.95, 2.15, 280, 0, 7000);
				implantEnergyAOQ[i][z].push_back(implantBetaHis2D);
			}
			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "ImplantBeta_AllDSSD";
			implantBetaHis = new TH1D(hisName.c_str(), "", 2e2, -10, 10);
			implantBetaAll[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "ImplantBeta1p_AllDSSD";
			implantBetaHis = new TH1D(hisName.c_str(), "", 2e2, -10, 10);
			implantBeta1pAll[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "Implant1p_AllDSSD";
			implantBetaHis = new TH1D(hisName.c_str(), "", 2e2, -10, 10);
			implant1pAll[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "DecayEnergy_AllDSSD";
			implantBetaHis = new TH1D(hisName.c_str(), "", 350, 0, 7000);
			decayEnergyAll[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "Delayed1pEnergy_AllDSSD";
			implantBetaHis = new TH1D(hisName.c_str(), "", 350, 0, 7000);
			delayed1pEnergy_AllDSSD[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "Delayed1pEnergyTotal_AllDSSD";
			implantBetaHis = new TH1D(hisName.c_str(), "", 350, 0, 7000);
			delayed1pEnergyAll_AllDSSD[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "Delayed1pEnergyRandom_AllDSSD";
			implantBetaHis = new TH1D(hisName.c_str(), "", 350, 0, 7000);
			delayed1pEnergyRandom_AllDSSD[i].push_back(implantBetaHis);

			//hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "ImplantE_AllDSSD";
			//implantBetaHis = new TH1D(hisName.c_str(), "", 350, 0, 7000);
			//implantEAll[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "EdT_AllDSSD_NoMultiGate_s";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 2e2, -10, 10, 700, 0, 7000);
			EdTAll_NoMultiGate[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "EdT_AllDSSD_NoMultigate_ms";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 4000, -1000, 1000, 700, 0, 7000);
			EdTAll_ms[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "EdT_AllDSSD_NoMultigate_us";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 4000, -20000, 20000, 700, 0, 7000);
			EdTAll_us[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "EdT_AllDSSD_11_s";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 2e2, -10, 10, 700, 0, 7000);
			EdTAll11[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "EdT_AllDSSD_12_s";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 2e2, -10, 10, 700, 0, 7000);
			EdTAll12[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "EdT_AllDSSD_21_s";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 2e2, -10, 10, 700, 0, 7000);
			EdTAll21[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "EdT_AllDSSD_22_s";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 2e2, -10, 10, 700, 0, 7000);
			EdTAll22[i].push_back(implantBetaHis2D);

			//Cluster position checks

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "ExEy_AllDSSD_11";
			implantBetaHis2D = new TH2D(hisName.c_str(), "",700, 0, 7000, 700, 0, 7000);
			ExEy11[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "ExEy_AllDSSD_12";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 700, 0, 7000, 700, 0, 7000);
			ExEy12[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "ExEy_AllDSSD_21";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 700, 0, 7000, 700, 0, 7000);
			ExEy21[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "ExEy_AllDSSD_22";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 700, 0, 7000, 700, 0, 7000);
			ExEy22[i].push_back(implantBetaHis2D);

			//////////////////////////

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "ImplantVelocityimplantE_AllDSSD";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 100, 0.6, 0.65, 280, 0, 7000);
			implantVelocityimplantEAll[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "ImplantZ";
			implantBetaHis = new TH1D(hisName.c_str(), "", 6, 0, 6);
			implantZ[i].push_back(implantBetaHis);
//////////////////////////////////////***gamma_EdT***//////////////////////////////////////////////////

			/*hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "beta_gamma_EdT_s_corr";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 1e2, 0, 10, 700, 0, 7000);
			beta_gamma_EdT_s_corr[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "beta_gamma_EdT_ms_corr";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 2000, 0, 1000, 700, 0, 7000);
			beta_gamma_EdT_ms_corr[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "beta_gamma_EdT_us_corr";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 2000, 0, 20000, 700, 0, 7000);
			beta_gamma_EdT_us_corr[i].push_back(implantBetaHis2D);*/

			/*hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_beta_gamma_EdT_s_corr";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 1e2, 0, 10, 500, 0, 10000);
			summed_beta_gamma_EdT_s_corr[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_beta_gamma_EdT_ms_corr";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 2000, 0, 1000, 500, 0, 10000);
			summed_beta_gamma_EdT_ms_corr[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_beta_gamma_EdT_us_corr";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 2000, 0, 20000, 500, 0, 10000);
			summed_beta_gamma_EdT_us_corr[i].push_back(implantBetaHis2D);*/

			/*hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "bp_gamma_EdT_s_corr";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 1e2, 0, 10, 700, 0, 7000);
			bp_gamma_EdT_s_corr[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "bp_gamma_EdT_ms_corr";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 200, 0, 1000, 700, 0, 7000);
			bp_gamma_EdT_ms_corr[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "bp_gamma_EdT_us_corr";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 2000, 0, 20000, 700, 0, 7000);
			bp_gamma_EdT_us_corr[i].push_back(implantBetaHis2D);*/

			/*hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_bp_gamma_EdT_s_corr";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 1e2, 0, 10, 500, 0, 10000);
			summed_bp_gamma_EdT_s_corr[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_bp_gamma_EdT_ms_corr";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 200, 0, 1000, 500, 0, 10000);
			summed_bp_gamma_EdT_ms_corr[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_bp_gamma_EdT_us_corr";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 2000, 0, 20000, 500, 0, 10000);
			summed_bp_gamma_EdT_us_corr[i].push_back(implantBetaHis2D);*/
			
//////////////////////////////////////***gamma-gated histos****////////////////////////////////////////
			/*hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "EdT_gammagate";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 2000, -1000, 1000, 700, 0, 7000);
			EdT_gammagate[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "EdT_gammagate_longer";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 8000, -400, 400, 700, 0, 7000);
			EdT_gammagate_longer[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "ExEy_gammagate";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 700, 0, 7000, 700, 0, 7000);
			ExEy_gammagate[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "NxNy_gammagate";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 20, 0, 20, 20, 0, 20);
			NxNy_gammagate[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "clustersizes_gammaGate";
			implantBetaHis2D = new TH2D(hisName.c_str(), "",  20, 0, 20, 20, 0, 20);
			clustersize_gammagate[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "EDiff_gammaGate";
			implantBetaHis = new TH1D(hisName.c_str(), "",  1000, -500, 500);
			EDiff_gammagate[i].push_back(implantBetaHis);*/
			
			/////////////////////////////////////////********implant stuff**********//////////////////////////////////////////

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "ImplantVelocityimplantZ";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 100, 0.6, 0.65, 6, 0, 6);
			implantVelocityimplantZ[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "ImplantVelocityAOQ_All";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 100, 1.95, 2.15, 100, 0.6, 0.7);
			implantVelocityAOQ_AllDSSD[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "ImplantEnergyAOQ_All";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 100, 1.95, 2.15, 350, 0, 7000);
			implantEnergyAOQ_AllDSSD[i].push_back(implantBetaHis2D);

			//hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "ExEyDiff_All";
			//implantBetaHis = new TH1D(hisName.c_str(), "", 100, -500, 500);
			//ExEyDiffAll[i].push_back(implantBetaHis);

			//hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "ExEy_All";
			//implantBetaHis2D = new TH2D(hisName.c_str(), "", 700, 0, 7000, 700, 0, 7000);
			//ExEyAll[i].push_back(implantBetaHis2D);

			//hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "EvsXStrip_All";
			//implantBetaHis2D = new TH2D(hisName.c_str(), "", 128, 0, 127, 350, 0, 7000);
			//EnergyXChannelAll[i].push_back(implantBetaHis2D);

			//hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "EvsYStrip_All";
			//implantBetaHis2D = new TH2D(hisName.c_str(), "", 128, 0, 127, 350, 0, 7000);
			//EnergyYChannelAll[i].push_back(implantBetaHis2D);
		}
	}
}

#endif