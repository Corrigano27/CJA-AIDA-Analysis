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

std::vector<int> isotopeDSSDStart[4];
std::vector<int> isotopeDSSDEnd[4];

//graphical cuts for isotopes
TCutG *particleCuts[4][8]; //indices denote elements and isotopes

//directories
std::vector<TDirectory *> ElementDir;
std::vector<TDirectory *> IsotopeDir;
std::vector<TDirectory *> DSSD_Dir;

TDirectory *ElDir;
TDirectory *IsoDir;
TDirectory *ZDir;


//histogram vectors
//decays
std::vector<TH1D *> implantBeta[4][6];
std::vector<TH1D *> implantBeta1p[4][6];
std::vector<TH1D *> implant1p[4][6];
std::vector<TH1D *> decayEnergy[4][6];
std::vector<TH1D *> delayed1pEnergy[4][6];
std::vector<TH1D *> delayed1pEnergyX[4][6];

std::vector<TH1D *> delayed1pEnergyX_0_1_ms[4][6];
std::vector<TH1D *> delayed1pEnergyX_0_100_ms[4][6];
std::vector<TH1D *> delayed1pEnergyX_0_1_s[4][6];

std::vector<TH1D *> delayed1pEnergyY[4][6];
std::vector<TH1D *> delayed1pEnergyRandom[4][6];
std::vector<TH1D *> delayed1pEnergyAll[4][6];
std::vector<TH2D *> EdT[4][6];
std::vector<TH2D *> EdT_bp[4][6];
std::vector<TH2D *> implantVelocityimplantZ[4];
std::vector<TH2D *> implantVelocityimplantE[4][6];

std::vector<TH2D *> EnergyXChannel[4][6];
std::vector<TH2D *> EnergyYChannel[4][6];
std::vector<TH2D *> ExEy[4][6];
std::vector<TH1D *> ExEyDiff[4][6];

std::vector<TH1D *> implantBetaAll[4];
std::vector<TH1D *> implantBeta1pAll[4];
std::vector<TH1D *> implant1pAll[4];
std::vector<TH1D *> decayEnergyAll[4];
std::vector<TH1D *> delayed1pEnergy_AllDSSD[4];

std::vector<TH1D *> delayed1pEnergyRandom_AllDSSD[4];
std::vector<TH1D *> delayed1pEnergyAll_AllDSSD[4];
std::vector<TH2D *> EdTAll_NoMultiGate[4];
std::vector<TH2D *> EdTAll[4];
std::vector<TH2D *> EdTAll11[4];
std::vector<TH2D *> EdTAll12[4];
std::vector<TH2D *> EdTAll21[4];
std::vector<TH2D *> EdTAll22[4];
std::vector<TH1D *> implantEAll[4];
std::vector<TH2D *> implantVelocityimplantEAll[4];

//all DSSD, gated on 511 summed spectra
std::vector<TH2D *> ExEyAll_gammaloop[4];
std::vector<TH2D *> EdT_gammagate[4];
std::vector<TH2D *> EdT_gammagate_longer[4];
std::vector<TH2D *> ExEy_gammagate[4];
std::vector<TH2D *> NxNy_gammagate[4];
std::vector<TH2D *> clustersize_gammagate[4];
std::vector<TH1D *> EDiff_gammagate[4];

//correcting beta-gamma singles
std::vector<TH1D *> beta_gamma_1[4];
std::vector<TH1D *> beta_gamma_2[4];
std::vector<TH1D *> beta_gamma_3[4];
std::vector<TH1D *> beta_gamma_4[4];

std::vector<TH1D *> beta_gamma_corr[4];

//correcting beta-gamma summed
std::vector<TH1D *> summed_beta_gamma_1[4];
std::vector<TH1D *> summed_beta_gamma_2[4];
std::vector<TH1D *> summed_beta_gamma_3[4];
std::vector<TH1D *> summed_beta_gamma_4[4];

std::vector<TH1D *> summed_beta_gamma_corr[4];

//correcting bp-gamma singles
std::vector<TH1D *> bp_gamma_1[4];
std::vector<TH1D *> bp_gamma_2[4];
std::vector<TH1D *> bp_gamma_3[4];
std::vector<TH1D *> bp_gamma_4[4];

std::vector<TH1D *> bp_gamma_corr[4];

//correcting bp-gamma summed
std::vector<TH1D *> summed_bp_gamma_1[4];
std::vector<TH1D *> summed_bp_gamma_2[4];
std::vector<TH1D *> summed_bp_gamma_3[4];
std::vector<TH1D *> summed_bp_gamma_4[4];

std::vector<TH1D *> summed_bp_gamma_corr[4];

std::vector<TH2D *> EnergyXChannelAll[4];
std::vector<TH2D *> EnergyYChannelAll[4];
std::vector<TH2D *> ExEyAll[4];
std::vector<TH1D *> ExEyDiffAll[4];

std::vector<TH2D *> implantVelocityAOQ[4][6];
std::vector<TH2D *> implantVelocityAOQ_AllDSSD[4];

std::vector<TH2D *> implantEnergyAOQ[4][6];
std::vector<TH2D *> implantEnergyAOQ_AllDSSD[4];

//implants
std::vector<TH1D *> implantZ[4];
std::vector<TH1D *> implantE[4][6];


//template histograms

TH1D *implantBetaHis;
TH2D *implantBetaHis2D;

//pid
TH2D *PID;

TH2D *PID_implant;

//global no gamma gate
TH2D *EdT_global;

TH2D *EdT_global_longer;

TH1D *EDiff_global;

TH2D *ExEy_global;

//global gamma gate

TH2D *EdT_global_longer_gammagate;

TH2D *EdT_global_gammagate;

TH1D *EDiff_global_gammagate;

TH2D *ExEy_global_gammagate;

TH2D *NxNy_global_gammagate;

TH2D *CxCy_global_gammagate;


//101SnSplitting
/*

TH1D *Sn101ImplantBeta1p_DSSD1_smallpeak;

TH1D *Sn101ImplantBeta1p_DSSD1_largepeak;

TH1D *Sn101_peak_GammaSingle;

TH1D *Sn101_peak_Gamma777;

TH1D *Sn101_peak_GammaSingle_Bg;

TH1D *Sn101_peak_Gamma777_Bg;

TH1D *Ag94ImplantBeta1p_DSSD2_smallpeak;

TH1D *Ag94ImplantBeta1p_DSSD2_largepeak;

TH1D *Ag94_peak_GammaSingle;

TH1D *Ag94_peak_Gamma777;

TH1D *Ag94_peak_GammaSingle_Bg;

TH1D *Ag94_peak_Gamma777_Bg;


//In97m gamma

TH1D *In97m_Gamma777;

TH1D *In97m_GammaSingle;

TH1D *In97m_Gamma777_Bg;

TH1D *In97m_GammaSingle_Bg;

TH1D *In97_GroundStateE;
*/

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

	PID = new TH2D("PID", "", 1e3, 1.95, 2.35, 1e3, 39, 56);

	PID_implant = new TH2D("PID_implant", "", 1e3, 1.95, 2.35, 1e3, 39, 56);

	/*Sn101ImplantBeta1p_DSSD1_smallpeak = new TH1D("Sn101ImplantBeta1p_DSSD1_smallpeak", "", 2e2, -10, 10);

	Sn101ImplantBeta1p_DSSD1_largepeak = new TH1D("Sn101ImplantBeta1p_DSSD1_largepeak", "", 2e2, -10, 10);

	Sn101_peak_GammaSingle = new TH1D("Sn101_peak_GammaSingle","",3500,0,7000);

	Sn101_peak_Gamma777 = new TH1D("Sn101_peak_Gamma777","",3500,0,7000);

	Sn101_peak_GammaSingle_Bg = new TH1D("Sn101_peak_GammaSingle_Bg","",3500,0,7000);

	Sn101_peak_Gamma777_Bg = new TH1D("Sn101_peak_Gamma777_Bg","",3500,0,7000);

	In97m_GammaSingle = new TH1D("In97m_GammaSingle","",3500,0,7000);

	In97m_Gamma777 = new TH1D("In97m_Gamma777","",3500,0,7000);

	In97m_GammaSingle_Bg = new TH1D("In97m_GammaSingle_Bg","",3500,0,7000);

	In97m_Gamma777_Bg = new TH1D("In97m_Gamma777_Bg","",3500,0,7000);

	In97_GroundStateE = new TH1D("In97_GroundStateE","",250,0,5000);	

	Ag94ImplantBeta1p_DSSD2_smallpeak = new TH1D("Ag94ImplantBeta1p_DSSD2_smallpeak", "", 2e2, -10, 10);

	Ag94ImplantBeta1p_DSSD2_largepeak = new TH1D("Ag94ImplantBeta1p_DSSD2_largepeak", "", 2e2, -10, 10);

	Ag94_peak_GammaSingle = new TH1D("Ag94_peak_GammaSingle","",3500,0,7000);

	Ag94_peak_Gamma777 = new TH1D("Ag94_peak_Gamma777","",3500,0,7000);

	Ag94_peak_GammaSingle_Bg = new TH1D("Ag94_peak_GammaSingle_Bg","",3500,0,7000);

	Ag94_peak_Gamma777_Bg = new TH1D("Ag94_peak_Gamma777_Bg","",3500,0,7000);*/

	EdT_global_longer_gammagate = new TH2D("EdT_global_longer_gammagate","",8000, -400, 400, 700, 0, 7000);

	EdT_global_gammagate = new TH2D("EdT_global_gammagate","",2000,-10000,10000,700,0,7000);

	EDiff_global_gammagate = new TH1D("EDiff_global_gammagate","",1000,-500,500);

	NxNy_global_gammagate = new TH2D("NxNy_global_gammagate","",20,0,20,20,0,20); 

	CxCy_global_gammagate = new TH2D("CxCy_global_gammagate","",20,0,20,20,0,20);

	EdT_global_longer = new TH2D("EdT_global_longer","",8000, -400, 400, 700, 0, 7000);

	EdT_global = new TH2D("EdT_global","",2000,-10000,10000,700,0,7000);

	EDiff_global = new TH1D("EDiff_global","",1000,-500,500);

	ExEy_global = new TH2D("ExEy_global","",1000,0,20000,1000,0,20000);

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

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "Delayed1pEnergyX_0_1_ms_DSSD" + std::to_string(z);
				implantBetaHis = new TH1D(hisName.c_str(), "", 150, 0, 3000);
				delayed1pEnergyX_0_1_ms[i][z].push_back(implantBetaHis);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "Delayed1pEnergyX_0_100_ms_DSSD" + std::to_string(z);
				implantBetaHis = new TH1D(hisName.c_str(), "", 150, 0, 3000);
				delayed1pEnergyX_0_100_ms[i][z].push_back(implantBetaHis);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "Delayed1pEnergyX_0_1_s_DSSD" + std::to_string(z);
				implantBetaHis = new TH1D(hisName.c_str(), "", 150, 0, 3000);
				delayed1pEnergyX_0_1_s[i][z].push_back(implantBetaHis);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "Delayed1pEnergyY_DSSD" + std::to_string(z);
				implantBetaHis = new TH1D(hisName.c_str(), "", 350, 0, 7000);
				delayed1pEnergyY[i][z].push_back(implantBetaHis);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "ExEyDiff_DSSD" + std::to_string(z);
				implantBetaHis = new TH1D(hisName.c_str(), "", 1000, -500, 500);
				ExEyDiff[i][z].push_back(implantBetaHis);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "Delayed1pEnergyTotal_DSSD" + std::to_string(z);
				implantBetaHis = new TH1D(hisName.c_str(), "", 350, 0, 7000);
				delayed1pEnergyAll[i][z].push_back(implantBetaHis);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "Delayed1pEnergyRandom_DSSD" + std::to_string(z);
				implantBetaHis = new TH1D(hisName.c_str(), "", 350, 0, 7000);
				delayed1pEnergyRandom[i][z].push_back(implantBetaHis);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "ImplantE_DSSD" + std::to_string(z);
				implantBetaHis = new TH1D(hisName.c_str(), "", 350, 0, 7000);
				implantE[i][z].push_back(implantBetaHis);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "EdT_DSSD" + std::to_string(z);
				implantBetaHis2D = new TH2D(hisName.c_str(), "", 200, -2000, 2000, 280, 0, 7000);
				EdT[i][z].push_back(implantBetaHis2D);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "EdT_bp_DSSD" + std::to_string(z);
				implantBetaHis2D = new TH2D(hisName.c_str(), "", 4000, -10000, 10000, 280, 0, 7000);
				EdT_bp[i][z].push_back(implantBetaHis2D);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "ExEy_DSSD" + std::to_string(z);
				implantBetaHis2D = new TH2D(hisName.c_str(), "", 280, 0, 7000, 280, 0, 7000);
				ExEy[i][z].push_back(implantBetaHis2D);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "EvsXStrip_DSSD" + std::to_string(z);
				implantBetaHis2D = new TH2D(hisName.c_str(), "", 128, 0, 127, 280, 0, 7000);
				EnergyXChannel[i][z].push_back(implantBetaHis2D);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "EvsYStrip_DSSD" + std::to_string(z);
				implantBetaHis2D = new TH2D(hisName.c_str(), "", 128, 0, 127, 280, 0, 7000);
				EnergyYChannel[i][z].push_back(implantBetaHis2D);

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

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "ImplantE_AllDSSD";
			implantBetaHis = new TH1D(hisName.c_str(), "", 350, 0, 7000);
			implantEAll[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "EdT_AllDSSD_NoMultiGate";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 8000, -400, 400, 700, 0, 7000);
			EdTAll_NoMultiGate[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "EdT_All_ShortLifetime";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 2000, -1000, 1000, 700, 0, 7000);
			EdTAll[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "EdT_AllDSSD_11";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 4000, -2000, 2000, 700, 0, 7000);
			EdTAll11[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "EdT_AllDSSD_12";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 4000, -2000, 2000, 700, 0, 7000);
			EdTAll12[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "EdT_AllDSSD_21";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 4000, -2000, 2000, 700, 0, 7000);
			EdTAll21[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "EdT_AllDSSD_22";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 4000, -2000, 2000, 700, 0, 7000);
			EdTAll22[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "ImplantVelocityimplantE_AllDSSD";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 100, 0.6, 0.65, 280, 0, 7000);
			implantVelocityimplantEAll[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "ImplantZ";
			implantBetaHis = new TH1D(hisName.c_str(), "", 6, 0, 6);
			implantZ[i].push_back(implantBetaHis);

//////////////////////////////////////***gamma-gated histos****////////////////////////////////////////
			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "EdT_gammagate";
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
			EDiff_gammagate[i].push_back(implantBetaHis);

//////////////////////////////////////***beta-gammas***///////////////////////////////////////////////
			//singles
			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "beta_gamma_1";
			implantBetaHis = new TH1D(hisName.c_str(), "", 3500,0,7000);
			beta_gamma_1[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "beta_gamma_2";
			implantBetaHis = new TH1D(hisName.c_str(), "", 3500,0,7000);
			beta_gamma_2[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "beta_gamma_3";
			implantBetaHis = new TH1D(hisName.c_str(), "", 3500,0,7000);
			beta_gamma_3[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "beta_gamma_4";
			implantBetaHis = new TH1D(hisName.c_str(), "", 3500,0,7000);
			beta_gamma_4[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "beta_gamma_corr";
			implantBetaHis = new TH1D(hisName.c_str(), "", 3500,0,7000);
			beta_gamma_corr[i].push_back(implantBetaHis);

			//summed
			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_beta_gamma_1";
			implantBetaHis = new TH1D(hisName.c_str(), "", 3500,0,7000);
			summed_beta_gamma_1[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_beta_gamma_2";
			implantBetaHis = new TH1D(hisName.c_str(), "", 3500,0,7000);
			summed_beta_gamma_2[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_beta_gamma_3";
			implantBetaHis = new TH1D(hisName.c_str(), "", 3500,0,7000);
			summed_beta_gamma_3[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_beta_gamma_4";
			implantBetaHis = new TH1D(hisName.c_str(), "", 3500,0,7000);
			summed_beta_gamma_4[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_beta_gamma_corr";
			implantBetaHis = new TH1D(hisName.c_str(), "", 3500,0,7000);
			summed_beta_gamma_corr[i].push_back(implantBetaHis);

///////////////////////////////////*****beta-delayed proton gammas******/////////////////////////////////////////////////////
			//singles
			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "bp_gamma_1";
			implantBetaHis = new TH1D(hisName.c_str(), "", 3500,0,7000);
			bp_gamma_1[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "bp_gamma_2";
			implantBetaHis = new TH1D(hisName.c_str(), "", 3500,0,7000);
			bp_gamma_2[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "bp_gamma_3";
			implantBetaHis = new TH1D(hisName.c_str(), "", 3500,0,7000);
			bp_gamma_3[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "bp_gamma_4";
			implantBetaHis = new TH1D(hisName.c_str(), "", 3500,0,7000);
			bp_gamma_4[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "bp_gamma_corr";
			implantBetaHis = new TH1D(hisName.c_str(), "", 3500,0,7000);
			bp_gamma_corr[i].push_back(implantBetaHis);

			//summed
			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_bp_gamma_1";
			implantBetaHis = new TH1D(hisName.c_str(), "", 3500,0,7000);
			summed_bp_gamma_1[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_bp_gamma_2";
			implantBetaHis = new TH1D(hisName.c_str(), "", 3500,0,7000);
			summed_bp_gamma_2[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_bp_gamma_3";
			implantBetaHis = new TH1D(hisName.c_str(), "", 3500,0,7000);
			summed_bp_gamma_3[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_bp_gamma_4";
			implantBetaHis = new TH1D(hisName.c_str(), "", 3500,0,7000);
			summed_bp_gamma_4[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_bp_gamma_corr";
			implantBetaHis = new TH1D(hisName.c_str(), "", 3500,0,7000);
			summed_bp_gamma_corr[i].push_back(implantBetaHis);

			///////////////////////////////////////////////////////////////////////////////////

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "ImplantVelocityimplantZ";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 100, 0.6, 0.65, 6, 0, 6);
			implantVelocityimplantZ[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "ImplantVelocityAOQ_All";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 100, 1.95, 2.15, 100, 0.6, 0.7);
			implantVelocityAOQ_AllDSSD[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "ImplantEnergyAOQ_All";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 100, 1.95, 2.15, 350, 0, 7000);
			implantEnergyAOQ_AllDSSD[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "ExEyDiff_All";
			implantBetaHis = new TH1D(hisName.c_str(), "", 100, -500, 500);
			ExEyDiffAll[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "ExEy_All";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 700, 0, 7000, 700, 0, 7000);
			ExEyAll[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "EvsXStrip_All";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 128, 0, 127, 350, 0, 7000);
			EnergyXChannelAll[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "EvsYStrip_All";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 128, 0, 127, 350, 0, 7000);
			EnergyYChannelAll[i].push_back(implantBetaHis2D);
		}
	}
}

#endif