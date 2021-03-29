#ifndef PARTICLECUTSSN100_PDRHRU_H
#define PARTICLECUTSSN100_PDRHRU_H

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

int numElements = 3;

//string declarations;
std::string particleCut; //PID cuts from file

std::string elements[3]; //isotope elements


//vectors of isotope mass numbers
int isotopeStart[3];
int isotopeEnd[3];

std::vector<int> isotopeDSSDStart[3];
std::vector<int> isotopeDSSDEnd[3];

//graphical cuts for isotopes
TCutG *particleCuts[3][8]; //indices denote elements and isotopes

//directories
std::vector<TDirectory *> ElementDir;
std::vector<TDirectory *> IsotopeDir;
std::vector<TDirectory *> DSSD_Dir;

TDirectory *ElDir;
TDirectory *IsoDir;
TDirectory *ZDir;


//histogram vectors
//decays
std::vector<TH1D *> implantBeta[3][6];
std::vector<TH1D *> implantBeta1p[3][6];
std::vector<TH1D *> implant1p[3][6];
std::vector<TH1D *> decayEnergy[3][6];
std::vector<TH1D *> delayed1pEnergy[3][6];
std::vector<TH1D *> delayed1pEnergyX[3][6];

std::vector<TH1D *> delayed1pEnergyX_0_1_ms[3][6];
std::vector<TH1D *> delayed1pEnergyX_0_100_ms[3][6];
std::vector<TH1D *> delayed1pEnergyX_0_1_s[3][6];

std::vector<TH1D *> delayed1pEnergyY[3][6];
std::vector<TH1D *> delayed1pEnergyRandom[3][6];
std::vector<TH1D *> delayed1pEnergyAll[3][6];
std::vector<TH2D *> EdT[3][6];
std::vector<TH2D *> EdT_bp[3][6];
std::vector<TH2D *> implantVelocityimplantZ[3];
std::vector<TH2D *> implantVelocityimplantE[3][6];

std::vector<TH2D *> EnergyXChannel[3][6];
std::vector<TH2D *> EnergyYChannel[3][6];
std::vector<TH2D *> ExEy[3][6];
std::vector<TH1D *> ExEyDiff[3][6];

std::vector<TH1D *> implantBetaAll[3];
std::vector<TH1D *> implantBeta1pAll[3];
std::vector<TH1D *> implant1pAll[3];
std::vector<TH1D *> decayEnergyAll[3];
std::vector<TH1D *> delayed1pEnergy_AllDSSD[3];

std::vector<TH1D *> delayed1pEnergyRandom_AllDSSD[3];
std::vector<TH1D *> delayed1pEnergyAll_AllDSSD[3];
std::vector<TH2D *> EdTAll_NoMultiGate[3];
std::vector<TH2D *> EdTAll_ms[3];
std::vector<TH2D *> EdTAll_us[3];
std::vector<TH2D *> EdTAll11[3];
std::vector<TH2D *> EdTAll12[3];
std::vector<TH2D *> EdTAll21[3];
std::vector<TH2D *> EdTAll22[3];
std::vector<TH1D *> implantEAll[3];
std::vector<TH2D *> implantVelocityimplantEAll[3];

//all DSSD, gated on 511 summed spectra
std::vector<TH2D *> ExEyAll_gammaloop[3];
std::vector<TH2D *> EdT_gammagate[3];
std::vector<TH2D *> EdT_gammagate_longer[3];
std::vector<TH2D *> ExEy_gammagate[3];
std::vector<TH2D *> NxNy_gammagate[3];
std::vector<TH2D *> clustersize_gammagate[3];
std::vector<TH1D *> EDiff_gammagate[3];

//correcting beta-gamma singles
std::vector<TH1D *> beta_gamma_1[4];
std::vector<TH1D *> beta_gamma_2[4];
std::vector<TH1D *> beta_gamma_3[4];
std::vector<TH1D *> beta_gamma_4[4];

std::vector<TH1D *> beta_gamma_corr[4];

std::vector<TH2D *> beta_gamma_EdT_s[3][4];
std::vector<TH2D *> beta_gamma_EdT_ms[3][4];
std::vector<TH2D *> beta_gamma_EdT_us[3][4];

std::vector<TH2D *> beta_gamma_EdT_s_corr[4];
std::vector<TH2D *> beta_gamma_EdT_ms_corr[4];
std::vector<TH2D *> beta_gamma_EdT_us_corr[4];

//correcting beta-gamma summed
std::vector<TH1D *> summed_beta_gamma_1[4];
std::vector<TH1D *> summed_beta_gamma_2[4];
std::vector<TH1D *> summed_beta_gamma_3[4];
std::vector<TH1D *> summed_beta_gamma_4[4];

std::vector<TH1D *> summed_beta_gamma_corr[4];

std::vector<TH2D *> summed_beta_gamma_EdT_s[3][4];
std::vector<TH2D *> summed_beta_gamma_EdT_ms[3][4];
std::vector<TH2D *> summed_beta_gamma_EdT_us[3][4];

std::vector<TH2D *> summed_beta_gamma_EdT_s_corr[4];
std::vector<TH2D *> summed_beta_gamma_EdT_ms_corr[4];
std::vector<TH2D *> summed_beta_gamma_EdT_us_corr[4];


//correcting bp-gamma singles
std::vector<TH1D *> bp_gamma_1[4];
std::vector<TH1D *> bp_gamma_2[4];
std::vector<TH1D *> bp_gamma_3[4];
std::vector<TH1D *> bp_gamma_4[4];

std::vector<TH1D *> bp_gamma_corr[4];

std::vector<TH2D *> bp_gamma_EdT_s[3][4];
std::vector<TH2D *> bp_gamma_EdT_ms[3][4];
std::vector<TH2D *> bp_gamma_EdT_us[3][4];

std::vector<TH2D *> bp_gamma_EdT_s_corr[4];
std::vector<TH2D *> bp_gamma_EdT_ms_corr[4];
std::vector<TH2D *> bp_gamma_EdT_us_corr[4];


//correcting bp-gamma summed
std::vector<TH1D *> summed_bp_gamma_1[4];
std::vector<TH1D *> summed_bp_gamma_2[4];
std::vector<TH1D *> summed_bp_gamma_3[4];
std::vector<TH1D *> summed_bp_gamma_4[4];

std::vector<TH1D *> summed_bp_gamma_corr[4];

std::vector<TH2D *> summed_bp_gamma_EdT_s[3][4];
std::vector<TH2D *> summed_bp_gamma_EdT_ms[3][4];
std::vector<TH2D *> summed_bp_gamma_EdT_us[3][4];

std::vector<TH2D *> summed_bp_gamma_EdT_s_corr[4];
std::vector<TH2D *> summed_bp_gamma_EdT_ms_corr[4];
std::vector<TH2D *> summed_bp_gamma_EdT_us_corr[4];

/*
//Sn101 specific
TH1D *Tin101_bp_gamma_peak[4];
TH1D *Tin101_bp_gamma_peak_corr;
TH1D *Tin101_bp_gamma_rest[4];
TH1D *Tin101_bp_gamma_rest_corr;
TH1D *Tin101_summed_bp_gamma_peak[4];
TH1D *Tin101_summed_bp_gamma_peak_corr;
TH1D *Tin101_summed_bp_gamma_rest[4];
TH1D *Tin101_summed_bp_gamma_rest_corr;
*/


//end of gammas

std::vector<TH2D *> EnergyXChannelAll[3];
std::vector<TH2D *> EnergyYChannelAll[3];
std::vector<TH2D *> ExEyAll[3];
std::vector<TH1D *> ExEyDiffAll[3];

std::vector<TH2D *> implantVelocityAOQ[3][6];
std::vector<TH2D *> implantVelocityAOQ_AllDSSD[3];

std::vector<TH2D *> implantEnergyAOQ[3][6];
std::vector<TH2D *> implantEnergyAOQ_AllDSSD[3];

//implants
std::vector<TH1D *> implantZ[3];
std::vector<TH1D *> implantE[3][6];


//template histograms

TH1D *implantBetaHis;
TH2D *implantBetaHis2D;

//pid
//TH2D *PID;

//TH2D *PID_noise;

TH2D *PID_implant;

//TH2D *PID_gamma;
/*
//global no gamma gate
TH2D *EdT_global;

TH2D *EdT_global_longer;

TH1D *EDiff_global;

TH2D *ExEy_global;

TH2D *XY_Hits;

TH1D *FastDSSD;

//global gamma gate

TH2D *EdT_global_longer_gammagate;

TH2D *EdT_global_gammagate;

TH1D *EDiff_global_gammagate;

TH2D *ExEy_global_gammagate;

TH2D *XY_Hits_gammagate;

TH1D *FastDSSD_gammagate;

TH2D *NxNy_global_gammagate;

TH2D *CxCy_global_gammagate;


//101SnSplitting

TH1D *Sn101_EdT;

TH1D *Sn101_peak_GammaSingle;

TH1D *Sn101_peak_Gamma777;

TH1D *Sn101_peak_GammaSingle_Bg;

TH1D *Sn101_peak_Gamma777_Bg;
*/

//TH2D *Indium97_gammaveto_EdT;

//range of isotopes for each element
void SetParticles()
{

	elements[0] = "Pd";
	isotopeStart[0] = 93;
	isotopeEnd[0] = 99;

	elements[1] = "Rh";
	isotopeStart[1] = 92;
	isotopeEnd[1] = 98;

	elements[2] = "Ru";
	isotopeStart[2] = 91;
	isotopeEnd[2] = 97;

}

//Set DSSD layers of expected implantation for each isotope
void SetImplantDSSD()
{

	//Pd isotopes 93-99
	isotopeDSSDStart[0].push_back(2);
	isotopeDSSDEnd[0].push_back(3);

	isotopeDSSDStart[0].push_back(1);
	isotopeDSSDEnd[0].push_back(3);

	isotopeDSSDStart[0].push_back(1);
	isotopeDSSDEnd[0].push_back(2);

	isotopeDSSDStart[0].push_back(1);
	isotopeDSSDEnd[0].push_back(2);

	isotopeDSSDStart[0].push_back(1);
	isotopeDSSDEnd[0].push_back(2);

	isotopeDSSDStart[0].push_back(1);
	isotopeDSSDEnd[0].push_back(2);

	isotopeDSSDStart[0].push_back(1);
	isotopeDSSDEnd[0].push_back(2);

	//Rh isotopes 92-98
	isotopeDSSDStart[1].push_back(2);
	isotopeDSSDEnd[1].push_back(3);

	isotopeDSSDStart[1].push_back(2);
	isotopeDSSDEnd[1].push_back(3);

	isotopeDSSDStart[1].push_back(1);
	isotopeDSSDEnd[1].push_back(2);

	isotopeDSSDStart[1].push_back(1);
	isotopeDSSDEnd[1].push_back(2);

	isotopeDSSDStart[1].push_back(1);
	isotopeDSSDEnd[1].push_back(2);

	isotopeDSSDStart[1].push_back(1);
	isotopeDSSDEnd[1].push_back(2);

	isotopeDSSDStart[1].push_back(1);
	isotopeDSSDEnd[1].push_back(2);

	//Ru isotopes 91-97
	isotopeDSSDStart[2].push_back(2);
	isotopeDSSDEnd[2].push_back(3);

	isotopeDSSDStart[2].push_back(1);
	isotopeDSSDEnd[2].push_back(3);

	isotopeDSSDStart[2].push_back(1);
	isotopeDSSDEnd[2].push_back(2);

	isotopeDSSDStart[2].push_back(1);
	isotopeDSSDEnd[2].push_back(2);

	isotopeDSSDStart[2].push_back(1);
	isotopeDSSDEnd[2].push_back(2);

	isotopeDSSDStart[2].push_back(1);
	isotopeDSSDEnd[2].push_back(2);

	isotopeDSSDStart[2].push_back(1);
	isotopeDSSDEnd[2].push_back(2);

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

	//PID = new TH2D("PID", "", 1e3, 1.95, 2.35, 1e3, 39, 56);

	//PID_noise = new TH2D("PID_noise", "", 1e3, 1.95, 2.35, 1e3, 39, 56);

	//PID_gamma = new TH2D("PID_gamma", "", 1e3, 1.95, 2.35, 1e3, 39, 56);

	PID_implant = new TH2D("PID_implant", "", 1e3, 1.95, 2.35, 1e3, 39, 56);
/*
	Sn101_peak_GammaSingle = new TH1D("Sn101_peak_GammaSingle","",3500,0,7000);

	Sn101_peak_Gamma777 = new TH1D("Sn101_peak_Gamma777","",3500,0,7000);

	Sn101_peak_GammaSingle_Bg = new TH1D("Sn101_peak_GammaSingle_Bg","",3500,0,7000);

	Sn101_peak_Gamma777_Bg = new TH1D("Sn101_peak_Gamma777_Bg","",3500,0,7000);	

	EdT_global_longer_gammagate = new TH2D("EdT_global_longer_gammagate","",8000, -400, 400, 700, 0, 7000);

	EdT_global_gammagate = new TH2D("EdT_global_gammagate","",2000,-10000,10000,700,0,7000);

	EDiff_global_gammagate = new TH1D("Ag95EDiff_global_gammagate","",1000,-500,500);

	NxNy_global_gammagate = new TH2D("Ag95NxNy_global_gammagate","",20,0,20,20,0,20); 

	CxCy_global_gammagate = new TH2D("Ag95CxCy_global_gammagate","",20,0,20,20,0,20);

	EdT_global_longer = new TH2D("EdT_global_longer","",8000, -400, 400, 700, 0, 7000);

	EdT_global = new TH2D("EdT_global","",2000,-10000,10000,700,0,7000);

	EDiff_global = new TH1D("EDiff_global","",1000,-500,500);

	ExEy_global = new TH2D("ExEy_global","",1000,0,20000,1000,0,20000);

	XY_Hits = new TH2D("XY_Hits","",128,0,128,128,0,128);

	XY_Hits_gammagate = new TH2D("Ag95XY_Hits_gammagate","",128,0,128,128,0,128);

	FastDSSD = new TH1D("FastDSSD","",6,0,5);

	FastDSSD_gammagate = new TH1D("Ag95FastDSSD_gammagate","",6,0,5); */

	//Indium97_gammaveto_EdT = new TH2D("Indium97_gammaveto_EdT_us","", 2000, -10, 10,700,0,7000);

	/*
	//Sn101
	for (int g = 0; g < 4; g++){
			
		Sn101histName = "Tin101_bp_gamma_peak_" + std::to_string(g);
		Tin101_bp_gamma_peak[g] = new TH1D(hisName.c_str(), "", 1000, 0, 7000);

		Sn101histName = "Tin101_summed_bp_gamma_peak_" + std::to_string(g);
		Tin101_summed_bp_gamma_peak[g] = new TH1D(hisName.c_str(), "", 1000, 0, 7000);

		Sn101histName = "Tin101_bp_gamma_rest_" + std::to_string(g);
		Tin101_bp_gamma_rest[g] = new TH1D(hisName.c_str(), "", 1000, 0, 7000);

		Sn101histName = "Tin101_summed_bp_gamma_rest_" + std::to_string(g);
		Tin101_summed_bp_gamma_rest[g] = new TH1D(hisName.c_str(), "", 1000, 0, 7000);

	}
	
	Tin101_bp_gamma_peak_corr = new TH1D("Tin101_bp_gamma_peak_corr","",1000, 0, 7000);
	Tin101_summed_bp_gamma_peak_corr = new TH1D("Tin101_summed_bp_gamma_peak_corr","",1000, 0, 7000);
	Tin101_bp_gamma_rest_corr = new TH1D("Tin101_bp_gamma_rest_corr","",1000, 0, 7000);
	Tin101_summed_bp_gamma_rest_corr = new TH1D("Tin101_summed_bp_gamma_rest_corr","",1000, 0, 7000);

	*/

	for (int i = 0; i < numElements; i++)
	{
		
		for (int j = 0; j <= isotopeEnd[i] - isotopeStart[i]; j++)
		{
			//EdT_gamma_histograms
			for (int g = 0; g < 4; g++){
			
				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "beta_gamma_EdT_s_" + std::to_string(g);
				implantBetaHis2D = new TH2D(hisName.c_str(), "", 1e2, 0, 10, 1000, 0, 7000);
				beta_gamma_EdT_s[i][g].push_back(implantBetaHis2D);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "beta_gamma_EdT_ms_" + std::to_string(g);
				implantBetaHis2D = new TH2D(hisName.c_str(), "", 100, 0, 1000, 1000, 0, 7000);
				beta_gamma_EdT_ms[i][g].push_back(implantBetaHis2D);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "beta_gamma_EdT_us_" + std::to_string(g);
				implantBetaHis2D = new TH2D(hisName.c_str(), "", 1000, 0, 10000, 1000, 0, 7000);
				beta_gamma_EdT_us[i][g].push_back(implantBetaHis2D);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_beta_gamma_EdT_s_" + std::to_string(g);
				implantBetaHis2D = new TH2D(hisName.c_str(), "", 1e2, 0, 10, 1000, 0, 7000);
				summed_beta_gamma_EdT_s[i][g].push_back(implantBetaHis2D);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_beta_gamma_EdT_ms_" + std::to_string(g);
				implantBetaHis2D = new TH2D(hisName.c_str(), "", 100, 0, 1000, 1000, 0, 7000);
				summed_beta_gamma_EdT_ms[i][g].push_back(implantBetaHis2D);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_beta_gamma_EdT_us_" + std::to_string(g);
				implantBetaHis2D = new TH2D(hisName.c_str(), "", 1000, 0, 10000, 1000, 0, 7000);
				summed_beta_gamma_EdT_us[i][g].push_back(implantBetaHis2D);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "bp_gamma_EdT_s_" + std::to_string(g);
				implantBetaHis2D = new TH2D(hisName.c_str(), "", 1e2, 0, 10, 1000, 0, 7000);
				bp_gamma_EdT_s[i][g].push_back(implantBetaHis2D);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "bp_gamma_EdT_ms_" + std::to_string(g);
				implantBetaHis2D = new TH2D(hisName.c_str(), "", 100, 0, 1000, 1000, 0, 7000);
				bp_gamma_EdT_ms[i][g].push_back(implantBetaHis2D);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "bp_gamma_EdT_us_" + std::to_string(g);
				implantBetaHis2D = new TH2D(hisName.c_str(), "", 1000, 0, 10000, 1000, 0, 7000);
				bp_gamma_EdT_us[i][g].push_back(implantBetaHis2D);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_bp_gamma_EdT_s_" + std::to_string(g);
				implantBetaHis2D = new TH2D(hisName.c_str(), "", 1e2, 0, 10, 1000, 0, 7000);
				summed_bp_gamma_EdT_s[i][g].push_back(implantBetaHis2D);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_bp_gamma_EdT_ms_" + std::to_string(g);
				implantBetaHis2D = new TH2D(hisName.c_str(), "", 100, 0, 1000, 1000, 0, 7000);
				summed_bp_gamma_EdT_ms[i][g].push_back(implantBetaHis2D);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_bp_gamma_EdT_us_" + std::to_string(g);
				implantBetaHis2D = new TH2D(hisName.c_str(), "", 1000, 0, 10000, 1000, 0, 7000);
				summed_bp_gamma_EdT_us[i][g].push_back(implantBetaHis2D);

			}

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

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "EdT_AllDSSD_NoMultiGate_s";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 2e2, -10, 10, 700, 0, 7000);
			EdTAll_NoMultiGate[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "EdT_AllDSSD_NoMultigate_ms";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 200, -1000, 1000, 700, 0, 7000);
			EdTAll_ms[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "EdT_AllDSSD_NoMultigate_us";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 2000, -10000, 10000, 700, 0, 7000);
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

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "ImplantVelocityimplantE_AllDSSD";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 100, 0.6, 0.65, 280, 0, 7000);
			implantVelocityimplantEAll[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "ImplantZ";
			implantBetaHis = new TH1D(hisName.c_str(), "", 6, 0, 6);
			implantZ[i].push_back(implantBetaHis);
//////////////////////////////////////***gamma_EdT***//////////////////////////////////////////////////

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "beta_gamma_EdT_s_corr";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 1e2, 0, 10, 1000, 0, 7000);
			beta_gamma_EdT_s_corr[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "beta_gamma_EdT_ms_corr";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 100, 0, 1000, 1000, 0, 7000);
			beta_gamma_EdT_ms_corr[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "beta_gamma_EdT_us_corr";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 1000, 0, 10000, 1000, 0, 7000);
			beta_gamma_EdT_us_corr[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_beta_gamma_EdT_s_corr";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 1e2, 0, 10, 1000, 0, 7000);
			summed_beta_gamma_EdT_s_corr[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_beta_gamma_EdT_ms_corr";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 100, 0, 1000, 1000, 0, 7000);
			summed_beta_gamma_EdT_ms_corr[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_beta_gamma_EdT_us_corr";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 1000, 0, 10000, 1000, 0, 7000);
			summed_beta_gamma_EdT_us_corr[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "bp_gamma_EdT_s_corr";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 1e2, 0, 10, 1000, 0, 7000);
			bp_gamma_EdT_s_corr[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "bp_gamma_EdT_ms_corr";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 100, 0, 1000, 1000, 0, 7000);
			bp_gamma_EdT_ms_corr[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "bp_gamma_EdT_us_corr";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 1000, 0, 10000, 1000, 0, 7000);
			bp_gamma_EdT_us_corr[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_bp_gamma_EdT_s_corr";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 1e2, 0, 10, 1000, 0, 7000);
			summed_bp_gamma_EdT_s_corr[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_bp_gamma_EdT_ms_corr";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 100, 0, 1000, 1000, 0, 7000);
			summed_bp_gamma_EdT_ms_corr[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_bp_gamma_EdT_us_corr";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 1000, 0, 10000, 1000, 0, 7000);
			summed_bp_gamma_EdT_us_corr[i].push_back(implantBetaHis2D);

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
			/*singles
			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "beta_gamma_1";
			implantBetaHis = new TH1D(hisName.c_str(), "", 1000,0,7000);
			beta_gamma_1[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "beta_gamma_2";
			implantBetaHis = new TH1D(hisName.c_str(), "", 1000,0,7000);
			beta_gamma_2[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "beta_gamma_3";
			implantBetaHis = new TH1D(hisName.c_str(), "", 1000,0,7000);
			beta_gamma_3[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "beta_gamma_4";
			implantBetaHis = new TH1D(hisName.c_str(), "", 1000,0,7000);
			beta_gamma_4[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "beta_gamma_corr";
			implantBetaHis = new TH1D(hisName.c_str(), "", 1000,0,7000);
			beta_gamma_corr[i].push_back(implantBetaHis);

			//summed
			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_beta_gamma_1";
			implantBetaHis = new TH1D(hisName.c_str(), "", 1000,0,7000);
			summed_beta_gamma_1[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_beta_gamma_2";
			implantBetaHis = new TH1D(hisName.c_str(), "", 1000,0,7000);
			summed_beta_gamma_2[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_beta_gamma_3";
			implantBetaHis = new TH1D(hisName.c_str(), "", 1000,0,7000);
			summed_beta_gamma_3[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_beta_gamma_4";
			implantBetaHis = new TH1D(hisName.c_str(), "", 1000,0,7000);
			summed_beta_gamma_4[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_beta_gamma_corr";
			implantBetaHis = new TH1D(hisName.c_str(), "", 1000,0,7000);
			summed_beta_gamma_corr[i].push_back(implantBetaHis);

///////////////////////////////////*****beta-delayed proton gammas******/////////////////////////////////////////////////////
			/*//singles
			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "bp_gamma_1";
			implantBetaHis = new TH1D(hisName.c_str(), "", 1000,0,7000);
			bp_gamma_1[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "bp_gamma_2";
			implantBetaHis = new TH1D(hisName.c_str(), "", 1000,0,7000);
			bp_gamma_2[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "bp_gamma_3";
			implantBetaHis = new TH1D(hisName.c_str(), "", 1000,0,7000);
			bp_gamma_3[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "bp_gamma_4";
			implantBetaHis = new TH1D(hisName.c_str(), "", 1000,0,7000);
			bp_gamma_4[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "bp_gamma_corr";
			implantBetaHis = new TH1D(hisName.c_str(), "", 1000,0,7000);
			bp_gamma_corr[i].push_back(implantBetaHis);

			//summed
			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_bp_gamma_1";
			implantBetaHis = new TH1D(hisName.c_str(), "", 1000,0,7000);
			summed_bp_gamma_1[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_bp_gamma_2";
			implantBetaHis = new TH1D(hisName.c_str(), "", 1000,0,7000);
			summed_bp_gamma_2[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_bp_gamma_3";
			implantBetaHis = new TH1D(hisName.c_str(), "", 1000,0,7000);
			summed_bp_gamma_3[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_bp_gamma_4";
			implantBetaHis = new TH1D(hisName.c_str(), "", 1000,0,7000);
			summed_bp_gamma_4[i].push_back(implantBetaHis);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_bp_gamma_corr";
			implantBetaHis = new TH1D(hisName.c_str(), "", 1000,0,7000);
			summed_bp_gamma_corr[i].push_back(implantBetaHis);
			
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