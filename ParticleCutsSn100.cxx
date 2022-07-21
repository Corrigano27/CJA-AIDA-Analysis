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

//std::vector<TH2D *> implantVelocityimplantZ[numElements];
//std::vector<TH2D *> implantVelocityimplantE[numElements][6];

std::vector<TH1D *> implantBetaAll[numElements];
std::vector<TH1D *> implantBeta1pAll[numElements];
std::vector<TH1D *> implant1pAll[numElements];
std::vector<TH1D *> decayEnergyAll[numElements];

std::vector<TH1D *> delayed1pEnergyAll_AllDSSD_Ex[numElements][2];
std::vector<TH1D *> delayed1pEnergyAll_AllDSSD_ExSumCorr[numElements][2];
std::vector<TH1D *> delayed1pEnergyAll_AllDSSD_ExMax[numElements][2];
std::vector<TH1D *> delayed1pEnergyAll_AllDSSD_ExMaxSumCorr[numElements][2];
std::vector<TH2D *> EdTAll_NoMultiGate[numElements];
std::vector<TH2D *> EdTAll_NoMultiGate_[numElements];
std::vector<TH2D *> EdTAll_ms[numElements];
std::vector<TH2D *> EdTAll_us[numElements];
//std::vector<TH2D *> implantVelocityimplantEAll[numElements];

//correcting beta-gamma singles
std::vector<TH2D *> beta_gamma_EdT_s[numElements][2];
std::vector<TH2D *> beta_gamma_EdT_ms[numElements][2];
std::vector<TH2D *> beta_gamma_EdT_us[numElements][2];

//correcting beta-gamma summed
std::vector<TH2D *> summed_beta_gamma_EdT_s[numElements][2];
std::vector<TH2D *> summed_beta_gamma_EdT_ms[numElements][2];
std::vector<TH2D *> summed_beta_gamma_EdT_us[numElements][2];


//correcting bp-gamma singles
std::vector<TH2D *> bp_gamma_EdT_s[numElements][2];
std::vector<TH2D *> bp_gamma_EdT_ms[numElements][2];
std::vector<TH2D *> bp_gamma_EdT_us[numElements][2];


//correcting bp-gamma summed
std::vector<TH2D *> summed_bp_gamma_EdT_s[numElements][2];
std::vector<TH2D *> summed_bp_gamma_EdT_ms[numElements][2];
std::vector<TH2D *> summed_bp_gamma_EdT_us[numElements][2];



//Sn101 specific
TH2D *Tin101_bp_gamma_peak[2];
TH2D *Tin101_bp_gamma_rest[2];

TH2D *Tin101_summed_bp_gamma_peak[2];
TH2D *Tin101_summed_bp_gamma_rest[2];

//95Ag specific

TH2D *Ag95_EdT_160keVgammaGated;
TH2D *Ag95_EdT_800_1000keVgammaGated;
TH2D *Ag95_EdT_160_800_1000keVgammaGated;
TH2D *Ag95_EdT_440keVgammaGated;
TH2D *Ag95_EdT_511keVgammaGated;
TH2D *Ag95_EdT_allpeaks_gammaGated;

TH2D *Ag95_EdT_randomcheck;

TH2D *Ag95_EdT_2104keVsummed_gammaGated;
TH2D *Ag95_EdT_2104keVsummed_gammaGated_back;
TH2D *Ag95_EdT_2104keVsummed_gammaGated_11;
TH2D *Ag95_Implant_EdT_2104keVsummed_gammaGated;

TH2D *Ag95_EdT_77keVsummed_gammaGated;

TH2D *Ag95_EDiff_dT_2104keVsummed_gammaGated;

TH2D *Ag95_single_vs_summed;
TH2D *Ag95_gamma_gamma;
TH2D *Ag95_single_vs_summed_shorter;
TH2D *Ag95_gamma_gamma_shorter;
//Ag94 specific

//isomeric states
TH2D *Cd96_gs_EdT;
TH2D *Cd96_m_EdT;

TH2D *Cd97_gs_EdT;
TH2D *Cd97_m_EdT;
TH2D *Cd97_n_EdT;

//end of gammas

//std::vector<TH2D *> EnergyXChannelAll[numElements];
//std::vector<TH2D *> EnergyYChannelAll[numElements];
//std::vector<TH2D *> ExEyAll[numElements];
//std::vector<TH1D *> ExEyDiffAll[numElements];

std::vector<TH2D *> implantVelocityAOQ[numElements][6];
std::vector<TH2D *> implantVelocityAOQ_AllDSSD[numElements];

//std::vector<TH2D *> implantEnergyAOQ[numElements][6];
//std::vector<TH2D *> implantEnergyAOQ_AllDSSD[numElements];

//implants
std::vector<TH1D *> implantZ[numElements];
//std::vector<TH1D *> implantE[numElements][6];


//template histograms

TH1D *implantBetaHis;
TH2D *implantBetaHis2D;

//pid
//TH2D *PID;

//TH2D *PID_noise;

TH2D *implant_multi;

TH2D *implantE_Ex;

TH2D *PID_implant;

//TH2D *Indium97_gammaveto_EdT;

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

	implant_multi = new TH2D("implant_clustersize", "", 10, 0, 10, 10, 0, 10);

	implantE_Ex = new TH2D("implantE_Ex","",1000,0,10000,1000,0,10000);

	Ag95_EdT_160keVgammaGated = new TH2D("Ag95_EdT_160keV","",200,-10000,10000,80,200,1000);
	Ag95_EdT_800_1000keVgammaGated = new TH2D("Ag95_EdT_800_1000keV","",200,-10000,10000,80,200,1000);
	Ag95_EdT_160_800_1000keVgammaGated = new TH2D("Ag95_EdT_160_800_1000keV","",200,-10000,10000,80,200,1000);
	Ag95_EdT_440keVgammaGated = new TH2D("Ag95_EdT_440keV","",200,-10000,10000,80,200,1000);
	Ag95_EdT_511keVgammaGated = new TH2D("Ag95_EdT_511keV","",200,-10000,10000,80,200,1000);
	Ag95_EdT_allpeaks_gammaGated = new TH2D("Ag95_EdT_allpeaks","",200,-10000,10000,80,200,1000);

	Ag95_EdT_randomcheck = new TH2D("Ag95_EdT_randomcheck","",200,-10000,10000,80,200,1000);

	Ag95_EdT_2104keVsummed_gammaGated = new TH2D("Ag95_2104keVsummed_300eec","",200,-10000,10000,80,200,1000);
	Ag95_EdT_2104keVsummed_gammaGated_back = new TH2D("Ag95_2104keVsummed_100eec","",200,-10000,10000,80,200,1000);
	Ag95_EdT_2104keVsummed_gammaGated_11 = new TH2D("Ag95_2104keVsummed_front_11","",200,-10000,10000,80,200,1000);
	Ag95_Implant_EdT_2104keVsummed_gammaGated = new TH2D("Ag95_2104keVsummed_implants","",200,-10000,10000,80,200,7000);

	Ag95_EdT_77keVsummed_gammaGated = new TH2D("Ag95_77keVsummed","",600,-30000,30000,80,200,1000);

	Ag95_EDiff_dT_2104keVsummed_gammaGated = new TH2D("Ag95_2104keVsummed_200eec","",200,-10000,10000,80,200,1000);

	//Ag96_GammaT_betaT_all3Peaks = new TH2D("Ag96_GammaT_betaT_all3Peaks","",4e2,-5000,5000,2e2,-10,30);

	//Ag96_E_correlatedGamma = new TH1D("Ag96_E_correlatedGamma","",200,0,1000);

	//Ag96_E_randomGamma = new TH1D("Ag96_E_randomGamma","",200,0,1000);

	Ag95_single_vs_summed = new TH2D("Ag95_DTAS_single_vs_summed","",400,0,4000,400,0,4000);
	Ag95_single_vs_summed_shorter = new TH2D("Ag95_DTAS_single_vs_summed_shorter","",400,0,4000,400,0,4000);
	//Ag96_single_vs_summed = new TH2D("Ag96_DTAS_single_vs_summed","",400,0,4000,400,0,4000);


	Ag95_gamma_gamma = new TH2D("Ag95_gamma_gamma","",400,0,4000,400,0,4000);
	Ag95_gamma_gamma_shorter = new TH2D("Ag95_gamma_gamma_shorter","",400,0,4000,400,0,4000);
	//Ag96_gamma_gamma = new TH2D("Ag96_gamma_gamma","",400,0,4000,400,0,4000);

	Cd96_gs_EdT = new TH2D("Cd96_groundstate_EdT","",2e2,-10,10,700,0,7000);
	Cd96_m_EdT = new TH2D("Cd96_m_EdT","",2e2,-10,10,700,0,7000);
	
	Cd97_gs_EdT = new TH2D("Cd97_groundstate_EdT","",2e2,-10,10,700,0,7000);
	Cd97_m_EdT = new TH2D("Cd97_m_EdT","",2e2,-10,10,700,0,7000);
	Cd97_n_EdT = new TH2D("Cd97_n_EdT","",2e2,-10,10,700,0,7000);

	for (int p = 0; p<2; p++){
		hisName = "Tin101_bp_gamma_peak_" + std::to_string(p);
		Tin101_bp_gamma_peak[p] = new TH2D(hisName.c_str(), "",2e2,-10,10, 700, 0, 7000);

		hisName = "Tin101_bp_gamma_rest_" + std::to_string(p);
		Tin101_bp_gamma_rest[p] = new TH2D(hisName.c_str(), "",2e2,-10,10, 700, 0, 7000);

		hisName = "Tin101_summed_bp_gamma_peak_" + std::to_string(p);
		Tin101_summed_bp_gamma_peak[p] = new TH2D(hisName.c_str(), "",2e2,-10,10, 700, 0, 7000);

		hisName = "Tin101_summed_bp_gamma_rest_" + std::to_string(p);
		Tin101_summed_bp_gamma_rest[p] = new TH2D(hisName.c_str(), "", 2e2,-10,10, 700, 0, 7000);
	}

	for (int i = 0; i < numElements; i++)
	{
		
		for (int j = 0; j <= isotopeEnd[i] - isotopeStart[i]; j++)
		{
		//EdT_gamma_histograms
			//singles
			for (int g = 0; g < 2; g++){
			
				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "beta_gamma_EdT_s_" + std::to_string(g);
				implantBetaHis2D = new TH2D(hisName.c_str(), "", 2e2, -10, 10, 700, 0, 7000);
				beta_gamma_EdT_s[i][g].push_back(implantBetaHis2D);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "beta_gamma_EdT_ms_" + std::to_string(g);
				implantBetaHis2D = new TH2D(hisName.c_str(), "", 4000, -1000, 1000, 700, 0, 7000);
				beta_gamma_EdT_ms[i][g].push_back(implantBetaHis2D);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "beta_gamma_EdT_us_" + std::to_string(g);
				implantBetaHis2D = new TH2D(hisName.c_str(), "", 4000, -20000, 20000, 700, 0, 7000);
				beta_gamma_EdT_us[i][g].push_back(implantBetaHis2D);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "bp_gamma_EdT_s_" + std::to_string(g);
				implantBetaHis2D = new TH2D(hisName.c_str(), "", 2e2, -10, 10, 700, 0, 7000);
				bp_gamma_EdT_s[i][g].push_back(implantBetaHis2D);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "bp_gamma_EdT_ms_" + std::to_string(g);
				implantBetaHis2D = new TH2D(hisName.c_str(), "", 400, -1000, 1000, 700, 0, 7000);
				bp_gamma_EdT_ms[i][g].push_back(implantBetaHis2D);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "bp_gamma_EdT_us_" + std::to_string(g);
				implantBetaHis2D = new TH2D(hisName.c_str(), "", 4000, -20000, 20000, 700, 0, 7000);
				bp_gamma_EdT_us[i][g].push_back(implantBetaHis2D);
					
				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_beta_gamma_EdT_s_" + std::to_string(g);
				implantBetaHis2D = new TH2D(hisName.c_str(), "", 2e2, -10, 10, 500, 0, 10000);
				summed_beta_gamma_EdT_s[i][g].push_back(implantBetaHis2D);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_beta_gamma_EdT_ms_" + std::to_string(g);
				implantBetaHis2D = new TH2D(hisName.c_str(), "", 4000, -1000, 1000, 500, 0, 10000);
				summed_beta_gamma_EdT_ms[i][g].push_back(implantBetaHis2D);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_beta_gamma_EdT_us_" + std::to_string(g);
				implantBetaHis2D = new TH2D(hisName.c_str(), "", 4000, -20000, 20000, 500, 0, 10000);
				summed_beta_gamma_EdT_us[i][g].push_back(implantBetaHis2D);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_bp_gamma_EdT_s_" + std::to_string(g);
				implantBetaHis2D = new TH2D(hisName.c_str(), "", 2e2, -10, 10, 500, 0, 10000);
				summed_bp_gamma_EdT_s[i][g].push_back(implantBetaHis2D);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_bp_gamma_EdT_ms_" + std::to_string(g);
				implantBetaHis2D = new TH2D(hisName.c_str(), "", 400, -1000, 1000, 500, 0, 10000);
				summed_bp_gamma_EdT_ms[i][g].push_back(implantBetaHis2D);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_bp_gamma_EdT_us_" + std::to_string(g);
				implantBetaHis2D = new TH2D(hisName.c_str(), "", 4000, -20000, 20000, 500, 0, 10000);
				summed_bp_gamma_EdT_us[i][g].push_back(implantBetaHis2D);

				//beta-p
				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "delayed1pEnergyTotal_AllDSSD_Ex" + std::to_string(g);
				implantBetaHis = new TH1D(hisName.c_str(), "", 250, 1, 6);
				delayed1pEnergyAll_AllDSSD_Ex[i][g].push_back(implantBetaHis);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "delayed1pEnergyTotal_AllDSSD_ExSumCorr" + std::to_string(g);
				implantBetaHis = new TH1D(hisName.c_str(), "", 250, 1, 6);
				delayed1pEnergyAll_AllDSSD_ExSumCorr[i][g].push_back(implantBetaHis);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "delayed1pEnergyTotal_AllDSSD_ExMax" + std::to_string(g);
				implantBetaHis = new TH1D(hisName.c_str(), "", 250, 1, 6);
				delayed1pEnergyAll_AllDSSD_ExMax[i][g].push_back(implantBetaHis);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "delayed1pEnergyTotal_AllDSSD_ExMaxSumCorr" + std::to_string(g);
				implantBetaHis = new TH1D(hisName.c_str(), "", 250, 1, 6);
				delayed1pEnergyAll_AllDSSD_ExMaxSumCorr[i][g].push_back(implantBetaHis);
				

			}

			//use dssd arrays here
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

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "EdT_AllDSSD_NoMultiGate_s";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 2e2, -10, 10, 700, 0, 7000);
			EdTAll_NoMultiGate[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "EdT_AllDSSD_NoMultigate_ms";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 4000, -1000, 1000, 700, 0, 7000);
			EdTAll_ms[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "EdT_AllDSSD_NoMultigate_us";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 4000, -20000, 20000, 700, 0, 7000);
			EdTAll_us[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "ImplantZ";
			implantBetaHis = new TH1D(hisName.c_str(), "", 6, 0, 6);
			implantZ[i].push_back(implantBetaHis);

//////////////////////////////////////***gamma_EdT***//////////////////////////////////////////////////
			
//////////////////////////////////////***gamma-gated histos****////////////////////////////////////////
			
			/////////////////////////////////////////********implant stuff**********//////////////////////////////////////////

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "ImplantVelocityAOQ_All";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 100, 1.95, 2.15, 100, 0.6, 0.7);
			implantVelocityAOQ_AllDSSD[i].push_back(implantBetaHis2D);

		}
	}
}

#endif