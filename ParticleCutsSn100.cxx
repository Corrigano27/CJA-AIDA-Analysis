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
std::vector<TH2D *> EdTAll_ms[4];
std::vector<TH2D *> EdTAll_us[4];
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
std::vector<TH2D *> beta_gamma_EdT_s[4][4];
std::vector<TH2D *> beta_gamma_EdT_ms[4][4];
std::vector<TH2D *> beta_gamma_EdT_us[4][4];

std::vector<TH2D *> beta_gamma_EdT_s_corr[4];
std::vector<TH2D *> beta_gamma_EdT_ms_corr[4];
std::vector<TH2D *> beta_gamma_EdT_us_corr[4];

//correcting beta-gamma summed
std::vector<TH2D *> summed_beta_gamma_EdT_s[4][2];
std::vector<TH2D *> summed_beta_gamma_EdT_ms[4][2];
std::vector<TH2D *> summed_beta_gamma_EdT_us[4][2];

std::vector<TH2D *> summed_beta_gamma_EdT_s_corr[4];
std::vector<TH2D *> summed_beta_gamma_EdT_ms_corr[4];
std::vector<TH2D *> summed_beta_gamma_EdT_us_corr[4];

std::vector<TH2D *> summed_beta_gamma_E_beta_E[4][2];


//correcting bp-gamma singles
std::vector<TH2D *> bp_gamma_EdT_s[4][4];
std::vector<TH2D *> bp_gamma_EdT_ms[4][4];
std::vector<TH2D *> bp_gamma_EdT_us[4][4];

std::vector<TH2D *> bp_gamma_EdT_s_corr[4];
std::vector<TH2D *> bp_gamma_EdT_ms_corr[4];
std::vector<TH2D *> bp_gamma_EdT_us_corr[4];


//correcting bp-gamma summed
std::vector<TH2D *> summed_bp_gamma_EdT_s[4][2];
std::vector<TH2D *> summed_bp_gamma_EdT_ms[4][2];
std::vector<TH2D *> summed_bp_gamma_EdT_us[4][2];

std::vector<TH2D *> summed_bp_gamma_EdT_s_corr[4];
std::vector<TH2D *> summed_bp_gamma_EdT_ms_corr[4];
std::vector<TH2D *> summed_bp_gamma_EdT_us_corr[4];

std::vector<TH2D *> summed_p_gamma_E_p_E[4][2];


//Sn101 specific
TH1D *Tin101_bp_gamma_peak[4];
TH1D *Tin101_bp_gamma_peak_corr;
TH1D *Tin101_bp_gamma_rest[4];
TH1D *Tin101_bp_gamma_rest_corr;
TH1D *Tin101_summed_bp_gamma_peak[2];
TH1D *Tin101_summed_bp_gamma_peak_corr;
TH1D *Tin101_summed_bp_gamma_rest[2];
TH1D *Tin101_summed_bp_gamma_rest_corr;

//96Ag specific
TH2D *Ag96_EdT_470keVgammaGated;
TH2D *Ag96_EdT_743keVgammaGated;
TH2D *Ag96_EdT_1249keVgammaGated;
TH2D *Ag96_EdT_all3Peaks_gammaGated;
TH2D *Ag96_GammaT_betaT_all3Peaks;

TH2D *Ag96_EdT_all3Peaks_Random_gammaGated;

TH1D *Ag96_E_randomGamma;
TH1D *Ag96_E_correlatedGamma;

TH2D *Ag96_single_vs_summed;
TH2D *Ag96_gamma_gamma;

TH1D *Ag96_sum_E1E2_diff_470_740;
TH1D *Ag96_sum_E1E2_diff_740_1249;

TH2D *Ag96_EdT_2461keVgammaGated;
TH2D *Ag96_EdT_summed_gammaGated;
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
TH1D *Ag94_1800_bp_DTASindy[4];
TH1D *Ag94_1800_bp_DTASsummed[2];

//end of gammas

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
//TH2D *PID;

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

	Indium97_gammaveto_EdT = new TH2D("Indium97_gammaveto_EdT_ms","", 2000, -80, 80,700,0,7000);

	//Sn101 & Ag95/94
	for (int g = 0; g < 4; g++){
			
		hisName = "Tin101_bp_gamma_peak_" + std::to_string(g);
		Tin101_bp_gamma_peak[g] = new TH1D(hisName.c_str(), "", 1000, 0, 7000);

		hisName = "Tin101_bp_gamma_rest_" + std::to_string(g);
		Tin101_bp_gamma_rest[g] = new TH1D(hisName.c_str(), "", 1000, 0, 7000);

		if (g < 2){
			hisName = "Tin101_summed_bp_gamma_peak_" + std::to_string(g);
			Tin101_summed_bp_gamma_peak[g] = new TH1D(hisName.c_str(), "", 1000, 0, 7000);

			hisName = "Tin101_summed_bp_gamma_rest_" + std::to_string(g);
			Tin101_summed_bp_gamma_rest[g] = new TH1D(hisName.c_str(), "", 1000, 0, 7000);

			hisName = "Ag94_1800_bp_DTASsummed_" + std::to_string(g);
			Ag94_1800_bp_DTASsummed[g] = new TH1D(hisName.c_str(), "", 1000, 0, 7000);
			
		}

		hisName = "Ag94_1800_bp_DTASindy_" + std::to_string(g);
		Ag94_1800_bp_DTASindy[g] = new TH1D(hisName.c_str(), "", 1000, 0, 7000);

	}

	Ag96_EdT_470keVgammaGated = new TH2D("Ag96_EdT_470keVgammaGated","",1000,-5000,5000,80,200,1000);
	Ag96_EdT_743keVgammaGated = new TH2D("Ag96_EdT_743keVgammaGated","",1000,-5000,5000,80,200,1000);
	Ag96_EdT_1249keVgammaGated = new TH2D("Ag96_EdT_1249keVgammaGated","",1000,-5000,5000,80,200,1000);
	Ag96_EdT_all3Peaks_gammaGated = new TH2D("Ag96_EdT_all3Peaks_gammaGated","",1000,-5000,5000,80,200,1000);

	Ag96_EdT_all3Peaks_Random_gammaGated = new TH2D("Ag96_EdT_Random_gammaGated","",1000,-5000,5000,80,200,1000);
	Ag96_EdT_2461keVgammaGated = new TH2D("Ag96_EdT_2461keVgammaGated","",1000,-5000,5000,80,200,1000);
	Ag96_EdT_summed_gammaGated = new TH2D("Ag96_EdT_summed_gammaGated","",1000,-5000,5000,80,200,1000);

	Ag95_EdT_160keVgammaGated = new TH2D("Ag95_EdT_160keV","",200,-10000,10000,80,200,1000);
	Ag95_EdT_800_1000keVgammaGated = new TH2D("Ag95_EdT_800_1000keV","",200,-10000,10000,80,200,1000);
	Ag95_EdT_160_800_1000keVgammaGated = new TH2D("Ag95_EdT_160_800_1000keV","",200,-10000,10000,80,200,1000);
	Ag95_EdT_440keVgammaGated = new TH2D("Ag95_EdT_440keV","",200,-10000,10000,80,200,1000);
	Ag95_EdT_511keVgammaGated = new TH2D("Ag95_EdT_511keV","",200,-10000,10000,80,200,1000);
	Ag95_EdT_allpeaks_gammaGated = new TH2D("Ag95_EdT_allpeaks","",200,-10000,10000,80,200,1000);

	Ag95_EdT_randomcheck = new TH2D("Ag95_EdT_randomcheck","",200,-10000,10000,80,200,1000);

	Ag95_EdT_2104keVsummed_gammaGated = new TH2D("Ag95_2104keVsummed_front","",200,-10000,10000,80,200,1000);
	Ag95_EdT_2104keVsummed_gammaGated_back = new TH2D("Ag95_2104keVsummed_back","",200,-10000,10000,80,200,1000);
	Ag95_EdT_2104keVsummed_gammaGated_11 = new TH2D("Ag95_2104keVsummed_front_11","",200,-10000,10000,80,200,1000);
	Ag95_Implant_EdT_2104keVsummed_gammaGated = new TH2D("Ag95_2104keVsummed_implants","",200,-10000,10000,80,200,7000);

	Ag95_EdT_77keVsummed_gammaGated = new TH2D("Ag95_77keVsummed","",600,-30000,30000,80,200,1000);

	Ag95_EDiff_dT_2104keVsummed_gammaGated = new TH2D("Ag95_EDiff_dT_2104keVsummed_gammaGated","",200,-10000,10000,100,-600,600);

	Ag96_GammaT_betaT_all3Peaks = new TH2D("Ag96_GammaT_betaT_all3Peaks","",4e2,-5000,5000,2e2,-10,30);

	Ag96_E_correlatedGamma = new TH1D("Ag96_E_correlatedGamma","",200,0,1000);

	Ag96_E_randomGamma = new TH1D("Ag96_E_randomGamma","",200,0,1000);

	Ag95_single_vs_summed = new TH2D("Ag95_DTAS_single_vs_summed","",400,0,4000,400,0,4000);
	Ag95_single_vs_summed_shorter = new TH2D("Ag95_DTAS_single_vs_summed_shorter","",400,0,4000,400,0,4000);
	Ag96_single_vs_summed = new TH2D("Ag96_DTAS_single_vs_summed","",400,0,4000,400,0,4000);

	Ag96_sum_E1E2_diff_470_740 = new TH1D("Ag96_sum_E1E2_diff_470_740","",500,0,5000);
	Ag96_sum_E1E2_diff_740_1249 = new TH1D("Ag96_sum_E1E2_diff_740_1249","",500,0,5000);

	Ag95_gamma_gamma = new TH2D("Ag95_gamma_gamma","",400,0,4000,400,0,4000);
	Ag95_gamma_gamma_shorter = new TH2D("Ag95_gamma_gamma_shorter","",400,0,4000,400,0,4000);
	Ag96_gamma_gamma = new TH2D("Ag96_gamma_gamma","",400,0,4000,400,0,4000);

	Tin101_bp_gamma_peak_corr = new TH1D("Tin101_bp_gamma_peak_corr","",1000, 0, 7000);
	Tin101_summed_bp_gamma_peak_corr = new TH1D("Tin101_summed_bp_gamma_peak_corr","",1000, 0, 7000);
	Tin101_bp_gamma_rest_corr = new TH1D("Tin101_bp_gamma_rest_corr","",1000, 0, 7000);
	Tin101_summed_bp_gamma_rest_corr = new TH1D("Tin101_summed_bp_gamma_rest_corr","",1000, 0, 7000);


	for (int i = 0; i < numElements; i++)
	{
		
		for (int j = 0; j <= isotopeEnd[i] - isotopeStart[i]; j++)
		{
		//EdT_gamma_histograms
			//singles
			for (int g = 0; g < 4; g++){
			
				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "beta_gamma_EdT_s_" + std::to_string(g);
				implantBetaHis2D = new TH2D(hisName.c_str(), "", 1e2, 0, 10, 1000, 0, 7000);
				beta_gamma_EdT_s[i][g].push_back(implantBetaHis2D);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "beta_gamma_EdT_ms_" + std::to_string(g);
				implantBetaHis2D = new TH2D(hisName.c_str(), "", 2000, 0, 1000, 1000, 0, 7000);
				beta_gamma_EdT_ms[i][g].push_back(implantBetaHis2D);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "beta_gamma_EdT_us_" + std::to_string(g);
				implantBetaHis2D = new TH2D(hisName.c_str(), "", 2000, 0, 20000, 1000, 0, 7000);
				beta_gamma_EdT_us[i][g].push_back(implantBetaHis2D);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "bp_gamma_EdT_s_" + std::to_string(g);
				implantBetaHis2D = new TH2D(hisName.c_str(), "", 1e2, 0, 10, 1000, 0, 7000);
				bp_gamma_EdT_s[i][g].push_back(implantBetaHis2D);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "bp_gamma_EdT_ms_" + std::to_string(g);
				implantBetaHis2D = new TH2D(hisName.c_str(), "", 200, 0, 1000, 1000, 0, 7000);
				bp_gamma_EdT_ms[i][g].push_back(implantBetaHis2D);

				hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "bp_gamma_EdT_us_" + std::to_string(g);
				implantBetaHis2D = new TH2D(hisName.c_str(), "", 2000, 0, 20000, 1000, 0, 7000);
				bp_gamma_EdT_us[i][g].push_back(implantBetaHis2D);

				if (g < 2){

					hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_beta_gamma_EdT_s_" + std::to_string(g);
					implantBetaHis2D = new TH2D(hisName.c_str(), "", 1e2, 0, 10, 1000, 0, 10000);
					summed_beta_gamma_EdT_s[i][g].push_back(implantBetaHis2D);

					hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_beta_gamma_EdT_ms_" + std::to_string(g);
					implantBetaHis2D = new TH2D(hisName.c_str(), "", 2000, 0, 1000, 1000, 0, 10000);
					summed_beta_gamma_EdT_ms[i][g].push_back(implantBetaHis2D);

					hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_beta_gamma_EdT_us_" + std::to_string(g);
					implantBetaHis2D = new TH2D(hisName.c_str(), "", 2000, 0, 20000, 1000, 0, 10000);
					summed_beta_gamma_EdT_us[i][g].push_back(implantBetaHis2D);

					hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_bp_gamma_EdT_s_" + std::to_string(g);
					implantBetaHis2D = new TH2D(hisName.c_str(), "", 1e2, 0, 10, 1000, 0, 10000);
					summed_bp_gamma_EdT_s[i][g].push_back(implantBetaHis2D);

					hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_bp_gamma_EdT_ms_" + std::to_string(g);
					implantBetaHis2D = new TH2D(hisName.c_str(), "", 200, 0, 1000, 1000, 0, 10000);
					summed_bp_gamma_EdT_ms[i][g].push_back(implantBetaHis2D);

					hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_bp_gamma_EdT_us_" + std::to_string(g);
					implantBetaHis2D = new TH2D(hisName.c_str(), "", 2000, 0, 20000, 1000, 0, 10000);
					summed_bp_gamma_EdT_us[i][g].push_back(implantBetaHis2D);

					hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_p_gamma_E_p_E_" + std::to_string(g);
					implantBetaHis2D = new TH2D(hisName.c_str(), "", 700, 0, 7000, 1000, 0, 10000);
					summed_p_gamma_E_p_E[i][g].push_back(implantBetaHis2D);

					hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_beta_gamma_E_beta_E_" + std::to_string(g);
					implantBetaHis2D = new TH2D(hisName.c_str(), "", 700, 0, 7000, 1000, 0, 10000);
					summed_beta_gamma_E_beta_E[i][g].push_back(implantBetaHis2D);


				}

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
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 2000, 0, 1000, 1000, 0, 7000);
			beta_gamma_EdT_ms_corr[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "beta_gamma_EdT_us_corr";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 2000, 0, 20000, 1000, 0, 7000);
			beta_gamma_EdT_us_corr[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_beta_gamma_EdT_s_corr";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 1e2, 0, 10, 1000, 0, 10000);
			summed_beta_gamma_EdT_s_corr[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_beta_gamma_EdT_ms_corr";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 2000, 0, 1000, 1000, 0, 10000);
			summed_beta_gamma_EdT_ms_corr[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_beta_gamma_EdT_us_corr";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 2000, 0, 20000, 1000, 0, 10000);
			summed_beta_gamma_EdT_us_corr[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "bp_gamma_EdT_s_corr";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 1e2, 0, 10, 1000, 0, 7000);
			bp_gamma_EdT_s_corr[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "bp_gamma_EdT_ms_corr";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 200, 0, 1000, 1000, 0, 7000);
			bp_gamma_EdT_ms_corr[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "bp_gamma_EdT_us_corr";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 2000, 0, 20000, 1000, 0, 7000);
			bp_gamma_EdT_us_corr[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_bp_gamma_EdT_s_corr";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 1e2, 0, 10, 1000, 0, 10000);
			summed_bp_gamma_EdT_s_corr[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_bp_gamma_EdT_ms_corr";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 200, 0, 1000, 1000, 0, 10000);
			summed_bp_gamma_EdT_ms_corr[i].push_back(implantBetaHis2D);

			hisName = elements[i] + std::to_string(isotopeStart[i] + j) + "summed_bp_gamma_EdT_us_corr";
			implantBetaHis2D = new TH2D(hisName.c_str(), "", 2000, 0, 20000, 1000, 0, 10000);
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