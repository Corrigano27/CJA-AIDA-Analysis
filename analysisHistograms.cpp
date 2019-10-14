#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCut.h"
#include "TCutG.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

#include <string>
#include <iostream>
#include <map>
#include <fstream>

#include "/Disk/ds-sopa-personal/s1333561/PhD/MergerSoftware/data2Tree.cxx"
//#include "/home/corrigan/AidaSoftware/MergerSoftware/data2Tree.cxx"
#include "ParticleCutsSn100.cxx"

int analysisHistograms(std::string iName, std::string cutFile){

	std::string rootFile;
	std::string oName;
	std::ifstream rootFiles;

	//input text file with list of merged root files to be chained together.

	rootFiles.open( iName );
	if (!rootFiles.is_open()){
		std::cerr << " Problem opening input file, ending program." << std::endl;
		return -1; 
	}
	else {
		std::cout << "File ''" << iName << "'' is open" << std::endl;
	}

	TChain chain("mergedData");

	while ( std::getline ( rootFiles, rootFile )){

		chain.Add( rootFile.c_str() );
		std::cout << "Added " << rootFile.c_str() << " to the chain." << std::endl;

	}

 	size_t lastindex = iName.find_last_of(".");
	oName = iName.substr(0, lastindex);
	oName+="_AnalysisHistograms.root";
  	//Open the tree and create the branch to write to
  	TFile * ofile = TFile::Open( oName.c_str(), "recreate");

	std::cout << "Input and output files open" << std::endl;

	SetParticles();
	ReadParticleCuts(cutFile);
	SetImplantDSSD();
	DefineHistograms();

	std::cout << "ParticleCuts.cxx methods implemented" << std::endl;

	TTreeReader aReader( &chain );
	TTreeReaderValue <brData2TTree>    bigrips  (aReader, "bigrips.");
	TTreeReaderValue <impData2TTree>   implant  (aReader, "implantation.");
	TTreeReaderValue <betaData2TTree>  beta     (aReader, "beta.");
	TTreeReaderValue <gammaData2TTree> gamma    (aReader, "gamma.");
	TTreeReaderValue <ancData2TTree>   ancillary(aReader, "ancillary.");

	std::cout << "Tree reader set up" << std::endl;

	
	//Files read, histograms filled
	while (aReader.Next()){

		if ((*beta).T){
			if ((*beta).Ey >= 0.0 && (*beta).Ex>=0){
				for ( auto imp:(*beta).vectorOfImp ){
					for (int i = 0; i < numElements; i++){
						for (int j = 0; j <= isotopeEnd[i]-isotopeStart[i]; j++){
							if(particleCuts[i][j]->IsInside((imp).AOQ,(imp).ZET)){
								if ((*beta).z >= isotopeDSSDStart[i].at(j) && (*beta).z <= isotopeDSSDEnd[i].at(j)){
									edT[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).E);
									decayEnergy[i].at(j)->Fill((*beta).E);
									EDiff[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).Ex-(*beta).Ey);
									if ((*beta).nx < 3 && (*beta).ny < 3){
										if ((*beta).E>1500){
											delayed1pEnergy[i].at(j)->Fill((*beta).E);
											implantBeta1p[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9);
										}
										
									}
									else if ((*beta).E<1500){
										implantBeta[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9);
									}

								}
							}

						}
					}
				}
			}

		}

		if ((*bigrips).T>0){
			PID->Fill((*bigrips).aoq, (*bigrips).zet);
		}

		if ((*implant).T){
			for (int i =0; i<numElements; i++){
				for (int j = 0; j <= isotopeEnd[i]-isotopeStart[i]; j++){	
					if (particleCuts[i][j]->IsInside((*implant).aoq, (*implant).zet)){	//seg fault this line
						implantZ[i].at(j)->Fill((*implant).z);
					}
				}

			}
		}


	}//end of loop through chain

	ofile->cd();

	std::cout << "Writing to file" << std::endl;

	PID->Write();

	for(int i = 0; i < numElements; i++){
		for(unsigned int k = 0; k < decayEnergy[i].size(); k++){
				decayEnergy[i].at(k)->Write();
		}
		for(unsigned int k = 0; k < delayed1pEnergy[i].size(); k++){
				delayed1pEnergy[i].at(k)->Write();
		}
		for(unsigned int k = 0; k < implantBeta[i].size(); k++){
				implantBeta[i].at(k)->Write();
		}
		for(unsigned int k = 0; k < implantBeta1p[i].size(); k++){
				implantBeta1p[i].at(k)->Write();
		}
		for(unsigned int k = 0; k < implantZ[i].size(); k++){
				implantZ[i].at(k)->Write();
		}
		for(unsigned int k = 0; k < edT[i].size(); k++){
				edT[i].at(k)->Write();
		}
		for(unsigned int k = 0; k < EDiff[i].size(); k++){
				EDiff[i].at(k)->Write();
		}

	}


	ofile->Write();
	ofile->Close();

}