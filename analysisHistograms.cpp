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

	//boolean variable for beta vetoes

	bool betaVeto;

	uint8_t multix = 0;
	uint8_t multiy = 0;
	
	//Files read, histograms filled
	while (aReader.Next()){

		if ((*beta).T){
			if ((*beta).Ey >= 0.0 && (*beta).Ex>=0.0){
				multix = (*beta).TFast & 0xFF;
				multiy = ((*beta).TFast >> 8) & 0xFF;
				for ( auto imp:(*beta).vectorOfImp ){ //if non-element gated histos needed, do here
				betaVeto = false;
					for(auto anc:(*beta).vectorOfAnc){
					//AIDA plastic veto
						if((*beta).T - anc.TIME < 20e3 && (anc.ID == 34)){
							if((*beta).T - anc.TIME > 10e3 && (anc.ID == 34)){

								betaVeto = true;
							}
						}
											
					}
					if (betaVeto == false){
						if (multix < 3 && multiy < 3){
							edT_All->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).E);
							dT_All->Fill(((*beta).T-(imp).TIME)/1.0e3);
							if ((*beta).E>1500 && (*beta).E<10000 && (*beta).T-(imp).TIME)/1e3 >0 && ((*beta).T-(imp).TIME)/1e3 <200){
								PIDfast->Fill((imp).AOQ, (imp).ZET);
							}
						}
					}
					if (multix < 3 && multiy < 3){
						edT_All_beforeVeto->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).E);
					}

					for (int i = 0; i < numElements; i++){
						for (int j = 0; j <= isotopeEnd[i]-isotopeStart[i]; j++){
							if(particleCuts[i][j]->IsInside((imp).AOQ,(imp).ZET)){
								if ((*beta).z >= isotopeDSSDStart[i].at(j) && (*beta).z <= isotopeDSSDEnd[i].at(j)){
									//start applying vetoes here
									betaVeto = false;
										//initialise veto as false, then set true when conditions are met. Fill histograms when false
										
									for(auto anc:(*beta).vectorOfAnc){
										//AIDA Plastic veto (beta)
										if((*beta).T - anc.TIME < 20e3 && (anc.ID == 34)){
											if((*beta).T - anc.TIME > 10e3 && (anc.ID == 34)){
												betaVeto = true;
											}
										}

										//F11 veto (beta)
										//if((*beta).T - anc.TIME < 40e3 && (anc.ID == 32 || anc.ID == 33)){
											//if((*beta).T - anc.TIME > 0 && (anc.ID == 32 || anc.ID == 33)){
												//betaVeto = true;
											//}
										//}
											
									}
									
									if (betaVeto == false){
										decayEnergy[i].at(j)->Fill((*beta).E);
										//EDiff[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).Ex-(*beta).Ey);
										//EDiffLong[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9, (*beta).Ex-(*beta).Ey);
										if (multix <3 && multiy <3){
											edT[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).E);
											edTLong[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).E);
											edTMid[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).E);
											edTMid2[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).E);

											//time projections at different energies
											if ((*beta).E<520){
												dTMidUnder520[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e3);
												dTMid2Under520[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e3);
											}
											
											if ((*beta).E>520){
												dTMidOver520[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e3);
												dTMid2Over520[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e3);
											}

											if ((*beta).E>1500){
												delayed1pEnergy[i].at(j)->Fill((*beta).E);
												implantBeta1p[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9);
											}//end of lower beta-p energy cut
										
										}
										else if (multix < 5 && multiy < 5 && (*beta).E<1500){
											implantBeta[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9);
										}//end of upper beta energy cut
									}//end of beta veto application
									

								}//end of stopping layer if statement
							}//end of particle cut if statement

						}//end of isotope for loop
					}//end of elements for loop
				}//end of loop over correlated events
			} //end of if beta events with positive energy

		}//end of loop through beta events

		if ((*bigrips).T>0){
			PID->Fill((*bigrips).aoq, (*bigrips).zet);
		}

		if ((*implant).T){
			for (int i =0; i<numElements; i++){
				for (int j = 0; j <= isotopeEnd[i]-isotopeStart[i]; j++){	
					if (particleCuts[i][j]->IsInside((*implant).aoq, (*implant).zet)){
						implantZ[i].at(j)->Fill((*implant).z);
					}
				}

			}
		}


	}//end of loop through chain

	ofile->cd();

	std::cout << "Writing to file" << std::endl;

	PID->Write();

	PIDfast->Write();

	edT_All->Write();

	dT_All->Write();

	edT_All_beforeVeto->Write();

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
		//for(unsigned int k = 0; k < EDiff[i].size(); k++){
		//		EDiff[i].at(k)->Write();
		//}
		for(unsigned int k = 0; k < edTLong[i].size(); k++){
				edTLong[i].at(k)->Write();
		}
		//for(unsigned int k = 0; k < EDiffLong[i].size(); k++){
		//		EDiffLong[i].at(k)->Write();
		//}
		for(unsigned int k = 0; k < edTMid[i].size(); k++){
				edTMid[i].at(k)->Write();
		}
		for(unsigned int k = 0; k < edTMid2[i].size(); k++){
				edTMid2[i].at(k)->Write();
		}
		for(unsigned int k = 0; k < dTMidUnder520[i].size(); k++){
				dTMidUnder520[i].at(k)->Write();
		}
		for(unsigned int k = 0; k < dTMidOver520[i].size(); k++){
				dTMidOver520[i].at(k)->Write();
		}
		for(unsigned int k = 0; k < dTMid2Under520[i].size(); k++){
				dTMid2Under520[i].at(k)->Write();
		}
		for(unsigned int k = 0; k < dTMid2Over520[i].size(); k++){
				dTMid2Over520[i].at(k)->Write();
		}
		




	}


	ofile->Write();
	ofile->Close();

}