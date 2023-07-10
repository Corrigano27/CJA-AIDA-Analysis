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

//data2tree file select
#include "/Disk/ds-sopa-personal/s1333561/PhD/MergerSoftware/data2Tree.cxx"
//#include "/home/corrigan/DTAS_Merger/merger/MergerSoft/data2Tree.cxx"
#include "ParticleCutsSn100.cxx"
#include "analysisHistograms.hpp"

int analysisHistograms(std::string iName, std::string cutFile){

	std::string rootFile;
	std::string oName;
	std::ifstream rootFiles;

	//input text file with list of merged root files to be chained together.

	ROOT::EnableImplicitMT(4);

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
	oName+="_SnInCdAg_AnalysisHistograms.root";
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
			if ((*beta).E>=300.0){
				if (abs((*beta).Ex-(*beta).Ey)>=0){
					multix = (*beta).TFast & 0xFF;
					multiy = ((*beta).TFast >> 8) & 0xFF;
					for ( auto imp:(*beta).vectorOfImp ){ //if non-element gated histos needed, do here
						if ((*beta).z == (imp).Z){
							gammaVeto = true;
							In97gammaVeto = false;
							
							for (int i = 0; i < numElements; i++){
								for (int j = 0; j <= isotopeEnd[i]-isotopeStart[i]; j++){
									if (particleCuts[i][j]->IsInside((imp).AOQ,(imp).ZET)){
										if ((*beta).z >= isotopeDSSDStart[i].at(j) && (*beta).z <= isotopeDSSDEnd[i].at(j)){// use these statements for the dssd loop later on
											//start applying vetoes here
											betaVeto = false;
											//initialise veto as false, then set true when conditions are met. Fill histograms when false
												
											for (auto anc:(*beta).vectorOfAnc){
												//AIDA Plastic veto (beta)
												if ((*beta).T - anc.TIME < 20e3 && (anc.ID == 34)){
													if ((*beta).T - anc.TIME > 10e3 && (anc.ID == 34)){
														betaVeto = true;
													}
												}

												//F11 veto (beta)
												if (anc.EN > 20){
													if((*beta).T - anc.TIME < 40e3 && (anc.ID == 32 || anc.ID == 33)){
														if((*beta).T - anc.TIME > 0 && (anc.ID == 32 || anc.ID == 33)){
															betaVeto = true;
														}
													}

												}
													
											}
											
											if (betaVeto == false){
												double idx = ((imp).TFAST & 0xFF)/2.0;

                        						double idy = (((imp).TFAST >> 8) & 0xFF)/2.0;

												double dx = multix/2.0;

												double dy = multiy/2.0;
												//use below to have variable dssd - will need to introduce further dssd vectors
												//if ((*beta).z >= isotopeDSSDStart[i].at(j) && (*beta).z <= isotopeDSSDEnd[i].at(j)){
												int DSSD = ((*beta).z);
											
												if (((imp).Y + (idy) >= (((*beta).y)-((dy)+1.0))) && ((imp).Y - (idy) <= (((*beta).y)+((dy)+1.0)))){
													if(((imp).X + (idx)>= (((*beta).x)-((dx)+1.0))) && ((imp).X - (idx)<= (((*beta).x)+((dx)+1.0)))){
														EdTAll_NoMultiGate[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9, (*beta).Ex);
														EdTAll_us[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).Ex);
														EdTAll_ms[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e6, (*beta).Ex);
													}
													//if ((*beta).Ex<1100){
														//implantBeta[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9);
													//}
												}

												//end of upper beta energy cut
												//end of dssd if
												//end of dssd for
												decayEnergyAll[i].at(j)->Fill((*beta).E);
												if (multix ==0 && multiy == 0){ //beta-delayed protons
													if ((*beta).Ex>1100){ //beta-delayed protons
														//beta-p gamma loop
														ProtonGammaSumTemp = 0;
														ProtonGammaSumTempBg = 0;

														for ( auto gamma:(*beta).vectorOfGamma ){ //loop over gamma events
															IndyGammaE = DTAS_SingleCalib(gamma.EN);													
															if (((*beta).T-(gamma).TIME) > 10000){ //forward gammas
																if (((*beta).T-(gamma).TIME) < 20000){
																	if ((gamma).ID<16){		
																		//bp_gamma_1[i].at(j)->Fill((gamma.EN));
																		bp_gamma_EdT_s[i][0].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9, IndyGammaE);
																		bp_gamma_EdT_ms[i][0].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e6, IndyGammaE);
																		bp_gamma_EdT_us[i][0].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e3, IndyGammaE);
																		ProtonGammaSumTemp+=(IndyGammaE);

																		if (elements[i] == "Sn" && isotopeStart[i]+j == 101){
																			if (((*beta).T-(imp).TIME) > 0){
																				Tin101_singles_bpE_gammaE[0]->Fill(((*beta).E), IndyGammaE);
																			}
																			
																		}
																		if (elements[i] == "In" && isotopeStart[i]+j == 98){
																			if (((*beta).T-(imp).TIME)>0){
																				In98_bpE_gammaE->Fill((*beta).E, IndyGammaE);
																			}
																		}


																		
																	}
																		
																}
																
															}
															if (((*beta).T-(gamma).TIME) > 20000){ //random gammas
																if(((*beta).T-(gamma).TIME) < 30000){
																	if ((gamma).ID<16){		
																		//bp_gamma_2[i].at(j)->Fill((IndyGammaE));
																		bp_gamma_EdT_s[i][1].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9, (IndyGammaE));
																		bp_gamma_EdT_ms[i][1].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e6, (IndyGammaE));
																		bp_gamma_EdT_us[i][1].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e3, (IndyGammaE));
																		ProtonGammaSumTempBg+=(IndyGammaE);

																		if (elements[i] == "Sn" && isotopeStart[i]+j == 101){
																			if (((*beta).T-(imp).TIME) > 0){
																				Tin101_singles_bpE_gammaE[0]->Fill(((*beta).E), IndyGammaE);
																			}
																			
																		}
																	}
																}
															}	
															
															
														}// end gamma loop
														
														if (ProtonGammaSumTemp != 0){
															summed_bp_gamma_EdT_s[i][0].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9, ProtonGammaSumTemp);
															summed_bp_gamma_EdT_ms[i][0].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e6, ProtonGammaSumTemp);
															summed_bp_gamma_EdT_us[i][0].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e3, ProtonGammaSumTemp);
															if (elements[i] == "Sn" && isotopeStart[i]+j == 101){
																Tin101_summed_bpE_gammaE[0]->Fill((*beta).E, ProtonGammaSumTemp);
															}

														}

														if (ProtonGammaSumTempBg != 0){
															summed_bp_gamma_EdT_s[i][1].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9, ProtonGammaSumTempBg);
															summed_bp_gamma_EdT_ms[i][1].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e6, ProtonGammaSumTempBg);
															summed_bp_gamma_EdT_us[i][1].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e3, ProtonGammaSumTempBg);
															if (elements[i] == "Sn" && isotopeStart[i]+j == 101){
																Tin101_summed_bpE_gammaE[1]->Fill((*beta).E, ProtonGammaSumTempBg);
															}
														}

													}//end of lower beta-p energy cut
												}//end of beta-p multiplicity cut
	

												//ExEyAll[i].at(j)->Fill((*beta).Ex, (*beta).Ey);

												//DTAS 511/1022 coincidence check
												//gamma gated spectra

												//beta - DTAS correlations
												
												///////BETAS//////////

												if ((*beta).Ex<1100){
													
													implantBetaAll[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9);
													
													GammaSumTemp=0;
													GammaSumTempBg=0;
													//ProtonGammaSumTemp=0;
													//ProtonGammaSumTempBg=0;
													for ( auto gamma:(*beta).vectorOfGamma){ //loop over gamma events
														IndyGammaE = DTAS_SingleCalib(gamma.EN);
														if (((*beta).T-(gamma).TIME) > 10000){ //forward gammas
															if(((*beta).T-(gamma).TIME) < 20000){
																//summed_beta_gamma_1[i].at(j)->Fill((IndyGammaE));

																
																if ((gamma).ID<16){

																	beta_gamma_EdT_s[i][0].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9, (IndyGammaE));
																	beta_gamma_EdT_ms[i][0].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e6, (IndyGammaE));
																	beta_gamma_EdT_us[i][0].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e3, (IndyGammaE));
																	GammaSumTemp+=(IndyGammaE);
																}
																	
															}
																
														}	

														if (((*beta).T-(gamma).TIME) > 20000){ //random gammas
															if(((*beta).T-(gamma).TIME) < 30000){
																if ((gamma).ID<16){	

																	beta_gamma_EdT_s[i][1].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9, (IndyGammaE));
																	beta_gamma_EdT_ms[i][1].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e6, (IndyGammaE));
																	beta_gamma_EdT_us[i][1].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e3, (IndyGammaE));
																	GammaSumTempBg+=(IndyGammaE);

																}
																
															}
														}
												
			
													}//end of gamma loop
													//fill tallied histograms
													gammaSubtract = 0;
													if (GammaSumTemp != 0){
														summed_beta_gamma_EdT_s[i][0].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9, GammaSumTemp);
														summed_beta_gamma_EdT_ms[i][0].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e6, GammaSumTemp);
														summed_beta_gamma_EdT_us[i][0].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e3, GammaSumTemp);

														//summed_beta_gamma_E_beta_E[i][0].at(j)->Fill((*beta).E, GammaSumTemp);

														
													}

													if (GammaSumTempBg != 0){
														summed_beta_gamma_EdT_s[i][1].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9, GammaSumTempBg);
														summed_beta_gamma_EdT_ms[i][1].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e6, GammaSumTempBg);
														summed_beta_gamma_EdT_us[i][1].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e3, GammaSumTempBg);
													}
													

														
												}//beta energy and multi-cut
												



											
											}//end of beta veto application
											

										}//end of stopping layer if statement
									}//end of particle cut if statement

								}//end of isotope for loop
							}//end of elements for loop
						}//end of imp = decay dssd if
					}//end of loop over correlated events
				}//end of equal energy cut
			} //end of if beta events with positive energy

		}//end of loop through beta events
		double imp_multix = 0;
		double imp_multiy = 0;
		if ((*implant).T){
			for ( auto pid:(*implant).vectorOfPid ){
				PID_implant->Fill((*implant).aoq, (*implant).zet);
				imp_multix = (*implant).TFast & 0xFF;
				imp_multiy = ((*implant).TFast >> 8) & 0xFF;
				implant_multi->Fill(imp_multix, imp_multiy);
				implantE_Ex->Fill((*implant).Ex, (*implant).E);
				int iDSSD = (*implant).z;
				for (int i =0; i<numElements; i++){
					for (int j = 0; j <= isotopeEnd[i]-isotopeStart[i]; j++){	
						if (particleCuts[i][j]->IsInside((*implant).aoq, (*implant).zet)){
							if (iDSSD >= isotopeDSSDStart[i].at(j) && iDSSD <= isotopeDSSDEnd[i].at(j)){
								implantZ[i].at(j)->Fill((*implant).z);
								implantVelocityAOQ_AllDSSD[i].at(j)->Fill((*implant).aoq,(pid).VELOCITY);
							}
														
						}
					}
				}
			}
			
		}


	}//end of loop through chain

	ofile->cd();

	std::cout << "Writing to file" << std::endl;

	implant_multi->Write();

	implantE_Ex->Write();

	PID_implant->Write();

	In98_bpE_gammaE->Write();

	//Tin101_summed_bp_gamma_peak[0]->Add(Tin101_summed_bp_gamma_peak[1],-1);
	//Tin101_summed_bp_gamma_rest[0]->Add(Tin101_summed_bp_gamma_rest[1],-1);

	//Tin101_bp_gamma_peak[0]->Add(Tin101_bp_gamma_peak[1],-1);
	//Tin101_bp_gamma_rest[0]->Add(Tin101_bp_gamma_rest[1],-1);

	//Tin101_bp_gamma_peak[0]->Write();
	//Tin101_summed_bp_gamma_peak[0]->Write();
	//Tin101_bp_gamma_rest[0]->Write();
	//Tin101_summed_bp_gamma_rest[0]->Write();

	Tin101_singles_bpE_gammaE[0]->Add(Tin101_singles_bpE_gammaE[1],-1);
	Tin101_summed_bpE_gammaE[0]->Add(Tin101_summed_bpE_gammaE[1],-1);

	Tin101_singles_bpE_gammaE[0]->Write();
	Tin101_summed_bpE_gammaE[0]->Write();

	std::string isoDirName;

	for(int i = 0; i < numElements; i++){

		ElDir = ofile->mkdir(elements[i].c_str());
		ElementDir.push_back(ElDir);
		//element-directory

		for (int k = 0; k <= isotopeEnd[i]-isotopeStart[i]; k++){

			isoDirName = elements[i].c_str() + std::to_string(isotopeStart[i] + k);
			IsoDir = ElDir->mkdir(isoDirName.c_str());
			IsotopeDir.push_back(IsoDir);
			//isotope-directory

			
			//combined DSSD
			IsoDir->Append(decayEnergyAll[i].at(k));			
			//IsoDir->Append(delayed1pEnergyAll_AllDSSD[i].at(k));
			IsoDir->Append(implantBetaAll[i].at(k));
			IsoDir->Append(implantBeta1pAll[i].at(k));
			IsoDir->Append(EdTAll_ms[i].at(k));
			IsoDir->Append(EdTAll_us[i].at(k));
			IsoDir->Append(EdTAll_NoMultiGate[i].at(k));
			IsoDir->Append(implantVelocityAOQ_AllDSSD[i].at(k));
			IsoDir->Append(implantZ[i].at(k));


			

			//gamma spectra correction
			beta_gamma_EdT_s[i][0].at(k)->Add(beta_gamma_EdT_s[i][1].at(k),-1);
			beta_gamma_EdT_ms[i][0].at(k)->Add(beta_gamma_EdT_ms[i][1].at(k),-1);
			beta_gamma_EdT_us[i][0].at(k)->Add(beta_gamma_EdT_us[i][1].at(k),-1);

			bp_gamma_EdT_s[i][0].at(k)->Add(bp_gamma_EdT_s[i][1].at(k),-1);
			bp_gamma_EdT_ms[i][0].at(k)->Add(bp_gamma_EdT_ms[i][1].at(k),-1);
			bp_gamma_EdT_us[i][0].at(k)->Add(bp_gamma_EdT_us[i][1].at(k),-1);

			summed_beta_gamma_EdT_s[i][0].at(k)->Add(summed_beta_gamma_EdT_s[i][1].at(k),-1);
			summed_beta_gamma_EdT_ms[i][0].at(k)->Add(summed_beta_gamma_EdT_ms[i][1].at(k),-1);
			summed_beta_gamma_EdT_us[i][0].at(k)->Add(summed_beta_gamma_EdT_us[i][1].at(k),-1);

			summed_bp_gamma_EdT_s[i][0].at(k)->Add(summed_bp_gamma_EdT_s[i][1].at(k),-1);
			summed_bp_gamma_EdT_ms[i][0].at(k)->Add(summed_bp_gamma_EdT_ms[i][1].at(k),-1);
			summed_bp_gamma_EdT_us[i][0].at(k)->Add(summed_bp_gamma_EdT_us[i][1].at(k),-1);


			IsoDir->Append(beta_gamma_EdT_us[i][0].at(k));
			IsoDir->Append(beta_gamma_EdT_ms[i][0].at(k));
			IsoDir->Append(beta_gamma_EdT_s[i][0].at(k));
			IsoDir->Append(bp_gamma_EdT_us[i][0].at(k));
			IsoDir->Append(bp_gamma_EdT_ms[i][0].at(k));
			IsoDir->Append(bp_gamma_EdT_s[i][0].at(k));

			IsoDir->Append(summed_beta_gamma_EdT_us[i][0].at(k));
			IsoDir->Append(summed_beta_gamma_EdT_ms[i][0].at(k));
			IsoDir->Append(summed_beta_gamma_EdT_s[i][0].at(k));
			IsoDir->Append(summed_bp_gamma_EdT_us[i][0].at(k));
			IsoDir->Append(summed_bp_gamma_EdT_ms[i][0].at(k));
			IsoDir->Append(summed_bp_gamma_EdT_s[i][0].at(k));
			
			
			
		}//isotope_loop?
	}

	//std::cout<<"seg fault after here" <<std::endl;
	ofile->Write();
	//ofile->Close();

}//end of program
