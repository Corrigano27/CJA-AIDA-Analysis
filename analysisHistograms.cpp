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
			if ((*beta).Ey >= 0.0 && (*beta).Ex>=0.0){
				if (abs((*beta).Ex-(*beta).Ey)>=0){
					multix = (*beta).TFast & 0xFF;
					multiy = ((*beta).TFast >> 8) & 0xFF;
					for ( auto imp:(*beta).vectorOfImp ){ //if non-element gated histos needed, do here
						if ((*beta).z == (imp).Z){
							gammaVeto = true;
							In97gammaVeto = false;
							Ag95_160 = false;
							Ag95_800_1000 = false;
							Ag95_440 = false;
							Ag95_511 = false;
							Ag95_randomCheck = false;
							Cd96_gs = false;
							Cd96_m = false;
							Cd97_gs = false;
							Cd97_m = false;
							Cd97_n = false;
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
												if((*beta).T - anc.TIME < 40e3 && (anc.ID == 32 || anc.ID == 33)){
													if((*beta).T - anc.TIME > 0 && (anc.ID == 32 || anc.ID == 33)){
														betaVeto = true;
													}
												}
													
											}
											
											if (betaVeto == false){
												//use below to have variable dssd - will need to introduce further dssd vectors
												//if ((*beta).z >= isotopeDSSDStart[i].at(j) && (*beta).z <= isotopeDSSDEnd[i].at(j)){
												int DSSD = ((*beta).z);
												//if ((*beta).nx < 4 && (*beta).ny < 4 && (*beta).E<1500){
													//if (abs((imp).X-(*beta).x) < (0.5*multix + 0.5*((imp).TFAST &0xFF) +1.0)){
														//implantBeta[i][DSSD].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9);
													//}
												//}//end of upper beta energy cut
												//end of dssd if
												//end of dssd for
												decayEnergyAll[i].at(j)->Fill((*beta).E);
												if (multix >=0 && multiy >= 0){ //beta-delayed protons
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

																		//Cd-bp-isomers
																		if (elements[i] == "Cd" && isotopeStart[i]+j == 97){
																		//ground-state
																			if ((IndyGammaE > 590 && IndyGammaE < 783)||(IndyGammaE > 1300 && IndyGammaE < 1600)||(IndyGammaE > 1864 && IndyGammaE < 1992)){
																				Cd97_gs = true;
																			}
																		//25/2+ (m) isomer
																		//beta-gammas
																			if ((IndyGammaE > 76 && IndyGammaE < 123)||(IndyGammaE > 251 && IndyGammaE < 379)||(IndyGammaE > 724 && IndyGammaE < 906)||(IndyGammaE > 1148 && IndyGammaE < 1310)||(IndyGammaE > 1310 && IndyGammaE < 1573)||(IndyGammaE > 1681 && IndyGammaE < 1857)){
																				Cd97_m = true;
																			}
																		//one of beta-p gammas
																			if (IndyGammaE > 2000 && IndyGammaE < 2500){
																				Cd97_m = true;
																			}
																		//1/2- (n) isomer
																			if (IndyGammaE > 1024 && IndyGammaE < 1128){
																				Cd97_n = true;
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
																	}
																}
															}	
															
															
														}// end gamma loop
														
														if (ProtonGammaSumTemp != 0){
															summed_bp_gamma_EdT_s[i][0].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9, ProtonGammaSumTemp);
															summed_bp_gamma_EdT_ms[i][0].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e6, ProtonGammaSumTemp);
															summed_bp_gamma_EdT_us[i][0].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e3, ProtonGammaSumTemp);
														}

														if (ProtonGammaSumTempBg != 0){
															summed_bp_gamma_EdT_s[i][1].at(j)->Fill(-((*beta).T-(imp).TIME)/1.0e9, ProtonGammaSumTempBg);
															summed_bp_gamma_EdT_ms[i][1].at(j)->Fill(-((*beta).T-(imp).TIME)/1.0e6, ProtonGammaSumTempBg);
															summed_bp_gamma_EdT_us[i][1].at(j)->Fill(-((*beta).T-(imp).TIME)/1.0e3, ProtonGammaSumTempBg);

														}
														

														implantBeta1pAll[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9);

														delayed1pEnergyAll_AllDSSD[i].at(j)->Fill((*beta).E);
													}//end of lower beta-p energy cut
												}//end of beta-p multiplicity cut
	
												EdTAll_NoMultiGate[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9, (*beta).Ex);
												EdTAll_us[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).Ex);
												EdTAll_ms[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e6, (*beta).Ex);

												//ExEyAll[i].at(j)->Fill((*beta).Ex, (*beta).Ey);

												//DTAS 511/1022 coincidence check
												//gamma gated spectra

												//beta - DTAS correlations
												
												///////BETAS//////////

												if ((*beta).Ex<1100){
													if ((*beta).nx>=0 && (*beta).ny>=0){
														implantBetaAll[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9);
													}
													GammaSumTemp=0;
													GammaSumTempBg=0;
													//ProtonGammaSumTemp=0;
													//ProtonGammaSumTempBg=0;
													for ( auto gamma:(*beta).vectorOfGamma){ //loop over gamma events
														IndyGammaE = DTAS_SingleCalib(gamma.EN);
														if (((*beta).T-(gamma).TIME) > 10000){ //forward gammas
															if(((*beta).T-(gamma).TIME) < 20000){
																//summed_beta_gamma_1[i].at(j)->Fill((IndyGammaE));

																if ((gamma).ID == 777 && elements[i] == "Ag" && isotopeStart[i]+j == 95){
																	if ((*beta).Ex<1000 && (*beta).Ey<1000){
																		if (IndyGammaE > 1950 && IndyGammaE < 2200){
																			Ag95_EdT_2104keVsummed_gammaGated->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).Ex);
																			Ag95_EdT_2104keVsummed_gammaGated_back->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).Ey);

																			if (multix==0 && multiy==0){
																				Ag95_EdT_2104keVsummed_gammaGated_11->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).Ex);
																			}

																			Ag95_Implant_EdT_2104keVsummed_gammaGated->Fill(((*beta).T-(imp).TIME)/1.0e3, (imp).EN);
																			Ag95_EDiff_dT_2104keVsummed_gammaGated->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).Ex - (*beta).Ey);
																		}
																		if (IndyGammaE > 66 && IndyGammaE < 86){
																			Ag95_EdT_77keVsummed_gammaGated->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).Ex);

																		}
																	}
																}
																
																if ((gamma).ID<16){

																	if (elements[i] == "Cd" && isotopeStart[i]+j == 96){
																		//ground-state
																		if ((IndyGammaE > 388 && IndyGammaE < 456)||(IndyGammaE > 792 && IndyGammaE < 928)){
																			Cd96_gs = true;
																		}
																		//16+ (m) isomer
																		//beta-gammas
																		if ((IndyGammaE > 205 && IndyGammaE < 266)||(IndyGammaE > 595 && IndyGammaE < 764)||(IndyGammaE > 1268 && IndyGammaE < 1568)){
																			Cd96_m = true;
																		}
																		//one of beta-p gammas
																		if (IndyGammaE > 82 && IndyGammaE < 150){
																			Cd96_m = true;
																		}
																	}

																	if (elements[i] == "Cd" && isotopeStart[i]+j == 97){
																		//ground-state
																		if ((IndyGammaE > 590 && IndyGammaE < 783)||(IndyGammaE > 1300 && IndyGammaE < 1600)||(IndyGammaE > 1864 && IndyGammaE < 1992)){
																			Cd97_gs = true;
																		}
																		//25/2+ (m) isomer
																		//beta-gammas
																		if ((IndyGammaE > 76 && IndyGammaE < 123)||(IndyGammaE > 251 && IndyGammaE < 379)||(IndyGammaE > 724 && IndyGammaE < 906)||(IndyGammaE > 1148 && IndyGammaE < 1310)||(IndyGammaE > 1310 && IndyGammaE < 1573)||(IndyGammaE > 1681 && IndyGammaE < 1857)){
																			Cd97_m = true;
																		}
																		//one of beta-p gammas
																		if (IndyGammaE > 2000 && IndyGammaE < 2500){
																			Cd97_m = true;
																		}
																		//1/2- (n) isomer
																		if (IndyGammaE > 1024 && IndyGammaE < 1128){
																			Cd97_n = true;
																		}
																	}

																	if (elements[i] == "Ag" && isotopeStart[i]+j == 95){
																		if (IndyGammaE > 132 && IndyGammaE < 178){
																			Ag95_160 = true;
																			//Ag95_EdT_160keVgammaGated->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).Ex);
																		}
																		if (IndyGammaE > 787 && IndyGammaE < 1196){
																			Ag95_800_1000 = true;
																			//Ag95_EdT_800_1000keVgammaGated->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).Ex);
																		}
																		if (IndyGammaE > 400 && IndyGammaE < 445){
																			Ag95_440 = true;
																			//Ag95_EdT_440keVgammaGated->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).Ex);
																		}
																		if (IndyGammaE > 480 && IndyGammaE < 525){
																			Ag95_511 = true;
																			//Ag95_EdT_440keVgammaGated->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).Ex);
																		}
																		if (IndyGammaE > 600 && IndyGammaE < 645){
																			Ag95_randomCheck = true;
																			//Ag95_EdT_440keVgammaGated->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).Ex);
																		}

																	}

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

														if ((*beta).Ex<1000 && (*beta).Ey<1000 && elements[i] == "Ag"){
															
															for ( int d=0; d<(*beta).vectorOfGamma.size(); d++ ){ //loop over gamma events
																if (((*beta).T-(*beta).vectorOfGamma.at(d).TIME) > 10000){ //forward gammas
																	if(((*beta).T-(*beta).vectorOfGamma.at(d).TIME) < 20000){
																		if ((*beta).vectorOfGamma.at(d).ID<16){
																			if (isotopeStart[i]+j == 95){
																				if (((*beta).T - (imp).TIME)/1e6 > 0){
																					if (((*beta).T - (imp).TIME)/1e6 < 5){
																						Ag95_single_vs_summed_shorter->Fill(GammaSumTemp,DTAS_SingleCalib((*beta).vectorOfGamma.at(d).EN));
																						for ( int e=0; e<(*beta).vectorOfGamma.size(); e++ ){
																							if ( e != d && (*beta).vectorOfGamma.at(e).ID<16){
																								if (((*beta).T-(*beta).vectorOfGamma.at(e).TIME) > 10000){ //forward gammas
																									if(((*beta).T-(*beta).vectorOfGamma.at(e).TIME) < 20000){
																										Ag95_gamma_gamma_shorter->Fill(DTAS_SingleCalib((*beta).vectorOfGamma.at(d).EN), DTAS_SingleCalib((*beta).vectorOfGamma.at(e).EN));
																									}
																								}
																							}
																						}

																					}
																				}

																				if (((*beta).T - (imp).TIME)/1e6 > 10){
																					if (((*beta).T - (imp).TIME)/1e6 < 300){
																						Ag95_single_vs_summed->Fill(GammaSumTemp,DTAS_SingleCalib((*beta).vectorOfGamma.at(d).EN));
																						for ( int e=0; e<(*beta).vectorOfGamma.size(); e++ ){
																							if ( e != d && (*beta).vectorOfGamma.at(e).ID<16){
																								if (((*beta).T-(*beta).vectorOfGamma.at(e).TIME) > 10000){ //forward gammas
																									if(((*beta).T-(*beta).vectorOfGamma.at(e).TIME) < 20000){
																										Ag95_gamma_gamma->Fill(DTAS_SingleCalib((*beta).vectorOfGamma.at(d).EN), DTAS_SingleCalib((*beta).vectorOfGamma.at(e).EN));
																									}
																								}
																							}
																						}
																					}
																				}
																			}
																			
																		}
																	}
																}
															}
														}
													}

													if (GammaSumTempBg != 0){
														summed_beta_gamma_EdT_s[i][1].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9, GammaSumTempBg);
														summed_beta_gamma_EdT_ms[i][1].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e6, GammaSumTempBg);
														summed_beta_gamma_EdT_us[i][1].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e3, GammaSumTempBg);
													}
													

													if(Ag95_160 == true){Ag95_EdT_160keVgammaGated->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).Ex);}
													if(Ag95_440 == true){Ag95_EdT_440keVgammaGated->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).Ex);}
													if(Ag95_511 == true){Ag95_EdT_511keVgammaGated->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).Ex);}
													if(Ag95_800_1000 == true){Ag95_EdT_800_1000keVgammaGated->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).Ex);}
													if(Ag95_160 == true && Ag95_800_1000 == true){
														if (GammaSumTemp > 2000){
															if (GammaSumTemp < 2200){
																Ag95_EdT_allpeaks_gammaGated->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).Ex);
															}
														}
													}
													if(Ag95_randomCheck == true){Ag95_EdT_randomcheck->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).Ex);}

														
												}//beta energy and multi-cut
												if (Cd96_gs == true){Cd96_gs_EdT->Fill(((*beta).T-(imp).TIME)/1.0e9, (*beta).Ex);}
												if (Cd96_m == true){Cd96_m_EdT->Fill(((*beta).T-(imp).TIME)/1.0e9, (*beta).Ex);}
												if (Cd97_gs == true){Cd97_gs_EdT->Fill(((*beta).T-(imp).TIME)/1.0e9, (*beta).Ex);}
												if (Cd97_m == true){Cd97_m_EdT->Fill(((*beta).T-(imp).TIME)/1.0e9, (*beta).Ex);}
												if (Cd97_n == true){Cd97_n_EdT->Fill(((*beta).T-(imp).TIME)/1.0e9, (*beta).Ex);}



											
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

		if ((*implant).T){
			for ( auto pid:(*implant).vectorOfPid ){
				PID_implant->Fill((*implant).aoq, (*implant).zet);
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

	PID_implant->Write();


	Ag95_EdT_160keVgammaGated->Write();
	Ag95_EdT_800_1000keVgammaGated->Write();
	Ag95_EdT_440keVgammaGated->Write();
	Ag95_EdT_511keVgammaGated->Write();
	Ag95_EdT_randomcheck->Write();

	Ag95_EdT_160_800_1000keVgammaGated->Add(Ag95_EdT_160keVgammaGated, 1);
	Ag95_EdT_160_800_1000keVgammaGated->Add(Ag95_EdT_800_1000keVgammaGated, 1);

	Ag95_EdT_160_800_1000keVgammaGated->Write();
	
	Ag95_EdT_allpeaks_gammaGated->Write();

	Ag95_EdT_2104keVsummed_gammaGated->Write();

	Ag95_EdT_2104keVsummed_gammaGated_back->Write();

	Ag95_EdT_2104keVsummed_gammaGated_11->Write();

	Ag95_EdT_77keVsummed_gammaGated->Write();

	Ag95_EDiff_dT_2104keVsummed_gammaGated->Write();

	Ag95_Implant_EdT_2104keVsummed_gammaGated->Write();

	Ag95_single_vs_summed->Write();

	Ag95_single_vs_summed_shorter->Write();

	Ag95_gamma_gamma->Write();

	Ag95_gamma_gamma_shorter->Write();

	Cd96_gs_EdT->Write();
	Cd96_m_EdT->Write();
	Cd97_gs_EdT->Write();
	Cd97_m_EdT->Write();
	Cd97_n_EdT->Write();

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
			IsoDir->Append(delayed1pEnergyAll_AllDSSD[i].at(k));
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