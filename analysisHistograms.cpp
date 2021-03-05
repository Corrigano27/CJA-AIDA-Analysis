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
//#include "/home/corrigan/DTAS_Merger/merger/MergerSoft/data2Tree.cxx"
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

	//boolean variable for beta vetoes

	bool betaVeto;

	bool gammaVeto;

	uint8_t multix = 0;
	uint8_t multiy = 0;
	
	//Files read, histograms filled
	while (aReader.Next()){

		if ((*beta).T){
			if ((*beta).Ey >= 200.0 && (*beta).Ex>=200.0){
				multix = (*beta).TFast & 0xFF;
				multiy = ((*beta).TFast >> 8) & 0xFF;
				for ( auto imp:(*beta).vectorOfImp ){ //if non-element gated histos needed, do here
					if ((*beta).z == (imp).Z){
						gammaVeto = true;
						//gamma gated spectra
						for ( auto gamma:(*beta).vectorOfGamma ){ //loop over gamma events
						//forward implant-decay events
						//ExEyAll_gammaloop[i].at(j)->Fill((*beta).Ex, (*beta).Ey);
							if (((*beta).T-(gamma).TIME) > 10000 && (gamma).ID==777){ //forward gammas
								if(((*beta).T-(gamma).TIME) < 20000 && (gamma).ID==777){
									if ((gamma).EN>1000){
										gammaVeto = false;
									}
								}
							}
						}
						EdT_global->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).E);	
						EdT_global_longer->Fill(((*beta).T-(imp).TIME)/1.0e6, (*beta).E);


						if (((*beta).T-(imp).TIME > 0)){
							PID->Fill((imp).AOQ, (imp).ZET);
							ExEy_global->Fill((*beta).Ex, (*beta).Ey);
							EDiff_global->Fill((*beta).Ex - (*beta).Ey);
							if (((*beta).T-(imp).TIME)/1e3 < 100){
								PID_noise->Fill((imp).AOQ, (imp).ZET);
								XY_Hits->Fill((*beta).x, (*beta).y);
								FastDSSD->Fill((*beta).z);
							}
						}

						if (gammaVeto == false){ 					
							EdT_global_gammagate->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).E);
							EdT_global_longer_gammagate->Fill(((*beta).T-(imp).TIME)/1.0e6, (*beta).E);
							if (((*beta).T-(imp).TIME > 0)){
								if (((*beta).T-(imp).TIME)/1e3 < 100){
									PID_gamma->Fill((imp).AOQ, (imp).ZET);
								//EDiff_global_gammagate->Fill((*beta).Ex - (*beta).Ey);
								}
							}
						}
						for (int i = 0; i < numElements; i++){
							for (int j = 0; j <= isotopeEnd[i]-isotopeStart[i]; j++){
								if(particleCuts[i][j]->IsInside((imp).AOQ,(imp).ZET)){
									if ((*beta).z >= isotopeDSSDStart[i].at(j) && (*beta).z <= isotopeDSSDEnd[i].at(j)){// use these statements for the dssd loop later on
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
											if((*beta).T - anc.TIME < 40e3 && (anc.ID == 32 || anc.ID == 33)){
												if((*beta).T - anc.TIME > 0 && (anc.ID == 32 || anc.ID == 33)){
													//betaVeto = true;
												}
											}
												
										} // start looping gammas here
										
										if (betaVeto == false){
											//use below to have variable dssd - will need to introduce further dssd vectors
											//if ((*beta).z >= isotopeDSSDStart[i].at(j) && (*beta).z <= isotopeDSSDEnd[i].at(j)){
											int DSSD = ((*beta).z);
											decayEnergy[i][DSSD].at(j)->Fill((*beta).E);

											if ((*beta).Ex>1400 && (*beta).Ey>1400){
												EdT_bp[i][DSSD].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e6, (*beta).E);
											}
											ExEyDiff[i][DSSD].at(j)->Fill((*beta).Ex - (*beta).Ey);

											if (elements[i] == "Ag" && isotopeStart[i]+j == 95){
												if ((*beta).Ex>300){
													if ((*beta).Ex<400){
														if (((*beta).T - (imp).TIME)/1e6 > 0){
															if (((*beta).T - (imp).TIME)/1e6 < 5){
																if (gammaVeto == false){
																	EDiff_global_gammagate->Fill((*beta).Ex - (*beta).Ey);
																	XY_Hits_gammagate->Fill((*beta).x, (*beta).y);
																	FastDSSD_gammagate->Fill((*beta).z);
																	NxNy_global_gammagate->Fill((*beta).nx, (*beta).ny);
																	CxCy_global_gammagate->Fill(multix, multiy);
																}
															}
														}
						
													}
												}
											}

											if (multix == 0 && multiy == 0){
												EdT[i][DSSD].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e6, (*beta).E);
												if ((*beta).Ex>0 && (*beta).Ey>0){
													if (((*beta).T-(imp).TIME > 0)){
														if ((((*beta).T-(imp).TIME)/1e6 < 1)){
															delayed1pEnergyX_0_1_ms[i][DSSD].at(j)->Fill((*beta).Ex);
														}
														if (((*beta).T-(imp).TIME)/1e9 < 0.1){
															delayed1pEnergyX_0_100_ms[i][DSSD].at(j)->Fill((*beta).Ex);
														}
														if ((((*beta).T-(imp).TIME)/1e9 < 1)){
															delayed1pEnergyX_0_1_s[i][DSSD].at(j)->Fill((*beta).Ex);
														}
													}
												}

												if ((*beta).Ex>1400 && (*beta).Ey>1400){
													if (((*beta).T-(imp).TIME > 0)){
														delayed1pEnergy[i][DSSD].at(j)->Fill((*beta).E);
														delayed1pEnergyX[i][DSSD].at(j)->Fill((*beta).Ex);
														delayed1pEnergyY[i][DSSD].at(j)->Fill((*beta).Ey);
														
														//ExEy[i][z].at(j)->Fill((*beta).Ex, (*beta).Ey);
														//EnergyXChannel[i][z].at(j)->Fill((*beta).x, (*beta).E);
														//EnergyYChannel[i][z].at(j)->Fill((*beta).y, (*beta).E);

													}

													implantBeta1p[i][DSSD].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9);

													/*if (elements[i] == "Ag" && isotopeStart[i]+j == 94 && (*beta).z==2){
														if ((*beta).Ex>1700){
															if ((*beta).Ex<2000){
																Ag94ImplantBeta1p_DSSD2_smallpeak->Fill(((*beta).T-(imp).TIME)/1.0e9);
																for ( auto gamma:(*beta).vectorOfGamma ){ //loop over gamma events
																	if (((*beta).T-(gamma).TIME) > 20000){
																		if(((*beta).T-(gamma).TIME) < 30000){
																			if ((gamma).ID==777){
																				Ag94_peak_Gamma777_Bg->Fill((gamma.EN));
																			}
																			if ((gamma).ID<16){		
																				Ag94_peak_GammaSingle_Bg->Fill((gamma.EN));
																			}

																		}
																	}


																	if (((*beta).T-(gamma).TIME) > 10000){
																		if(((*beta).T-(gamma).TIME) < 20000){
																			
															
																			if ((gamma).ID==777){
																				Ag94_peak_Gamma777->Fill((gamma.EN));
																			}
																			else if ((gamma).ID<16){		
																				Ag94_peak_GammaSingle->Fill((gamma.EN));
																			}
																			
																		}
																	}
																}//end gamma loop




															}
														}
														if ((*beta).Ex>2000){
															Ag94ImplantBeta1p_DSSD2_largepeak->Fill(((*beta).T-(imp).TIME)/1.0e9);
														}
													
													}
													if (elements[i] == "Sn" && isotopeStart[i]+j == 101 && (*beta).z==1){
														if ((*beta).Ex>2150){
															if ((*beta).Ex<2400){
																Sn101ImplantBeta1p_DSSD1_smallpeak->Fill(((*beta).T-(imp).TIME)/1.0e9);
																for ( auto gamma:(*beta).vectorOfGamma ){ //loop over gamma events
																	if (((*beta).T-(gamma).TIME) > 20000){
																		if(((*beta).T-(gamma).TIME) < 30000){
																			if ((gamma).ID==777){
																				Sn101_peak_Gamma777_Bg->Fill((gamma.EN));
																			}
																			if ((gamma).ID<16){		
																				Sn101_peak_GammaSingle_Bg->Fill((gamma.EN));
																			}

																		}
																	}
																	

																	if (((*beta).T-(gamma).TIME) > 10000){
																		if(((*beta).T-(gamma).TIME) < 20000){
																			
															
																			if ((gamma).ID==777){
																				Sn101_peak_Gamma777->Fill((gamma.EN));
																			}
																			if ((gamma).ID<16){		
																				Sn101_peak_GammaSingle->Fill((gamma.EN));
																			}
																			
																		}
																	}
																}//end gamma loop




															}
														}
														if ((*beta).Ex>2400){
															Sn101ImplantBeta1p_DSSD1_largepeak->Fill(((*beta).T-(imp).TIME)/1.0e9);
														}
													
													}*/
													if ((*beta).T-(imp).TIME < 0){
														delayed1pEnergyRandom[i][DSSD].at(j)->Fill((*beta).E);
													}

													delayed1pEnergyAll[i][DSSD].at(j)->Fill((*beta).E);
												}//end of lower beta-p energy cut
											
											}
											else if (multix < 3 && multiy < 3 && (*beta).E<1500){
												implantBeta[i][DSSD].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9);
										
											}//end of upper beta energy cut
											//end of dssd if
											//end of dssd for
											decayEnergyAll[i].at(j)->Fill((*beta).E);
											if (multix == 0 && multiy == 0){
												EdTAll11[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e6, (*beta).E);

												/*if (elements[i] == "In" && isotopeStart[i]+j == 97){
													if (((*beta).T-(imp).TIME)/1e6 > 0){
														if (((*beta).T-(imp).TIME)/1e6 < 100){
															In97_GroundStateE->Fill((*beta).Ex);
															for ( auto gamma:(*beta).vectorOfGamma ){
																if (((*beta).T-(gamma).TIME) > 20000){
																	if(((*beta).T-(gamma).TIME) < 30000){
																		if ((*beta).Ex<5000 && (*beta).Ey<5000){

																			if ((gamma).ID==777){
																				In97m_Gamma777_Bg->Fill((gamma.EN));
																			}
																			if ((gamma).ID<16){		
																				In97m_GammaSingle_Bg->Fill((gamma.EN));
																		
																			}
																		}

																	}
																}

																if (((*beta).T-(gamma).TIME) > 10000){
																	if(((*beta).T-(gamma).TIME) < 20000){
																		if ((*beta).Ex<5000 && (*beta).Ey<5000){
																		
																			if ((gamma).ID==777){
																				In97m_Gamma777->Fill((gamma.EN));
																			}
																			if ((gamma).ID<16){		
																				In97m_GammaSingle->Fill((gamma.EN));
																			}
																		}
																	}
																}
															}//end of gamma loop
														}


													}
												}*/

												if ((*beta).Ex>1400 && (*beta).Ey>1400){
													if (((*beta).T-(imp).TIME > 0)){
														delayed1pEnergy_AllDSSD[i].at(j)->Fill((*beta).E);
													}
														//ExEyAll[i].at(j)->Fill((*beta).Ex, (*beta).Ey);
														//ExEyDiffAll[i].at(j)->Fill((*beta).Ex - (*beta).Ey);
														//EnergyXChannelAll[i].at(j)->Fill((*beta).x, (*beta).E);
														//EnergyYChannelAll[i].at(j)->Fill((*beta).y, (*beta).E);
													for ( auto gamma:(*beta).vectorOfGamma ){ //loop over gamma events
														if (((*beta).T-(imp).TIME > 0)){ //forward implant-decay events
															if (((*beta).T-(gamma).TIME) > 10000){ //forward gammas
																if(((*beta).T-(gamma).TIME) < 20000){
																	if ((gamma).ID==777){
																		summed_bp_gamma_1[i].at(j)->Fill((gamma.EN));
																	}
																	else if ((gamma).ID<16){		
																		bp_gamma_1[i].at(j)->Fill((gamma.EN));
																	}
																}
															}
															if (((*beta).T-(gamma).TIME) > 20000){ //random gammas
																if(((*beta).T-(gamma).TIME) < 30000){
																	if ((gamma).ID==777){
																		summed_bp_gamma_2[i].at(j)->Fill((gamma.EN));
																	}
																	if ((gamma).ID<16){		
																		bp_gamma_2[i].at(j)->Fill((gamma.EN));
																	}
																}
															}	
														}//end forward implant-time if
														if (((*beta).T-(imp).TIME < 0)){ //backward implant-decay events
															if (((*beta).T-(gamma).TIME) > 10000){ //forward gammas
																if(((*beta).T-(gamma).TIME) < 20000){
																	if ((gamma).ID==777){
																		summed_bp_gamma_3[i].at(j)->Fill((gamma.EN));
																	}
																	else if ((gamma).ID<16){		
																		bp_gamma_3[i].at(j)->Fill((gamma.EN));
																	}
																}
															}
															if (((*beta).T-(gamma).TIME) > 20000){ //random gammas
																if(((*beta).T-(gamma).TIME) < 30000){
																	if ((gamma).ID==777){
																		summed_bp_gamma_4[i].at(j)->Fill((gamma.EN));
																	}
																	if ((gamma).ID<16){		
																		bp_gamma_4[i].at(j)->Fill((gamma.EN));
																	}
																}
															}	
														}//end backward implant-time if
													}// end gamma loop		
													

													implantBeta1pAll[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9);

													if ((*beta).T-(imp).TIME < 0){
														delayed1pEnergyRandom_AllDSSD[i].at(j)->Fill((*beta).E);
													}

													delayed1pEnergyAll_AllDSSD[i].at(j)->Fill((*beta).E);
												}//end of lower beta-p energy cut
											}//end of beta-p multiplicity cut

											if (multix == 0 && multiy == 1){
												EdTAll12[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e6, (*beta).E);
											}
											if (multix == 1 && multiy == 0){
												EdTAll21[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e6, (*beta).E);
											}
											if (multix == 1 && multiy == 1){
												EdTAll22[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e6, (*beta).E);
											}		
											EdTAll_NoMultiGate[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e6, (*beta).E);
											EdTAll[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).E);
											ExEyAll[i].at(j)->Fill((*beta).Ex, (*beta).Ey);

											gammaVeto = true;
											//gamma gated spectra
											for ( auto gamma:(*beta).vectorOfGamma ){ //loop over gamma events
												//forward implant-decay events
												//ExEyAll_gammaloop[i].at(j)->Fill((*beta).Ex, (*beta).Ey);
												if (((*beta).T-(gamma).TIME) > 10000 && (gamma).ID==777){ //forward gammas
													if(((*beta).T-(gamma).TIME) < 20000 && (gamma).ID==777){
														if ((gamma).EN>1000){
															gammaVeto = false;
														}
													}
												}
											}	
											if (gammaVeto == false){ 					
												EdT_gammagate[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).E);
												EdT_gammagate_longer[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e6, (*beta).E);
												if (((*beta).T-(imp).TIME > 0)){
													ExEy_gammagate[i].at(j)->Fill((*beta).Ex, (*beta).Ey);
													EDiff_gammagate[i].at(j)->Fill((*beta).Ex - (*beta).Ey);
													NxNy_gammagate[i].at(j)->Fill((*beta).nx, (*beta).ny);
													clustersize_gammagate[i].at(j)->Fill(multix, multiy);
												}
											}
																
			


											if (multix < 3 && multiy < 3 && (*beta).Ex>300){
												if ((*beta).Ex<400){
													implantBetaAll[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9);
													for ( auto gamma:(*beta).vectorOfGamma ){ //loop over gamma events
														if (((*beta).T-(imp).TIME > 0)){ //forward implant-decay events
															if ((((*beta).T-(imp).TIME)/1e6 < 5)){
																if (((*beta).T-(gamma).TIME) > 10000){ //forward gammas
																	if(((*beta).T-(gamma).TIME) < 20000){
																		if ((gamma).ID==777){
																			summed_beta_gamma_1[i].at(j)->Fill((gamma.EN));
																		}
																		else if ((gamma).ID<16){		
																			beta_gamma_1[i].at(j)->Fill((gamma.EN));
																		}
																	}
																}
																if (((*beta).T-(gamma).TIME) > 20000){ //random gammas
																	if(((*beta).T-(gamma).TIME) < 30000){
																		if ((gamma).ID==777){
																			summed_beta_gamma_2[i].at(j)->Fill((gamma.EN));
																		}
																		if ((gamma).ID<16){		
																			beta_gamma_2[i].at(j)->Fill((gamma.EN));
																		}
																	}
																}
															}	
														}//end forward implant-time if
														if (((*beta).T-(imp).TIME < 0)){ //backward implant-decay events
															if ((((*beta).T-(imp).TIME)/1e6 > -5)){
																if (((*beta).T-(gamma).TIME) > 10000){ //forward gammas
																	if(((*beta).T-(gamma).TIME) < 20000){
																		if ((gamma).ID==777){
																			summed_beta_gamma_3[i].at(j)->Fill((gamma.EN));
																		}
																		else if ((gamma).ID<16){		
																			beta_gamma_3[i].at(j)->Fill((gamma.EN));
																		}
																	}
																}
																if (((*beta).T-(gamma).TIME) > 20000){ //random gammas
																	if(((*beta).T-(gamma).TIME) < 30000){
																		if ((gamma).ID==777){
																			summed_beta_gamma_4[i].at(j)->Fill((gamma.EN));
																		}
																		if ((gamma).ID<16){		
																			beta_gamma_4[i].at(j)->Fill((gamma.EN));
																		}
																	}
																}
															}	
														}//end backward implant-time if
			
													}//end of gamma loop		
												}

													
											}
										
										}//end of beta veto application
										

									}//end of stopping layer if statement
								}//end of particle cut if statement

							}//end of isotope for loop
						}//end of elements for loop
					}//end of imp = decay dssd if
				}//end of loop over correlated events
			} //end of if beta events with positive energy

		}//end of loop through beta events

		if ((*implant).T){
			for ( auto pid:(*implant).vectorOfPid ){
				PID_implant->Fill((*implant).aoq, (*implant).zet);
				int iDSSD = (*implant).z;
				for (int i =0; i<numElements; i++){
					for (int j = 0; j <= isotopeEnd[i]-isotopeStart[i]; j++){	
						if (particleCuts[i][j]->IsInside((*implant).aoq, (*implant).zet)){
							implantZ[i].at(j)->Fill((*implant).z);
							implantVelocityimplantZ[i].at(j)->Fill((pid).VELOCITY, (*implant).z);
							implantEAll[i].at(j)->Fill((*implant).E);
							implantVelocityAOQ_AllDSSD[i].at(j)->Fill((*implant).aoq,(pid).VELOCITY);
							implantEnergyAOQ_AllDSSD[i].at(j)->Fill((*implant).aoq,(*implant).E);
							implantE[i][iDSSD].at(j)->Fill((*implant).E);
							implantVelocityimplantE[i][iDSSD].at(j)->Fill((pid).VELOCITY, (*implant).E);
							implantVelocityAOQ[i][iDSSD].at(j)->Fill((*implant).aoq,(pid).VELOCITY);
							implantEnergyAOQ[i][iDSSD].at(j)->Fill((*implant).aoq,(*implant).E);
								
							
						}
					}
				}
			}
			
		}


	}//end of loop through chain

	ofile->cd();

	std::cout << "Writing to file" << std::endl;

	PID->Write();

	PID_noise->Write();

	PID_gamma->Write();

	PID_implant->Write();

	/*Sn101ImplantBeta1p_DSSD1_smallpeak->Write();

	Sn101ImplantBeta1p_DSSD1_largepeak->Write();

	Sn101_peak_GammaSingle->Write();

	Sn101_peak_Gamma777->Write();

	Sn101_peak_GammaSingle_Bg->Write();

	Sn101_peak_Gamma777_Bg->Write();

	Ag94ImplantBeta1p_DSSD2_smallpeak->Write();

	Ag94ImplantBeta1p_DSSD2_largepeak->Write();

	Ag94_peak_GammaSingle->Write();

	Ag94_peak_Gamma777->Write();

	Ag94_peak_GammaSingle_Bg->Write();

	Ag94_peak_Gamma777_Bg->Write();

	In97_GroundStateE->Write();

	In97m_GammaSingle->Write();

	In97m_Gamma777->Write();

	In97m_GammaSingle_Bg->Write();

	In97m_Gamma777_Bg->Write();*/

	XY_Hits->Write();

	XY_Hits_gammagate->Write();

	FastDSSD->Write();

	FastDSSD_gammagate->Write();

	EdT_global_longer_gammagate->Write();

	EdT_global_gammagate->Write();

	EDiff_global_gammagate->Write();

	NxNy_global_gammagate->Write();

	CxCy_global_gammagate->Write();

	EdT_global_longer->Write();

	EdT_global->Write();

	EDiff_global->Write();

	ExEy_global->Write();

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

			for(int z = 0; z < 6; z++){

				IsoDir->Append(decayEnergy[i][z].at(k));
				IsoDir->Append(delayed1pEnergy[i][z].at(k));
				IsoDir->Append(delayed1pEnergyX[i][z].at(k));
				IsoDir->Append(delayed1pEnergyX_0_1_ms[i][z].at(k));
				IsoDir->Append(delayed1pEnergyX_0_100_ms[i][z].at(k));
				IsoDir->Append(delayed1pEnergyX_0_1_s[i][z].at(k));
				IsoDir->Append(delayed1pEnergyY[i][z].at(k));
				IsoDir->Append(delayed1pEnergyRandom[i][z].at(k));
				IsoDir->Append(delayed1pEnergyAll[i][z].at(k));
			//	ExEy[i][z].at(k)->Write(
			//	ExEyDiff[i][z].at(k)->Write();		
			//	EnergyXChannel[i][z].at(k)->Write();
			//	EnergyYChannel[i][z].at(k)->Write()
				IsoDir->Append(implantBeta[i][z].at(k));
				IsoDir->Append(implantBeta1p[i][z].at(k));
				IsoDir->Append(implantE[i][z].at(k));
				IsoDir->Append(implantVelocityimplantE[i][z].at(k));
				IsoDir->Append(EdT[i][z].at(k));
				IsoDir->Append(EdT_bp[i][z].at(k));
				IsoDir->Append(ExEyDiff[i][z].at(k));
				IsoDir->Append(implantVelocityAOQ[i][z].at(k));
				IsoDir->Append(implantEnergyAOQ[i][z].at(k));
			
			}
			//combined DSSD
			IsoDir->Append(decayEnergyAll[i].at(k));			
			IsoDir->Append(delayed1pEnergy_AllDSSD[i].at(k));						
			IsoDir->Append(delayed1pEnergyRandom_AllDSSD[i].at(k));
			IsoDir->Append(delayed1pEnergyAll_AllDSSD[i].at(k));
			//ExEyAll[i].at(k)->Writ
			//ExEyDiffAll[i].at(k)->Write();
			//EnergyXChannelAll[i].at(k)->Write();
			//EnergyYChannelAll[i].at(k)->Write
			IsoDir->Append(implantBetaAll[i].at(k));
			IsoDir->Append(implantBeta1pAll[i].at(k));
			IsoDir->Append(implantEAll[i].at(k));
			IsoDir->Append(implantVelocityimplantEAll[i].at(k));
			IsoDir->Append(EdTAll[i].at(k));
			IsoDir->Append(EdTAll_NoMultiGate[i].at(k));
			IsoDir->Append(EdTAll11[i].at(k));
			IsoDir->Append(EdTAll12[i].at(k));
			IsoDir->Append(EdTAll21[i].at(k));
			IsoDir->Append(EdTAll22[i].at(k));
			IsoDir->Append(implantVelocityAOQ_AllDSSD[i].at(k));
			IsoDir->Append(implantZ[i].at(k));
			IsoDir->Append(implantEnergyAOQ_AllDSSD[i].at(k));

			IsoDir->Append(ExEyAll[i].at(k));
			//IsoDir->Append(ExEyAll_gammaloop[i].at(k));
			IsoDir->Append(EdT_gammagate[i].at(k));
			IsoDir->Append(EdT_gammagate_longer[i].at(k));
			IsoDir->Append(EDiff_gammagate[i].at(k));
			IsoDir->Append(ExEy_gammagate[i].at(k));
			IsoDir->Append(NxNy_gammagate[i].at(k));
			IsoDir->Append(clustersize_gammagate[i].at(k));

			//gamma spectra correction

			/**beta_gamma_corr[i].at(k) = (*beta_gamma_1[i].at(k)-*beta_gamma_2[i].at(k)) - (*beta_gamma_3[i].at(k)-*beta_gamma_4[i].at(k));
			*summed_beta_gamma_corr[i].at(k) = (*summed_beta_gamma_1[i].at(k)-*summed_beta_gamma_2[i].at(k)) - (*summed_beta_gamma_3[i].at(k)-*summed_beta_gamma_4[i].at(k));
			*bp_gamma_corr[i].at(k) = (*bp_gamma_1[i].at(k)-*bp_gamma_2[i].at(k)) - (*bp_gamma_3[i].at(k)-*bp_gamma_4[i].at(k));
			*summed_bp_gamma_corr[i].at(k) = (*summed_bp_gamma_1[i].at(k)-*summed_bp_gamma_2[i].at(k)) - (*summed_bp_gamma_3[i].at(k)-*summed_bp_gamma_4[i].at(k));
			*/
			beta_gamma_corr[i].at(k)->Add(beta_gamma_1[i].at(k),1);
			beta_gamma_corr[i].at(k)->Add(beta_gamma_2[i].at(k),-1);
			beta_gamma_corr[i].at(k)->Add(beta_gamma_3[i].at(k),-1);
			beta_gamma_corr[i].at(k)->Add(beta_gamma_4[i].at(k),1);

			summed_beta_gamma_corr[i].at(k)->Add(summed_beta_gamma_1[i].at(k),1);
			summed_beta_gamma_corr[i].at(k)->Add(summed_beta_gamma_2[i].at(k),-1);
			summed_beta_gamma_corr[i].at(k)->Add(summed_beta_gamma_3[i].at(k),-1);
			summed_beta_gamma_corr[i].at(k)->Add(summed_beta_gamma_4[i].at(k),1);

			bp_gamma_corr[i].at(k)->Add(bp_gamma_1[i].at(k),1);
			bp_gamma_corr[i].at(k)->Add(bp_gamma_2[i].at(k),-1);
			bp_gamma_corr[i].at(k)->Add(bp_gamma_3[i].at(k),-1);
			bp_gamma_corr[i].at(k)->Add(bp_gamma_4[i].at(k),1);

			summed_bp_gamma_corr[i].at(k)->Add(summed_bp_gamma_1[i].at(k),1);
			summed_bp_gamma_corr[i].at(k)->Add(summed_bp_gamma_2[i].at(k),-1);
			summed_bp_gamma_corr[i].at(k)->Add(summed_bp_gamma_3[i].at(k),-1);
			summed_bp_gamma_corr[i].at(k)->Add(summed_bp_gamma_4[i].at(k),1);

			IsoDir->Append(beta_gamma_corr[i].at(k));
			IsoDir->Append(summed_beta_gamma_corr[i].at(k));
			IsoDir->Append(bp_gamma_corr[i].at(k));
			IsoDir->Append(summed_bp_gamma_corr[i].at(k));

			IsoDir->Append(implantVelocityimplantZ[i].at(k));
			
		}//isotope_loop?
	}


	ofile->Write();
	ofile->Close();

}