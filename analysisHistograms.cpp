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

//#include "/Disk/ds-sopa-personal/s1333561/PhD/MergerSoftware/data2Tree.cxx"
#include "/home/corrigan/AidaSoftware/MergerSoftware/data2Tree.cxx"
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
												betaVeto = true;
											}
										}
											
									} // start looping gammas here
									
									if (betaVeto == false){
										//use below to have variable dssd - will need to introduce further dssd vectors
										//if ((*beta).z >= isotopeDSSDStart[i].at(j) && (*beta).z <= isotopeDSSDEnd[i].at(j)){
										int DSSD = ((*beta).z);
										decayEnergy[i][DSSD].at(j)->Fill((*beta).E);
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
													//ExEyDiff[i][z].at(j)->Fill((*beta).Ex - (*beta).Ey);
													//EnergyXChannel[i][z].at(j)->Fill((*beta).x, (*beta).E);
													//EnergyYChannel[i][z].at(j)->Fill((*beta).y, (*beta).E);

												}

												implantBeta1p[i][DSSD].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9);

												if (elements[i] == "Ag" && isotopeStart[i]+j == 94 && (*beta).z==2){
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
												
												}
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
											EdTAll[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e6, (*beta).E);

											if (elements[i] == "In" && isotopeStart[i]+j == 97){
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
											}

											if ((*beta).Ex>1400 && (*beta).Ey>1400){
												if (((*beta).T-(imp).TIME > 0)){
													delayed1pEnergy_AllDSSD[i].at(j)->Fill((*beta).E);
													//ExEyAll[i].at(j)->Fill((*beta).Ex, (*beta).Ey);
													//ExEyDiffAll[i].at(j)->Fill((*beta).Ex - (*beta).Ey);
													//EnergyXChannelAll[i].at(j)->Fill((*beta).x, (*beta).E);
													//EnergyYChannelAll[i].at(j)->Fill((*beta).y, (*beta).E);
													for ( auto gamma:(*beta).vectorOfGamma ){ //loop over gamma events
														if (((*beta).T-(gamma).TIME) > 20000){
															if(((*beta).T-(gamma).TIME) < 30000){
																if ((gamma).ID==777){
																	bp_gammaTest777_Bg[i].at(j)->Fill((gamma.EN));
																}
																if ((gamma).ID<16){		
																	bp_gammaTest_Bg[i].at(j)->Fill((gamma.EN));
																}

															}
														}

														if (((*beta).T-(gamma).TIME) > 10000){
															if(((*beta).T-(gamma).TIME) < 20000){
																if ((gamma).ID==777){
																	bp_gammaTest777[i].at(j)->Fill((gamma.EN));
																}
																else if ((gamma).ID<16){		
																	bp_gammaTest[i].at(j)->Fill((gamma.EN));
																}
															}
														}


													}//end gamma loop

													

														
												}

												implantBeta1pAll[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9);

												if ((*beta).T-(imp).TIME < 0){
													delayed1pEnergyRandom_AllDSSD[i].at(j)->Fill((*beta).E);
												}

												delayed1pEnergyAll_AllDSSD[i].at(j)->Fill((*beta).E);
											}//end of lower beta-p energy cut
										}//end of beta-p multiplicity cut		
										
										else if (multix < 3 && multiy < 3 && (*beta).E<1400){
											implantBetaAll[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9);
											for ( auto gamma:(*beta).vectorOfGamma ){ //loop over gamma events
												if (((*beta).T-(gamma).TIME) > 20000){
													if(((*beta).T-(gamma).TIME) < 30000){
														if ((gamma).ID==777){
															gammaTest777_Bg[i].at(j)->Fill((gamma.EN));
														}
														if ((gamma).ID<16){		
															gammaTest_Bg[i].at(j)->Fill((gamma.EN));
														}
														
													}
												}


												if (((*beta).T-(gamma).TIME) > 10000){
													if(((*beta).T-(gamma).TIME) < 20000){
														if ((gamma).ID==777){
															gammaTest777[i].at(j)->Fill((gamma.EN));
														}
														if ((gamma).ID<16){		
															gammaTest[i].at(j)->Fill((gamma.EN));
														}
														
													}
												}
												
											}//end of gamma loop		
										

											
										}
									
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
			for ( auto pid:(*implant).vectorOfPid ){
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

	Sn101ImplantBeta1p_DSSD1_smallpeak->Write();

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

	In97m_Gamma777_Bg->Write();

	std::string isoDirName;

	for(int i = 0; i < numElements; i++){

		ElDir = ofile->mkdir(elements[i].c_str());
		ElementDir.push_back(ElDir);
		ofile->cd(elements[i].c_str());
		//element-directory

		for (int k = 0; k <= isotopeEnd[i]-isotopeStart[i]; k++){

			isoDirName = elements[i].c_str() + std::to_string(isotopeStart[i] + k);
			IsoDir = ofile->mkdir(isoDirName.c_str());
			IsotopeDir.push_back(IsoDir);

			for(int z = 0; z < 6; z++){

				//for(unsigned int k = 0; k < decayEnergy[i][z].size(); k++){
				ElementDir.at(i)->Append(decayEnergy[i][z].at(k));
				//}
				/*for(unsigned int k = 0; k < delayed1pEnergy[i][z].size(); k++){
					ElementDir.at(i)->Append(delayed1pEnergy[i][z].at(k));
				}
				for(unsigned int k = 0; k < delayed1pEnergyX[i][z].size(); k++){
					ElementDir.at(i)->Append(delayed1pEnergyX[i][z].at(k));
				}
				for(unsigned int k = 0; k < delayed1pEnergyX_0_1_ms[i][z].size(); k++){
					ElementDir.at(i)->Append(delayed1pEnergyX_0_1_ms[i][z].at(k));
				}
				for(unsigned int k = 0; k < delayed1pEnergyX_0_100_ms[i][z].size(); k++){
					ElementDir.at(i)->Append(delayed1pEnergyX_0_100_ms[i][z].at(k));
				}
				for(unsigned int k = 0; k < delayed1pEnergyX_0_1_s[i][z].size(); k++){
					ElementDir.at(i)->Append(delayed1pEnergyX_0_1_s[i][z].at(k));
				}
				for(unsigned int k = 0; k < delayed1pEnergyY[i][z].size(); k++){
					ElementDir.at(i)->Append(delayed1pEnergyY[i][z].at(k));
				}
				for(unsigned int k = 0; k < delayed1pEnergyRandom[i][z].size(); k++){
					ElementDir.at(i)->Append(delayed1pEnergyRandom[i][z].at(k));
				}
				for(unsigned int k = 0; k < delayed1pEnergyAll[i][z].size(); k++){
					ElementDir.at(i)->Append(delayed1pEnergyAll[i][z].at(k));
				}
				//for(unsigned int k = 0; k < ExEy[i][z].size(); k++){
				//	ExEy[i][z].at(k)->Write();
				//}
				//for(unsigned int k = 0; k < ExEyDiff[i][z].size(); k++){
				//	ExEyDiff[i][z].at(k)->Write();
				//}
				//for(unsigned int k = 0; k < EnergyXChannel[i][z].size(); k++){
				//	EnergyXChannel[i][z].at(k)->Write();
				//}
				//for(unsigned int k = 0; k < EnergyYChannel[i][z].size(); k++){
				//	EnergyYChannel[i][z].at(k)->Write();
				//}
				for(unsigned int k = 0; k < implantBeta[i][z].size(); k++){
					ElementDir.at(i)->Append(implantBeta[i][z].at(k));
				}
				for(unsigned int k = 0; k < implantBeta1p[i][z].size(); k++){
					ElementDir.at(i)->Append(implantBeta1p[i][z].at(k));
				}
				for(unsigned int k = 0; k < implantE[i][z].size(); k++){
					ElementDir.at(i)->Append(implantE[i][z].at(k));
				}
				for(unsigned int k = 0; k < implantVelocityimplantE[i][z].size(); k++){
					ElementDir.at(i)->Append(implantVelocityimplantE[i][z].at(k));
				}
				for(unsigned int k = 0; k < EdT[i][z].size(); k++){
					ElementDir.at(i)->Append(EdT[i][z].at(k));
				}
				for(unsigned int k = 0; k < implantVelocityAOQ[i][z].size(); k++){
					ElementDir.at(i)->Append(implantVelocityAOQ[i][z].at(k));
				}
				for(unsigned int k = 0; k < implantEnergyAOQ[i][z].size(); k++){
					ElementDir.at(i)->Append(implantEnergyAOQ[i][z].at(k));
				}
			}
			//combined dssds
			for(unsigned int k = 0; k < decayEnergyAll[i].size(); k++){
				ElementDir.at(i)->Append(decayEnergyAll[i].at(k));
			}
			for(unsigned int k = 0; k < delayed1pEnergy_AllDSSD[i].size(); k++){
				ElementDir.at(i)->Append(delayed1pEnergy_AllDSSD[i].at(k));
			}
			for(unsigned int k = 0; k < delayed1pEnergyRandom_AllDSSD[i].size(); k++){
				ElementDir.at(i)->Append(delayed1pEnergyRandom_AllDSSD[i].at(k));
			}
			for(unsigned int k = 0; k < delayed1pEnergyAll_AllDSSD[i].size(); k++){
				ElementDir.at(i)->Append(delayed1pEnergyAll_AllDSSD[i].at(k));
			}
			//for(unsigned int k = 0; k < ExEyAll[i].size(); k++){
				//ExEyAll[i].at(k)->Write();
			//}
			//for(unsigned int k = 0; k < ExEyDiffAll[i].size(); k++){
				//ExEyDiffAll[i].at(k)->Write();
			//}
			//for(unsigned int k = 0; k < EnergyXChannelAll[i].size(); k++){
			//	EnergyXChannelAll[i].at(k)->Write();
			//}
			//for(unsigned int k = 0; k < EnergyYChannelAll[i].size(); k++){
			//	EnergyYChannelAll[i].at(k)->Write();
			//}
			for(unsigned int k = 0; k < implantBetaAll[i].size(); k++){
				ElementDir.at(i)->Append(implantBetaAll[i].at(k));
			}
			for(unsigned int k = 0; k < implantBeta1pAll[i].size(); k++){
				ElementDir.at(i)->Append(implantBeta1pAll[i].at(k));
			}
			for(unsigned int k = 0; k < implantEAll[i].size(); k++){
				ElementDir.at(i)->Append(implantEAll[i].at(k));
			}
			for(unsigned int k = 0; k < implantVelocityimplantEAll[i].size(); k++){
				ElementDir.at(i)->Append(implantVelocityimplantEAll[i].at(k));
			}
			for(unsigned int k = 0; k < EdTAll[i].size(); k++){
				ElementDir.at(i)->Append(EdTAll[i].at(k));
			}
			for(unsigned int k = 0; k < implantVelocityAOQ_AllDSSD[i].size(); k++){
				ElementDir.at(i)->Append(implantVelocityAOQ_AllDSSD[i].at(k));
			}
			for(unsigned int k = 0; k < implantEnergyAOQ_AllDSSD[i].size(); k++){
				ElementDir.at(i)->Append(implantEnergyAOQ_AllDSSD[i].at(k));
			}
			//gamma stuff
			for(unsigned int k = 0; k < gammaTest[i].size(); k++){
				ElementDir.at(i)->Append(gammaTest[i].at(k));
			}
			for(unsigned int k = 0; k < gammaTest777[i].size(); k++){
				ElementDir.at(i)->Append(gammaTest777[i].at(k));
			}
			for(unsigned int k = 0; k < bp_gammaTest[i].size(); k++){
				ElementDir.at(i)->Append(bp_gammaTest[i].at(k));
			}
			for(unsigned int k = 0; k < bp_gammaTest777[i].size(); k++){
				ElementDir.at(i)->Append(bp_gammaTest777[i].at(k));
			}
			for(unsigned int k = 0; k < gammaTest_Bg[i].size(); k++){
				ElementDir.at(i)->Append(gammaTest_Bg[i].at(k));
			}
			for(unsigned int k = 0; k < gammaTest777_Bg[i].size(); k++){
				ElementDir.at(i)->Append(gammaTest777_Bg[i].at(k));
			}
			for(unsigned int k = 0; k < bp_gammaTest_Bg[i].size(); k++){
				ElementDir.at(i)->Append(bp_gammaTest_Bg[i].at(k));
			}
			for(unsigned int k = 0; k < bp_gammaTest777_Bg[i].size(); k++){
				ElementDir.at(i)->Append(bp_gammaTest777_Bg[i].at(k));
			}

			// implant stuff
			for(unsigned int k = 0; k < implantZ[i].size(); k++){
				ElementDir.at(i)->Append(implantZ[i].at(k));
			}
			for(unsigned int k = 0; k < implantVelocityimplantZ[i].size(); k++){
				ElementDir.at(i)->Append(implantVelocityimplantZ[i].at(k));
			}*/
			}
		}//isotope_loop?
	}


	ofile->Write();
	ofile->Close();

}