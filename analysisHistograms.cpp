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
			if ((*beta).Ey >= 200.0 && (*beta).Ex>=200.0){
				if (abs((*beta).Ex-(*beta).Ey)<120){
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
							/*Ag96_470 = false;
							Ag96_743 = false;
							Ag96_1249 = false;
							Ag96_Random = false;
							Ag96_subtraction = false;
							Ag96_isomer = false;*/
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
														//betaVeto = true;
													}
												}

												//F11 veto (beta)
												if((*beta).T - anc.TIME < 40e3 && (anc.ID == 32 || anc.ID == 33)){
													if((*beta).T - anc.TIME > 0 && (anc.ID == 32 || anc.ID == 33)){
														//betaVeto = true;
													}
												}
													
											}
											
											if (betaVeto == false){
												//use below to have variable dssd - will need to introduce further dssd vectors
												//if ((*beta).z >= isotopeDSSDStart[i].at(j) && (*beta).z <= isotopeDSSDEnd[i].at(j)){
												int DSSD = ((*beta).z);
												decayEnergy[i][DSSD].at(j)->Fill((*beta).E);
												//ExEyDiff[i][DSSD].at(j)->Fill((*beta).Ex - (*beta).Ey);

												if (multix == 0 && multiy == 0){
													EdT[i][DSSD].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e6, (*beta).E);
													
													if ((*beta).Ex>1400 && (*beta).Ey>1400){
														if (((*beta).T-(imp).TIME > 0)){
															//delayed1pEnergy[i][DSSD].at(j)->Fill((*beta).E);
															delayed1pEnergyX[i][DSSD].at(j)->Fill((*beta).Ex);
															//delayed1pEnergyY[i][DSSD].at(j)->Fill((*beta).Ey);
															
															//ExEy[i][z].at(j)->Fill((*beta).Ex, (*beta).Ey);
															//EnergyXChannel[i][z].at(j)->Fill((*beta).x, (*beta).E);
															//EnergyYChannel[i][z].at(j)->Fill((*beta).y, (*beta).E);

														}

														implantBeta1p[i][DSSD].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9);

														if ((*beta).T-(imp).TIME < 0){
															delayed1pEnergyRandom[i][DSSD].at(j)->Fill((*beta).E);
														}

														delayed1pEnergyAll[i][DSSD].at(j)->Fill((*beta).E);
													}//end of lower beta-p energy cut
												
												}
												if ((*beta).nx < 4 && (*beta).ny < 4 && (*beta).E<1500){
													if (abs((imp).X-(*beta).x) < (0.5*multix + 0.5*((imp).TFAST &0xFF) +1.0)){
														implantBeta[i][DSSD].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9);
													}
												}//end of upper beta energy cut
												//end of dssd if
												//end of dssd for
												decayEnergyAll[i].at(j)->Fill((*beta).E);
												if (multix == 0 && multiy == 0){
													if ((*beta).Ex>1400 && (*beta).Ey>1400){
														if (((*beta).T-(imp).TIME > 0)){
															delayed1pEnergy_AllDSSD[i].at(j)->Fill((*beta).E);
														}
														//beta-p gamma loop
														ProtonGammaSumTemp = 0;
														//ProtonGammaSumTempBg = 0;
														Sn101Counter = 0;
														Sn101CounterBg = 0;
														Sn101Counter_Pk = 0;
														Sn101CounterBg_Pk = 0;
														Ag94_Peak_Counter = 0;
														Ag94_Peak_CounterBg = 0;
														for ( auto gamma:(*beta).vectorOfGamma ){ //loop over gamma events
															IndyGammaE = DTAS_SingleCalib(gamma.EN);
															if (((*beta).T-(imp).TIME > 0)){ //forward implant-decay events
																if (((*beta).T-(gamma).TIME) > 10000){ //forward gammas
																	if (((*beta).T-(gamma).TIME) < 20000){
																		if ((gamma).ID<16){		
																			//bp_gamma_1[i].at(j)->Fill((gamma.EN));
																			bp_gamma_EdT_s[i][0].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9, IndyGammaE);
																			bp_gamma_EdT_ms[i][0].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e6, IndyGammaE);
																			bp_gamma_EdT_us[i][0].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e3, IndyGammaE);
																			ProtonGammaSumTemp+=(IndyGammaE);
																			if (elements[i] == "Sn" && isotopeStart[i]+j == 101 && (*beta).z==1){
																				if ((((*beta).T-(imp).TIME)/1e9 < 5)){
																					if ((*beta).Ex>2150){
																						if ((*beta).Ex<2400){
																							Tin101_bp_gamma_peak[0]->Fill(IndyGammaE);
																							Sn101Counter_Pk += IndyGammaE;
																							
																						}
																					}
																					if ((*beta).Ex>2400){
																						Tin101_bp_gamma_rest[0]->Fill(IndyGammaE);
																						Sn101Counter += IndyGammaE;
																					}
																				}
																			}
																			if (elements[i] == "Ag" && isotopeStart[i]+j == 94){
																				if ((*beta).Ex>1711){
																					if ((*beta).Ex<1940){
																						Ag94_1800_bp_DTASindy[0]->Fill(IndyGammaE);
																						Ag94_Peak_Counter += IndyGammaE;
																					}
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
																			if (elements[i] == "Sn" && isotopeStart[i]+j == 101 && (*beta).z==1){
																				if ((((*beta).T-(imp).TIME)/1e9 < 5)){
																					if ((*beta).Ex>2150){
																						if ((*beta).Ex<2400){
																							Tin101_bp_gamma_peak[1]->Fill((IndyGammaE));
																						}
																					}
																					if ((*beta).Ex>2400){
																						Tin101_bp_gamma_rest[1]->Fill((IndyGammaE));
																					}
																				}
																			}
																			if (elements[i] == "Ag" && isotopeStart[i]+j == 94){
																				if ((*beta).Ex>1711){
																					if ((*beta).Ex<1940){
																						Ag94_1800_bp_DTASindy[1]->Fill((IndyGammaE));
																					}
																				}
																			}
																		}
																	}
																}	
															}//end forward implant-time if
															if (((*beta).T-(imp).TIME < 0)){ //backward implant-decay events
																if (((*beta).T-(gamma).TIME) > 10000){ //forward gammas
																	if(((*beta).T-(gamma).TIME) < 20000){
																		if ((gamma).ID<16){		
																			//bp_gamma_3[i].at(j)->Fill((IndyGammaE));
																			bp_gamma_EdT_s[i][2].at(j)->Fill(-((*beta).T-(imp).TIME)/1.0e9, (IndyGammaE));
																			bp_gamma_EdT_ms[i][2].at(j)->Fill(-((*beta).T-(imp).TIME)/1.0e6, (IndyGammaE));
																			bp_gamma_EdT_us[i][2].at(j)->Fill(-((*beta).T-(imp).TIME)/1.0e3, (IndyGammaE));
																			ProtonGammaSumTemp+=(IndyGammaE);

																			if (elements[i] == "Sn" && isotopeStart[i]+j == 101 && (*beta).z==1){
																				if ((((*beta).T-(imp).TIME)/1e9 > -5)){
																					if ((*beta).Ex>2150){
																						if ((*beta).Ex<2400){
																							Tin101_bp_gamma_peak[2]->Fill((IndyGammaE));
																							Sn101CounterBg_Pk += (IndyGammaE);
																						}
																					}
																					if ((*beta).Ex>2400){
																						Tin101_bp_gamma_rest[2]->Fill((IndyGammaE));
																						Sn101CounterBg += (IndyGammaE);
																					}
																				}
																			}
																			if (elements[i] == "Ag" && isotopeStart[i]+j == 94){
																				if ((*beta).Ex>1711){
																					if ((*beta).Ex<1940){
																						Ag94_1800_bp_DTASindy[2]->Fill((IndyGammaE));
																						Ag94_Peak_CounterBg += (IndyGammaE);
																					}
																				}
																			}
																		}
																	}
																}
																if (((*beta).T-(gamma).TIME) > 20000){ //random gammas
																	if(((*beta).T-(gamma).TIME) < 30000){
																		if ((gamma).ID<16){		
																			//bp_gamma_4[i].at(j)->Fill((IndyGammaE));
																			bp_gamma_EdT_s[i][3].at(j)->Fill(-((*beta).T-(imp).TIME)/1.0e9, (IndyGammaE));
																			bp_gamma_EdT_ms[i][3].at(j)->Fill(-((*beta).T-(imp).TIME)/1.0e6, (IndyGammaE));
																			bp_gamma_EdT_us[i][3].at(j)->Fill(-((*beta).T-(imp).TIME)/1.0e3, (IndyGammaE));
																			if (elements[i] == "Sn" && isotopeStart[i]+j == 101 && (*beta).z==1){
																				if ((((*beta).T-(imp).TIME)/1e9 > -5)){
																					if ((*beta).Ex>2150){
																						if ((*beta).Ex<2400){
																							Tin101_bp_gamma_peak[3]->Fill((IndyGammaE));
																						}
																					}
																					if ((*beta).Ex>2400){
																						Tin101_bp_gamma_rest[3]->Fill((IndyGammaE));
																					}
																				
																				}
																			}
																			if (elements[i] == "Ag" && isotopeStart[i]+j == 94){
																				if ((*beta).Ex>1711){
																					if ((*beta).Ex<1940){
																						Ag94_1800_bp_DTASindy[3]->Fill((IndyGammaE));
																					}
																				}
																			}
																															
																		}
																	}
																}		
															}//end backward implant-time if
														}// end gamma loop
														
														if (ProtonGammaSumTemp != 0){
															summed_bp_gamma_EdT_s[i][0].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9, ProtonGammaSumTemp);
															summed_bp_gamma_EdT_ms[i][0].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e6, ProtonGammaSumTemp);
															summed_bp_gamma_EdT_us[i][0].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e3, ProtonGammaSumTemp);

															//summed_p_gamma_E_p_E[i][0].at(j)->Fill((*beta).Ex, ProtonGammaSumTemp);
														
			
															//summed_bp_gamma_EdT_s[i][0].at(j)->Fill(-((*beta).T-(imp).TIME)/1.0e9, ProtonGammaSumTemp);
															//summed_bp_gamma_EdT_ms[i][0].at(j)->Fill(-((*beta).T-(imp).TIME)/1.0e6, ProtonGammaSumTemp);
															//summed_bp_gamma_EdT_us[i][0].at(j)->Fill(-((*beta).T-(imp).TIME)/1.0e3, ProtonGammaSumTemp);

															//summed_p_gamma_E_p_E[i][1].at(j)->Fill((*beta).Ex, ProtonGammaSumTempBg);

														}
														
														if (Sn101Counter_Pk != 0){
															Tin101_summed_bp_gamma_peak[0]->Fill(Sn101Counter_Pk);
														}
														if (Sn101CounterBg_Pk != 0){
															Tin101_summed_bp_gamma_peak[1]->Fill(Sn101CounterBg_Pk);
														}
														if (Sn101Counter != 0){
															Tin101_summed_bp_gamma_rest[0]->Fill(Sn101Counter);
														}
														if (Sn101CounterBg != 0){
															Tin101_summed_bp_gamma_rest[1]->Fill(Sn101CounterBg);
														}
														if (Ag94_Peak_Counter != 0){
															Ag94_1800_bp_DTASsummed[0]->Fill(Ag94_Peak_Counter);
														}
														if (Ag94_Peak_CounterBg != 0){
															Ag94_1800_bp_DTASsummed[1]->Fill(Ag94_Peak_CounterBg);
														}

														implantBeta1pAll[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9);

														if ((*beta).T-(imp).TIME < 0){
															delayed1pEnergyRandom_AllDSSD[i].at(j)->Fill((*beta).E);
														}

														delayed1pEnergyAll_AllDSSD[i].at(j)->Fill((*beta).E);
													}//end of lower beta-p energy cut
												}//end of beta-p multiplicity cut

												if ((*beta).nx == 1 && (*beta).ny == 1){
													EdTAll11[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9, (*beta).E);
												}
												if ((*beta).nx == 1 && (*beta).ny == 2){
													EdTAll12[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9, (*beta).Ex);
												}
												if ((*beta).nx == 2 && (*beta).ny == 1){
													EdTAll21[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9, (*beta).Ex);
												}
												if ((*beta).nx == 2 && (*beta).ny == 2){
													EdTAll22[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9, (*beta).Ex);
												}		
												EdTAll_NoMultiGate[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9, (*beta).Ex);
												EdTAll_us[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).Ex);
												EdTAll_ms[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e6, (*beta).Ex);

												if ((*beta).T - (imp).TIME > 0){
													EdTAll_NoMultiGate_corr[i][0].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9, (*beta).Ex);
													EdTAll_ms_corr[i][0].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e6, (*beta).Ex);
													EdTAll_us_corr[i][0].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).Ex);
												}

												if ((*beta).T - (imp).TIME < 0){
													EdTAll_NoMultiGate_corr[i][1].at(j)->Fill(-((*beta).T-(imp).TIME)/1.0e9, (*beta).Ex);
													EdTAll_ms_corr[i][1].at(j)->Fill(-((*beta).T-(imp).TIME)/1.0e6, (*beta).Ex);
													EdTAll_us_corr[i][1].at(j)->Fill(-((*beta).T-(imp).TIME)/1.0e3, (*beta).Ex);
												}
												

												//ExEyAll[i].at(j)->Fill((*beta).Ex, (*beta).Ey);

												//DTAS 511/1022 coincidence check
												gammaVeto = true;
												//gamma gated spectra
												for ( auto gamma:(*beta).vectorOfGamma){ //loop over gamma events
													IndyGammaE = DTAS_SingleCalib(gamma.EN);;
													/*if (elements[i] == "Ag" && isotopeStart[i]+j == 96){
														if ((gamma).ID<16){		
															if ((*beta).Ex<1000 && (*beta).Ey<1000){
																if ((IndyGammaE > 450 && IndyGammaE < 550)||(IndyGammaE > 720 && IndyGammaE < 850)||(IndyGammaE > 1200 && IndyGammaE < 1400)){
																	Ag96_GammaT_betaT_all3Peaks->Fill(((*beta).T-(imp).TIME)/1e3, ((*beta).T-(gamma).TIME)/1e3);
																}
															}
														}
													}*/

													//forward implant-decay events
													//ExEyAll_gammaloop[i].at(j)->Fill((*beta).Ex, (*beta).Ey);
													if (((*beta).T-(gamma).TIME) > 10000){ //forward gammas
														if(((*beta).T-(gamma).TIME) < 20000){
															if ((gamma).ID==777){
																In97gammaVeto = true;
																//if (IndyGammaE>1000){
																//	gammaVeto = false;
																//}
															}
														}
													}
												}	
												/*if (gammaVeto == false){ 					
													EdT_gammagate[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).E);
													EdT_gammagate_longer[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e6, (*beta).E);
													if (((*beta).T-(imp).TIME > 0)){
														ExEy_gammagate[i].at(j)->Fill((*beta).Ex, (*beta).Ey);
														EDiff_gammagate[i].at(j)->Fill((*beta).Ex - (*beta).Ey);
														NxNy_gammagate[i].at(j)->Fill((*beta).nx, (*beta).ny);
														clustersize_gammagate[i].at(j)->Fill(multix, multiy);
													}
												}*/

												//beta - DTAS correlations
												if (elements[i] == "In" && isotopeStart[i]+j == 97 && In97gammaVeto == true){
													Indium97_gammaveto_EdT->Fill(((*beta).T-(imp).TIME)/1.0e6, (*beta).E);
												}

												if ((*beta).Ex<1500 && (*beta).Ey<1500){
													implantBetaAll[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9);
													GammaSumTemp=0;
													//GammaSumTempBg=0;
													ProtonGammaSumTemp=0;
													//ProtonGammaSumTempBg=0;
													for ( auto gamma:(*beta).vectorOfGamma ){ //loop over gamma events
													IndyGammaE = DTAS_SingleCalib(gamma.EN);
														if (((*beta).T-(imp).TIME > 0)){ //forward implant-decay events
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

																		if (elements[i] == "Ag" && isotopeStart[i]+j == 95){
																			if ((*beta).Ex<1000 && (*beta).Ey<1000){
																				beta_gamma_EdT_s[i][0].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9, (IndyGammaE));
																				beta_gamma_EdT_ms[i][0].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e6, (IndyGammaE));
																				beta_gamma_EdT_us[i][0].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e3, (IndyGammaE));
																				GammaSumTemp+=(IndyGammaE);
																			}
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

																		else{
																			//beta_gamma_1[i].at(j)->Fill((IndyGammaE));
																			beta_gamma_EdT_s[i][0].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9, (IndyGammaE));
																			beta_gamma_EdT_ms[i][0].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e6, (IndyGammaE));
																			beta_gamma_EdT_us[i][0].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e3, (IndyGammaE));
																			GammaSumTemp+=(IndyGammaE);
																			if (multix==0 && multiy==0 && (*beta).Ex<1400 && (*beta).Ey<1400){
																				ProtonGammaSumTemp+=(IndyGammaE);
																			}
																		}
																		/*if (elements[i] == "Ag" && isotopeStart[i]+j == 96){
																			if ((*beta).Ex<1000 && (*beta).Ey<1000){
																				Ag96_isomer = true;
																				if (IndyGammaE > 450 && IndyGammaE < 500){
																					Ag96_470 = true;
																				}
																				if (IndyGammaE > 720 && IndyGammaE < 760){
																					Ag96_743 = true;
																				}
																				if (IndyGammaE > 1200 && IndyGammaE < 1300){
																					Ag96_1249 = true;
																				}
																			}
													
																		}*/

																	}
																	
																}
															}
															if (((*beta).T-(gamma).TIME) > 20000){ //random gammas
																if(((*beta).T-(gamma).TIME) < 30000){
																	if ((gamma).ID<16){	

																		if (elements[i] == "Ag" && isotopeStart[i]+j == 95){
																			if ((*beta).Ex<1000 && (*beta).Ey<1000){
																				beta_gamma_EdT_s[i][1].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9, (IndyGammaE));
																				beta_gamma_EdT_ms[i][1].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e6, (IndyGammaE));
																				beta_gamma_EdT_us[i][1].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e3, (IndyGammaE));
																			}

																		}

																		else{
																			//beta_gamma_2[i].at(j)->Fill((IndyGammaE));
																			beta_gamma_EdT_s[i][1].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9, (IndyGammaE));
																			beta_gamma_EdT_ms[i][1].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e6, (IndyGammaE));
																			beta_gamma_EdT_us[i][1].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e3, (IndyGammaE));
																		}
																		/*if (elements[i] == "Ag" && isotopeStart[i]+j == 96){
																			if ((*beta).Ex<1000 && (*beta).Ey<1000){
																				Ag96_Random = true;
																			}
																		}*/

																		

																	}
																	
																}
															}
														
														}//end forward implant-time if
														if (((*beta).T-(imp).TIME < 0)){ //backward implant-decay events
															if (((*beta).T-(gamma).TIME) > 10000){ //forward gammas
																if(((*beta).T-(gamma).TIME) < 20000){

																	if ((gamma).ID == 777 && elements[i] == "Ag" && isotopeStart[i]+j == 95){
																		if ((*beta).Ex<1000 && (*beta).Ey<1000){
																			if (IndyGammaE > 1950 && IndyGammaE < 2200){
																				Ag95_EdT_2104keVsummed_gammaGated->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).Ex);
																				Ag95_EdT_2104keVsummed_gammaGated_back->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).Ey);

																				if ((*beta).nx==1 && (*beta).ny==1){
																					Ag95_EdT_2104keVsummed_gammaGated_11->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).Ex);
																				}
																				
																				Ag95_Implant_EdT_2104keVsummed_gammaGated->Fill(((*beta).T-(imp).TIME)/1.0e3, (imp).EN);
																				Ag95_EDiff_dT_2104keVsummed_gammaGated->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).Ex - (*beta).Ey);
																			}
																			if (IndyGammaE > 67 && IndyGammaE < 87){
																				Ag95_EdT_77keVsummed_gammaGated->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).Ex);

																			}
																		}
																	}

																	if ((gamma).ID<16){
																		if (elements[i] == "Ag" && isotopeStart[i]+j == 95){

																			if ((*beta).Ex<1000 && (*beta).Ey<1000){
																				beta_gamma_EdT_s[i][2].at(j)->Fill(-((*beta).T-(imp).TIME)/1.0e9, (IndyGammaE));
																				beta_gamma_EdT_ms[i][2].at(j)->Fill(-((*beta).T-(imp).TIME)/1.0e6, (IndyGammaE));
																				beta_gamma_EdT_us[i][2].at(j)->Fill(-((*beta).T-(imp).TIME)/1.0e3, (IndyGammaE));
																				GammaSumTemp += (IndyGammaE);
																			}
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
																		else{	
																			//beta_gamma_3[i].at(j)->Fill((IndyGammaE));
																			beta_gamma_EdT_s[i][2].at(j)->Fill(-((*beta).T-(imp).TIME)/1.0e9, (IndyGammaE));
																			beta_gamma_EdT_ms[i][2].at(j)->Fill(-((*beta).T-(imp).TIME)/1.0e6, (IndyGammaE));
																			beta_gamma_EdT_us[i][2].at(j)->Fill(-((*beta).T-(imp).TIME)/1.0e3, (IndyGammaE));
																			GammaSumTemp += (IndyGammaE);
																			if (multix==0 && multiy==0 && (*beta).Ex<1400 && (*beta).Ey<1400){
																				ProtonGammaSumTemp+=(IndyGammaE);
																			}
																		}
																		/*if (elements[i] == "Ag" && isotopeStart[i]+j == 96){
																			if ((*beta).Ex<1000 && (*beta).Ey<1000){
																				Ag96_isomer = true;
																				if (IndyGammaE > 450 && IndyGammaE < 500){
																					Ag96_470 = true;
																				}
																				if (IndyGammaE > 720 && IndyGammaE < 760){
																					Ag96_743 = true;
																				}
																				if (IndyGammaE > 1200 && IndyGammaE < 1300){
																					Ag96_1249 = true;
																				}
																			}
													
																		}*/

																	}
																	
																}
															}
															if (((*beta).T-(gamma).TIME) > 20000){ //random gammas
																if(((*beta).T-(gamma).TIME) < 30000){
																	if ((gamma).ID<16){

																		if (elements[i] == "Ag" && isotopeStart[i]+j == 95){
																			if ((*beta).Ex<1000 && (*beta).Ey<1000){

																				beta_gamma_EdT_s[i][3].at(j)->Fill(-((*beta).T-(imp).TIME)/1.0e9, (IndyGammaE));
																				beta_gamma_EdT_ms[i][3].at(j)->Fill(-((*beta).T-(imp).TIME)/1.0e6, (IndyGammaE));
																				beta_gamma_EdT_us[i][3].at(j)->Fill(-((*beta).T-(imp).TIME)/1.0e3, (IndyGammaE));
																			}
																		}
																		else{		
																		//beta_gamma_4[i].at(j)->Fill((IndyGammaE));
																			beta_gamma_EdT_s[i][3].at(j)->Fill(-((*beta).T-(imp).TIME)/1.0e9, (IndyGammaE));
																			beta_gamma_EdT_ms[i][3].at(j)->Fill(-((*beta).T-(imp).TIME)/1.0e6, (IndyGammaE));
																			beta_gamma_EdT_us[i][3].at(j)->Fill(-((*beta).T-(imp).TIME)/1.0e3, (IndyGammaE));
																		}
																		/*if (elements[i] == "Ag" && isotopeStart[i]+j == 96){
																			if ((*beta).Ex<1000 && (*beta).Ey<1000){
																				Ag96_Random = true;
																			}
																		}*/

																	}
																}
															}	
														}//end backward implant-time if
														
												
			
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
																			/*if (isotopeStart[i]+j == 96){
																				if (((*beta).T - (imp).TIME)/1e3 > 0){
																					if (((*beta).T - (imp).TIME)/1e3 < 500){
																						Ag96_single_vs_summed->Fill(GammaSumTemp,DTAS_SingleCalib((*beta).vectorOfGamma.at(d).EN));
																						for ( int e=0; e<(*beta).vectorOfGamma.size(); e++ ){
																							if ( e != d && (*beta).vectorOfGamma.at(e).ID<16){		
																								if (((*beta).T-(*beta).vectorOfGamma.at(e).TIME) > 10000){ //forward gammas
																									if(((*beta).T-(*beta).vectorOfGamma.at(e).TIME) < 20000){																	
																										Ag96_gamma_gamma->Fill(DTAS_SingleCalib((*beta).vectorOfGamma.at(d).EN), DTAS_SingleCalib((*beta).vectorOfGamma.at(e).EN));
																										gammaSubtract = DTAS_SingleCalib((*beta).vectorOfGamma.at(d).EN) + DTAS_SingleCalib((*beta).vectorOfGamma.at(e).EN);
																										if (GammaSumTemp != 0){
																											if (DTAS_SingleCalib((*beta).vectorOfGamma.at(d).EN)>430 && DTAS_SingleCalib((*beta).vectorOfGamma.at(d).EN)<490){
																												if (DTAS_SingleCalib((*beta).vectorOfGamma.at(e).EN)>700 && DTAS_SingleCalib((*beta).vectorOfGamma.at(e).EN)<806){
																													Ag96_sum_E1E2_diff_470_740->Fill(GammaSumTemp - gammaSubtract + 100);
																												}
																											}
																											if (DTAS_SingleCalib((*beta).vectorOfGamma.at(d).EN)>700 && DTAS_SingleCalib((*beta).vectorOfGamma.at(d).EN)<770){
																												if (DTAS_SingleCalib((*beta).vectorOfGamma.at(e).EN)>1176 && DTAS_SingleCalib((*beta).vectorOfGamma.at(e).EN)<1319){
																													Ag96_sum_E1E2_diff_740_1249->Fill(GammaSumTemp - gammaSubtract + 100);
																												}
																											}
																										}
																										
																									}
																								}
																							}
																						}
																					}
																				}
																			}*/
																		}
																	}
																}
															}
														}
													}
													//if(ProtonGammaSumTemp != 0){
														//summed_p_gamma_E_p_E[i][0].at(j)->Fill((*beta).Ex, ProtonGammaSumTemp);
													//}

													if (GammaSumTemp != 0){
														summed_beta_gamma_EdT_s[i][0].at(j)->Fill(-((*beta).T-(imp).TIME)/1.0e9, GammaSumTemp);
														summed_beta_gamma_EdT_ms[i][0].at(j)->Fill(-((*beta).T-(imp).TIME)/1.0e6, GammaSumTemp);
														summed_beta_gamma_EdT_us[i][0].at(j)->Fill(-((*beta).T-(imp).TIME)/1.0e3, GammaSumTemp);

														//summed_beta_gamma_E_beta_E[i][1].at(j)->Fill((*beta).E, GammaSumTempBg);
													}
													//if(ProtonGammaSumTempBg != 0){
														//summed_p_gamma_E_p_E[i][1].at(j)->Fill((*beta).Ex, ProtonGammaSumTempBg);
													//}

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

													/*if (Ag96_470 == true){Ag96_EdT_470keVgammaGated->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).Ex);}
													if (Ag96_743 == true){Ag96_EdT_743keVgammaGated->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).Ex);}
													if (Ag96_1249 == true){Ag96_EdT_1249keVgammaGated->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).Ex);}
													if (Ag96_470 == true || Ag96_743 == true || Ag96_1249 == true){
														Ag96_EdT_all3Peaks_gammaGated->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).Ex);
														Ag96_E_correlatedGamma->Fill((*beta).E);
													}
													if(Ag96_Random == true){
														Ag96_EdT_all3Peaks_Random_gammaGated->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).Ex);
														Ag96_E_randomGamma->Fill((*beta).E);
													}

													if(Ag96_isomer == true){
														if ((GammaSumTemp>0 && GammaSumTemp<946) || (GammaSumTemp>1053 && GammaSumTemp<1171) || (GammaSumTemp>1284 && GammaSumTemp<1628) || (GammaSumTemp>1773 && GammaSumTemp<1914) || (GammaSumTemp>2042 && GammaSumTemp<2319)){
															Ag96_EdT_summed_gammaGated->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).Ex);
														}
														if (GammaSumTemp>2400 && GammaSumTemp<2550){
															Ag96_EdT_2461keVgammaGated->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).Ex);
														}
													}*/
														
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

	PID_implant->Write();

	Indium97_gammaveto_EdT->Write();

	for(int g = 0; g < 4; g++){
		if (g == 0 || g == 3){
			Tin101_bp_gamma_peak_corr->Add(Tin101_bp_gamma_peak[g],1);
			//Tin101_summed_bp_gamma_peak_corr->Add(Tin101_summed_bp_gamma_peak[g],1);
			Tin101_bp_gamma_rest_corr->Add(Tin101_bp_gamma_rest[g],1);
			//Tin101_summed_bp_gamma_rest_corr->Add(Tin101_summed_bp_gamma_rest[g],1);
			
			
		}

		if (g == 1 || g == 2){
			Tin101_bp_gamma_peak_corr->Add(Tin101_bp_gamma_peak[g],-1);
			//Tin101_summed_bp_gamma_peak_corr->Add(Tin101_summed_bp_gamma_peak[g],-1);
			Tin101_bp_gamma_rest_corr->Add(Tin101_bp_gamma_rest[g],-1);
			//Tin101_summed_bp_gamma_rest_corr->Add(Tin101_summed_bp_gamma_rest[g],-1);
			
		}
		Ag94_1800_bp_DTASindy[g]->Write();
		if (g < 2){
			Ag94_1800_bp_DTASsummed[g]->Write();
		}
				
	}

	Tin101_summed_bp_gamma_peak_corr->Add(Tin101_summed_bp_gamma_peak[0],1);
	Tin101_summed_bp_gamma_peak_corr->Add(Tin101_summed_bp_gamma_peak[1],-1);
	Tin101_summed_bp_gamma_rest_corr->Add(Tin101_summed_bp_gamma_rest[0],1);
	Tin101_summed_bp_gamma_rest_corr->Add(Tin101_summed_bp_gamma_rest[1],-1);

	Tin101_bp_gamma_peak_corr->Write();
	Tin101_summed_bp_gamma_peak_corr->Write();
	Tin101_bp_gamma_rest_corr->Write();
	Tin101_summed_bp_gamma_rest_corr->Write();

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

	/*Ag96_single_vs_summed->Write();

	Ag96_gamma_gamma->Write();

	Ag96_sum_E1E2_diff_470_740->Write();

	Ag96_sum_E1E2_diff_740_1249->Write();

	Ag96_EdT_470keVgammaGated->Write();

	Ag96_EdT_743keVgammaGated->Write();

	Ag96_EdT_1249keVgammaGated->Write();

	Ag96_EdT_all3Peaks_gammaGated->Write();

	Ag96_EdT_all3Peaks_Random_gammaGated->Write();

	Ag96_GammaT_betaT_all3Peaks->Write();

	Ag96_E_randomGamma->Write();

	Ag96_E_correlatedGamma->Write();

	Ag96_EdT_2461keVgammaGated->Write();

	Ag96_EdT_summed_gammaGated->Write();*/

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
				//IsoDir->Append(delayed1pEnergy[i][z].at(k));
				IsoDir->Append(delayed1pEnergyX[i][z].at(k));
				//IsoDir->Append(delayed1pEnergyY[i][z].at(k));
				IsoDir->Append(delayed1pEnergyRandom[i][z].at(k));
				IsoDir->Append(delayed1pEnergyAll[i][z].at(k));
				IsoDir->Append(implantBeta[i][z].at(k));
				IsoDir->Append(implantBeta1p[i][z].at(k));
				IsoDir->Append(implantE[i][z].at(k));
				IsoDir->Append(implantVelocityimplantE[i][z].at(k));
				IsoDir->Append(EdT[i][z].at(k));
				//IsoDir->Append(ExEyDiff[i][z].at(k));
				IsoDir->Append(implantVelocityAOQ[i][z].at(k));
				IsoDir->Append(implantEnergyAOQ[i][z].at(k));
			
			}
			
			//combined DSSD
			IsoDir->Append(decayEnergyAll[i].at(k));			
			IsoDir->Append(delayed1pEnergy_AllDSSD[i].at(k));						
			IsoDir->Append(delayed1pEnergyRandom_AllDSSD[i].at(k));
			IsoDir->Append(delayed1pEnergyAll_AllDSSD[i].at(k));
			IsoDir->Append(implantBetaAll[i].at(k));
			IsoDir->Append(implantBeta1pAll[i].at(k));
			IsoDir->Append(implantEAll[i].at(k));
			IsoDir->Append(implantVelocityimplantEAll[i].at(k));
			IsoDir->Append(EdTAll_ms[i].at(k));
			IsoDir->Append(EdTAll_us[i].at(k));
			IsoDir->Append(EdTAll_NoMultiGate[i].at(k));
			IsoDir->Append(EdTAll11[i].at(k));
			IsoDir->Append(EdTAll12[i].at(k));
			IsoDir->Append(EdTAll21[i].at(k));
			IsoDir->Append(EdTAll22[i].at(k));
			IsoDir->Append(implantVelocityAOQ_AllDSSD[i].at(k));
			IsoDir->Append(implantZ[i].at(k));
			IsoDir->Append(implantEnergyAOQ_AllDSSD[i].at(k));

			/*IsoDir->Append(ExEyAll[i].at(k));
			//IsoDir->Append(ExEyAll_gammaloop[i].at(k));
			IsoDir->Append(EdT_gammagate[i].at(k));
			IsoDir->Append(EdT_gammagate_longer[i].at(k));
			IsoDir->Append(EDiff_gammagate[i].at(k));
			IsoDir->Append(ExEy_gammagate[i].at(k));
			IsoDir->Append(NxNy_gammagate[i].at(k));
			IsoDir->Append(clustersize_gammagate[i].at(k));*/

			

			//gamma spectra correction
			for(int g = 0; g < 4; g++){
				if (g == 0 || g == 3){
					beta_gamma_EdT_s_corr[i].at(k)->Add(beta_gamma_EdT_s[i][g].at(k),1);
					beta_gamma_EdT_ms_corr[i].at(k)->Add(beta_gamma_EdT_ms[i][g].at(k),1);
					beta_gamma_EdT_us_corr[i].at(k)->Add(beta_gamma_EdT_us[i][g].at(k),1);
					//summed_beta_gamma_EdT_s_corr[i].at(k)->Add(summed_beta_gamma_EdT_s[i][g].at(k),1);
					//summed_beta_gamma_EdT_ms_corr[i].at(k)->Add(summed_beta_gamma_EdT_ms[i][g].at(k),1);
					//summed_beta_gamma_EdT_us_corr[i].at(k)->Add(summed_beta_gamma_EdT_us[i][g].at(k),1);
					bp_gamma_EdT_s_corr[i].at(k)->Add(bp_gamma_EdT_s[i][g].at(k),1);
					bp_gamma_EdT_ms_corr[i].at(k)->Add(bp_gamma_EdT_ms[i][g].at(k),1);
					bp_gamma_EdT_us_corr[i].at(k)->Add(bp_gamma_EdT_us[i][g].at(k),1);
					//summed_bp_gamma_EdT_s_corr[i].at(k)->Add(summed_bp_gamma_EdT_s[i][g].at(k),1);
					//summed_bp_gamma_EdT_ms_corr[i].at(k)->Add(summed_bp_gamma_EdT_ms[i][g].at(k),1);
					//summed_bp_gamma_EdT_us_corr[i].at(k)->Add(summed_bp_gamma_EdT_us[i][g].at(k),1);
				}

				if (g == 1 || g == 2){
					beta_gamma_EdT_s_corr[i].at(k)->Add(beta_gamma_EdT_s[i][g].at(k),-1);
					beta_gamma_EdT_ms_corr[i].at(k)->Add(beta_gamma_EdT_ms[i][g].at(k),-1);
					beta_gamma_EdT_us_corr[i].at(k)->Add(beta_gamma_EdT_us[i][g].at(k),-1);
					//summed_beta_gamma_EdT_s_corr[i].at(k)->Add(summed_beta_gamma_EdT_s[i][g].at(k),-1);
					//summed_beta_gamma_EdT_ms_corr[i].at(k)->Add(summed_beta_gamma_EdT_ms[i][g].at(k),-1);
					//summed_beta_gamma_EdT_us_corr[i].at(k)->Add(summed_beta_gamma_EdT_us[i][g].at(k),-1);
					bp_gamma_EdT_s_corr[i].at(k)->Add(bp_gamma_EdT_s[i][g].at(k),-1);
					bp_gamma_EdT_ms_corr[i].at(k)->Add(bp_gamma_EdT_ms[i][g].at(k),-1);
					bp_gamma_EdT_us_corr[i].at(k)->Add(bp_gamma_EdT_us[i][g].at(k),-1);
					//summed_bp_gamma_EdT_s_corr[i].at(k)->Add(summed_bp_gamma_EdT_s[i][g].at(k),-1);
					//summed_bp_gamma_EdT_ms_corr[i].at(k)->Add(summed_bp_gamma_EdT_ms[i][g].at(k),-1);
					//summed_bp_gamma_EdT_us_corr[i].at(k)->Add(summed_bp_gamma_EdT_us[i][g].at(k),-1);
				}

				//IsoDir->Append(summed_beta_gamma_EdT_s[i][g].at(k));
				//IsoDir->Append(summed_bp_gamma_EdT_s[i][g].at(k));

			}
			summed_beta_gamma_EdT_s_corr[i].at(k)->Add(summed_beta_gamma_EdT_s[i][0].at(k),1);
			summed_beta_gamma_EdT_ms_corr[i].at(k)->Add(summed_beta_gamma_EdT_ms[i][0].at(k),1);
			summed_beta_gamma_EdT_us_corr[i].at(k)->Add(summed_beta_gamma_EdT_us[i][0].at(k),1);

			summed_beta_gamma_EdT_s_corr[i].at(k)->Add(summed_beta_gamma_EdT_s[i][1].at(k),-1);
			summed_beta_gamma_EdT_ms_corr[i].at(k)->Add(summed_beta_gamma_EdT_ms[i][1].at(k),-1);
			summed_beta_gamma_EdT_us_corr[i].at(k)->Add(summed_beta_gamma_EdT_us[i][1].at(k),-1);

			summed_bp_gamma_EdT_s_corr[i].at(k)->Add(summed_bp_gamma_EdT_s[i][0].at(k),1);
			summed_bp_gamma_EdT_ms_corr[i].at(k)->Add(summed_bp_gamma_EdT_ms[i][0].at(k),1);
			summed_bp_gamma_EdT_us_corr[i].at(k)->Add(summed_bp_gamma_EdT_us[i][0].at(k),1);

			summed_bp_gamma_EdT_s_corr[i].at(k)->Add(summed_bp_gamma_EdT_s[i][1].at(k),-1);
			summed_bp_gamma_EdT_ms_corr[i].at(k)->Add(summed_bp_gamma_EdT_ms[i][1].at(k),-1);
			summed_bp_gamma_EdT_us_corr[i].at(k)->Add(summed_bp_gamma_EdT_us[i][1].at(k),-1);

			EdTAll_NoMultiGate_corr[i][0].at(k)->Add(EdTAll_NoMultiGate_corr[i][1].at(k),-1);
			EdTAll_ms_corr[i][0].at(k)->Add(EdTAll_ms_corr[i][1].at(k),-1);
			EdTAll_us_corr[i][0].at(k)->Add(EdTAll_us_corr[i][1].at(k),-1);
			//summed_p_gamma_E_p_E[i][0].at(k)->Add(summed_p_gamma_E_p_E[i][1].at(k),-1);
			//summed_beta_gamma_E_beta_E[i][0].at(k)->Add(summed_beta_gamma_E_beta_E[i][1].at(k),-1);

			//corrected
			IsoDir->Append(EdTAll_NoMultiGate_corr[i][0].at(k));
			IsoDir->Append(EdTAll_ms_corr[i][0].at(k));
			IsoDir->Append(EdTAll_us_corr[i][0].at(k));


			IsoDir->Append(beta_gamma_EdT_s_corr[i].at(k));
			IsoDir->Append(beta_gamma_EdT_ms_corr[i].at(k));
			IsoDir->Append(beta_gamma_EdT_us_corr[i].at(k));
			IsoDir->Append(summed_beta_gamma_EdT_s_corr[i].at(k));
			IsoDir->Append(summed_beta_gamma_EdT_ms_corr[i].at(k));
			IsoDir->Append(summed_beta_gamma_EdT_us_corr[i].at(k));
			IsoDir->Append(bp_gamma_EdT_s_corr[i].at(k));
			IsoDir->Append(bp_gamma_EdT_ms_corr[i].at(k));
			IsoDir->Append(bp_gamma_EdT_us_corr[i].at(k));
			IsoDir->Append(summed_bp_gamma_EdT_s_corr[i].at(k));
			IsoDir->Append(summed_bp_gamma_EdT_ms_corr[i].at(k));
			IsoDir->Append(summed_bp_gamma_EdT_us_corr[i].at(k));

			


			//bgrnd components
			for (int j=0; j<4; j++){
				//IsoDir->Append(beta_gamma_EdT_us[i][j].at(k));
				//IsoDir->Append(beta_gamma_EdT_ms[i][j].at(k));
				//IsoDir->Append(beta_gamma_EdT_s[i][j].at(k));

				if (j < 1){
					IsoDir->Append(summed_beta_gamma_EdT_us[i][j].at(k));
					IsoDir->Append(summed_beta_gamma_EdT_ms[i][j].at(k));
					IsoDir->Append(summed_beta_gamma_EdT_s[i][j].at(k));
				}
			}
			

			//IsoDir->Append(summed_p_gamma_E_p_E[i][0].at(k));
			//IsoDir->Append(summed_beta_gamma_E_beta_E[i][0].at(k));

			
			//IsoDir->Append(summed_beta_gamma_EdT_us[i][0].at(k));
			//IsoDir->Append(summed_beta_gamma_EdT_us[i][1].at(k));


			//IsoDir->Append(implantVelocityimplantZ[i].at(k));
			
		}//isotope_loop?
	}

	//std::cout<<"seg fault after here" <<std::endl;
	ofile->Write();
	//ofile->Close();

}//end of program