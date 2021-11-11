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
				if (abs((*beta).Ex-(*beta).Ey)<120){
					multix = (*beta).TFast & 0xFF;
					multiy = ((*beta).TFast >> 8) & 0xFF;
					for ( auto imp:(*beta).vectorOfImp ){ //if non-element gated histos needed, do here
						if ((*beta).z == (imp).Z){
							gammaVeto = true;
							//In97gammaVeto = false;
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
												decayEnergy[i][DSSD].at(j)->Fill((*beta).E);
												//ExEyDiff[i][DSSD].at(j)->Fill((*beta).Ex - (*beta).Ey);

												if (multix >= 0 && multiy >= 0){
													EdT[i][DSSD].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9, (*beta).E);
													EdT_ms[i][DSSD].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e6, (*beta).E);
													EdT_us[i][DSSD].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).E);
													
													if ((*beta).Ex>1500 && (*beta).Ey>1500){
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
												isotopeSumEnergy->Fill((*beta).E);
												decayEnergyAll[i].at(j)->Fill((*beta).E);
												if (multix >= 0 && multiy >= 0){ //beta-delayed protons
													if ((*beta).Ex>1500 && (*beta).Ey>1500){
														if (((*beta).T-(imp).TIME > 0)){
															delayed1pEnergy_AllDSSD[i].at(j)->Fill((*beta).E);
														}
														PID->Fill((imp).AOQ, (imp).ZET);
														//beta-p gamma loop
														//ProtonGammaSumTemp = 0;
														//ProtonGammaSumTempBg = 0;
														/*Sn101Counter = 0;
														Sn101CounterBg = 0;
														Sn101Counter_Pk = 0;
														Sn101CounterBg_Pk = 0;
														Ag94_Peak_Counter = 0;
														Ag94_Peak_CounterBg = 0;*/
														

														implantBeta1pAll[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9);

														if ((*beta).T-(imp).TIME < 0){
															delayed1pEnergyRandom_AllDSSD[i].at(j)->Fill((*beta).E);
														}

														delayed1pEnergyAll_AllDSSD[i].at(j)->Fill((*beta).E);
													}//end of lower beta-p energy cut
												}//end of beta-p multiplicity cut

												if ((*beta).nx == 1 && (*beta).ny == 1){
													EdTAll11[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9, (*beta).Ex);
													ExEy11[i].at(j)->Fill((*beta).Ex, (*beta).Ey);
												}
												if ((*beta).nx == 1 && (*beta).ny == 2){
													EdTAll12[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9, (*beta).Ex);
													ExEy12[i].at(j)->Fill((*beta).Ex, (*beta).Ey);
												}
												if ((*beta).nx == 2 && (*beta).ny == 1){
													EdTAll21[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9, (*beta).Ex);
													ExEy21[i].at(j)->Fill((*beta).Ex, (*beta).Ey);
												}
												if ((*beta).nx == 2 && (*beta).ny == 2){
													EdTAll22[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9, (*beta).Ex);
													ExEy22[i].at(j)->Fill((*beta).Ex, (*beta).Ey);
												}		
												EdTAll_NoMultiGate[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9, (*beta).Ex);
												EdTAll_us[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e3, (*beta).Ex);
												EdTAll_ms[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e6, (*beta).Ex);

											
												//ExEyAll[i].at(j)->Fill((*beta).Ex, (*beta).Ey);

												//DTAS 511/1022 coincidence check
												gammaVeto = true;
												//gamma gated spectra
												

												//beta - DTAS correlations
												

												if ((*beta).Ex<1500 && (*beta).Ey<1500){
													if ((*beta).nx<4 && (*beta).ny<4){
														implantBetaAll[i].at(j)->Fill(((*beta).T-(imp).TIME)/1.0e9);
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

		if ((*implant).T){
			for ( auto pid:(*implant).vectorOfPid ){
				PID_implant->Fill((*implant).aoq, (*implant).zet);
				int iDSSD = (*implant).z;
				for (int i =0; i<numElements; i++){
					for (int j = 0; j <= isotopeEnd[i]-isotopeStart[i]; j++){	
						if (particleCuts[i][j]->IsInside((*implant).aoq, (*implant).zet)){
							implantZ[i].at(j)->Fill((*implant).z);
							implantVelocityimplantZ[i].at(j)->Fill((pid).VELOCITY, (*implant).z);
							//implantEAll[i].at(j)->Fill((*implant).E);
							implantVelocityAOQ_AllDSSD[i].at(j)->Fill((*implant).aoq,(pid).VELOCITY);
							implantEnergyAOQ_AllDSSD[i].at(j)->Fill((*implant).aoq,(*implant).E);
							//implantE[i][iDSSD].at(j)->Fill((*implant).E);
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

	PID->Write();

	globalEnergy->Write();

	isotopeSumEnergy->Write();

	//Indium97_gammaveto_EdT->Write();

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
				IsoDir->Append(implantVelocityimplantE[i][z].at(k));
				IsoDir->Append(EdT[i][z].at(k));
				IsoDir->Append(EdT_ms[i][z].at(k));
				IsoDir->Append(EdT_us[i][z].at(k));
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
			IsoDir->Append(implantVelocityimplantEAll[i].at(k));
			IsoDir->Append(EdTAll_ms[i].at(k));
			IsoDir->Append(EdTAll_us[i].at(k));
			IsoDir->Append(EdTAll_NoMultiGate[i].at(k));
			IsoDir->Append(EdTAll11[i].at(k));
			IsoDir->Append(EdTAll12[i].at(k));
			IsoDir->Append(EdTAll21[i].at(k));
			IsoDir->Append(EdTAll22[i].at(k));

			IsoDir->Append(ExEy11[i].at(k));
			IsoDir->Append(ExEy12[i].at(k));
			IsoDir->Append(ExEy21[i].at(k));
			IsoDir->Append(ExEy22[i].at(k));

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
			/*for(int g = 0; g < 4; g++){
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

			}*/

			
			//summed_p_gamma_E_p_E[i][0].at(k)->Add(summed_p_gamma_E_p_E[i][1].at(k),-1);
			//summed_beta_gamma_E_beta_E[i][0].at(k)->Add(summed_beta_gamma_E_beta_E[i][1].at(k),-1);

			//corrected


			//IsoDir->Append(beta_gamma_EdT_s_corr[i].at(k));
			//IsoDir->Append(beta_gamma_EdT_ms_corr[i].at(k));
			//IsoDir->Append(beta_gamma_EdT_us_corr[i].at(k));
			//IsoDir->Append(summed_beta_gamma_EdT_s_corr[i].at(k));
			//IsoDir->Append(summed_beta_gamma_EdT_ms_corr[i].at(k));
			//IsoDir->Append(summed_beta_gamma_EdT_us_corr[i].at(k));
			//IsoDir->Append(bp_gamma_EdT_s_corr[i].at(k));
			//IsoDir->Append(bp_gamma_EdT_ms_corr[i].at(k));
			//IsoDir->Append(bp_gamma_EdT_us_corr[i].at(k));
			//IsoDir->Append(summed_bp_gamma_EdT_s_corr[i].at(k));
			//IsoDir->Append(summed_bp_gamma_EdT_ms_corr[i].at(k));
			//IsoDir->Append(summed_bp_gamma_EdT_us_corr[i].at(k));

			


			//bgrnd components
			/*for (int j=0; j<4; j++){
				//IsoDir->Append(beta_gamma_EdT_us[i][j].at(k));
				//IsoDir->Append(beta_gamma_EdT_ms[i][j].at(k));
				//IsoDir->Append(beta_gamma_EdT_s[i][j].at(k));

				if (j < 2){
					IsoDir->Append(summed_beta_gamma_EdT_us[i][j].at(k));
					IsoDir->Append(summed_beta_gamma_EdT_ms[i][j].at(k));
					IsoDir->Append(summed_beta_gamma_EdT_s[i][j].at(k));
					IsoDir->Append(summed_bp_gamma_EdT_us[i][j].at(k));
					IsoDir->Append(summed_bp_gamma_EdT_ms[i][j].at(k));
					IsoDir->Append(summed_bp_gamma_EdT_s[i][j].at(k));
				}
			}*/
			

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