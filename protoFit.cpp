#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TH1.h"
#include "TList.h"
#include "TMath.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include "TCanvas.h"
#include "TStyle.h"

#include <string>
#include <iostream>
#include <map>
#include <fstream>
#include <functional>

#include "ParticleCutsSn100.cxx"
#include "FitIsotopes.cpp"

//
double BatemanEquation(double *x, double *par, int histogramType, std::vector<std::list<isotope>> decayChain, std::vector<isotopeParameters> decayChainParameters){

    double returnValue;
	double bckgSlope = par[0];
	double bckgIntercept = par[1];
	double n0 = par[2];
	double betaEff = par[3];
	double halflife = par[4];
	double halflifeErr = par[5];
    double bphalflife = par[6];
    double bphalflifeErr = par[7];
	double bp = par[8];
	double bpErr = par[9];
	double betaEffRel = par[10];
	double protonEff = par[11];
	double aN1 = 0;
	double chainA = 0;
	double t = x[0];

    //define iterators here
    std::list<isotope>::iterator isotopeIt;
	double prod1 = 1;
	double prod2 = 1;
	double sum1 = 0;
	double sum1top = 0;
	double posInPath = 0;  
	int previousPath = -1;

    int currentIsoPos = -1;	
	std::vector<isotopeParameters>::iterator currentIsoPar;
	std::list<isotope>::iterator isotopeSumi;
	std::vector<isotopeParameters>::iterator isotopeSumiPar;
	std::list<isotope>::iterator isotopeProdj; 
	std::vector<isotopeParameters>::iterator isotopeProdjPar;
	std::list<isotope>::iterator isotopeProdk;
	std::vector<isotopeParameters>::iterator isotopeProdkPar;

	std::list<isotope>::iterator currentPathCheckIt;
	std::list<isotope>::iterator previousPathCheckIt;
	bool pathContributionStatus = true;

	returnValue = 0;
    //linear background
    if(t<0){
		returnValue = bckgSlope*t+bckgIntercept;
	}
    //beta vs bp decay curves
	if(t>0){
		returnValue = bckgSlope*(-1.0)*t+bckgIntercept;
		if(histogramType == 0){
			aN1 = betaEff * (TMath::Log(2)/halflife) * n0 * TMath::Exp(-1.0*(TMath::Log(2)/halflife)*t);
		}
		if(histogramType == 1){
			aN1 = protonEff * bp * (TMath::Log(2)/bphalflife) * n0 * TMath::Exp(-1.0*(TMath::Log(2)/bphalflife)*t);
		}
    }
    
    for(int path = 0; path < decayChain.size(); path++){
		//Loop over the different decay paths
		prod1 = 1;
		for(isotopeIt = decayChain.at(path).begin(); isotopeIt != decayChain.at(path).end(); isotopeIt++ ){
			//Loop through the decay chain
			currentIsoPar = GetIsotopeParamPos(&decayChainParameters, isotopeIt->name);
			/*if(isotopeIt == decayChain.at(path).begin() ){
					//If first skip over calculating N as it is a parameter in the first part
					//Does calculate the first product now as it will be used in the next iteration before being calculated again
					if(isotopeIt->decayType == 0){
						prod1 = (1.0-par[currentIsoPar->pn] - par[currentIsoPar->p2n]) * par[currentIsoPar->lambda];
					}
					if(isotopeIt->decayType == 1){
						prod1 = par[currentIsoPar->pn]* par[currentIsoPar->lambda];
					}
					if(isotopeIt->decayType == 2){
						prod1 = par[currentIsoPar->p2n]*par[currentIsoPar->lambda];
					}
					continue;
				}//End if on first item in the decay chain*/
			//Calculate the sum for N_k(t)
			//Sums exp(-lambda_i*t)/prod(lamda_j-lamda_i)
			prod1 = 1.0;
			if(isotopeIt != decayChain.at(path).begin()){
				for(isotopeProdk = decayChain.at(path).begin(); isotopeProdk != isotopeIt; isotopeProdk++ ){	
					isotopeProdkPar = GetIsotopeParamPos(&decayChainParameters, isotopeProdk->name);	
					if(isotopeProdk->decayType == 0){
						prod1 = prod1 * (1.0-par[isotopeProdkPar->bp]) * (TMath::Log(2)/par[isotopeProdkPar->halflife]);
					}
					if(isotopeProdk->decayType == 1){
						prod1 = prod1 * par[isotopeProdkPar->bp]* (TMath::Log(2)/par[isotopeProdkPar->bphalflife]);
					}
				}
			}
			sum1 = 0;
			sum1top = 0;
			for(isotopeSumi = decayChain.at(path).begin();isotopeSumi != std::next(isotopeIt); isotopeSumi++ ){
				isotopeSumiPar = GetIsotopeParamPos(&decayChainParameters, isotopeSumi->name);
				sum1top = TMath::Exp(-1.0 * (TMath::Log(2)/par[isotopeSumiPar->halflife]) * t);
				prod2 = 1;
				for(isotopeProdj = decayChain.at(path).begin(); isotopeProdj != std::next(isotopeIt); isotopeProdj++ ){
					if(isotopeProdj->name == isotopeSumi->name){
						continue;
					}
					isotopeProdjPar = GetIsotopeParamPos(&decayChainParameters, isotopeProdj->name);
					prod2 = prod2*((TMath::Log(2)/par[isotopeProdjPar->halflife])- (TMath::Log(2)/par[isotopeSumiPar->halflife]));
				}
				sum1 = sum1 + sum1top/prod2;
			}

				//Check if the current part of the path has already been included.
			if(path == 0 || (path == previousPath && pathContributionStatus)){
				pathContributionStatus = true;
			}
			else{
				pathContributionStatus = true;
				for(int i = 0; i < path; i++){
					currentPathCheckIt = decayChain.at(path).begin();
					for(previousPathCheckIt = decayChain.at(i).begin(); previousPathCheckIt != decayChain.at(i).end();previousPathCheckIt++){	

						if(previousPathCheckIt->name == currentPathCheckIt->name){
							if(currentPathCheckIt == isotopeIt){
								//Current isotope already had contribution along this path accounted for
								//This isotope should not be accounted for
								pathContributionStatus = false;
								break;
								}
								if(currentPathCheckIt != isotopeIt){
									currentPathCheckIt++;
									if(currentPathCheckIt == decayChain.at(path).end()){
										continue;
									}
								}
							}
						else if(previousPathCheckIt->name != currentPathCheckIt->name){
							break;
						}
					}
					if (pathContributionStatus == false){
						break;
					}
				}
				previousPath = path;

				if(isotopeIt != decayChain.at(path).begin() && histogramType == 0 && pathContributionStatus){
					chainA = chainA + ((TMath::Log(2)/par[currentIsoPar->halflife]) * par[currentIsoPar->betaEffRel]* betaEff * n0 * sum1 * prod1);
				}
				else if(isotopeIt != decayChain.at(path).begin() && histogramType == 1 && isotopeIt->decayType == 1 && pathContributionStatus){
					chainA = chainA + ((TMath::Log(2)/par[currentIsoPar->bphalflife]) * par[currentIsoPar->betaEffRel] * protonEff * par[currentIsoPar->bp] * n0 * sum1 * prod1);
				}
				/*if(isotopeIt->decayType == 0){
					prod1 = prod1*(1.0-par[currentIsoPar->pn] - par[currentIsoPar->p2n]) * par[currentIsoPar->lambda];
				}
				if(isotopeIt->decayType == 1){
					prod1 = prod1*par[currentIsoPar->pn]* par[currentIsoPar->lambda];
				}
				if(isotopeIt->decayType == 2){
					prod1 = prod1*par[currentIsoPar->p2n]*par[currentIsoPar->lambda];
				}*/

			}//End loop over chain
		}//End loop over paths

		returnValue = returnValue + aN1 + chainA;

    }// end if on positive time    
    
	return returnValue;

}

int decayFit(std::string iName, std::string isotopeList){
	
	std::vector< isotope > isotopes;
	std::vector< std::list<isotope>> decayChain;
	std::vector< isotopeParameters > decayChainParameters;
	std::vector< std::list<isotope>> * decayChainPoint = &decayChain;
	std::vector< isotopeParameters > * decayChainParametersPoint = &decayChainParameters;
	int decayChainStatus = 0;

	TFile * ifile = TFile::Open( iName.c_str(), "read");
	TH1F *fitHistogram;	

	fitHistogram = (TH1F*) ifile->Get("Sn101ImplantBeta");//histogram to be fitted
	fitHistogram->Rebin(10);


	CreateIsotopeList(&isotopes, isotopeList);

	std::cout << isotopes.size() <<std::endl;
	std::cout << isotopes.front().name <<std::endl;

	decayChainStatus = DecayChain(&decayChain, &isotopes, 50, 101, -1);//z, a of histogram nucleus
	if(decayChainStatus == -1){
		return -1;
	}
	PrintDecayChain(&decayChain);
	DeclareParameters(&decayChain, &decayChainParameters);
	ReadParameters(&decayChainParameters);

	auto BatemanFunc = std::bind(BatemanEquation,std::placeholders::_1,std::placeholders::_2,0,decayChain,decayChainParameters); //partial application of bateman function with equation, decay chain and parameters. T and par as placeholder parameters
	std::cout << "Decay chain size " << decayChain.size() << std::endl;

	TF1 * background = new TF1("background","[0]*x+[1]",-9.75,-0.1);
	TF1 * genBateman = new TF1("genBateman",BatemanFunc,0.1,10,decayChainParameters.back().protonEff); //this line causes compilation error

	
	SetFitParameters(genBateman,&isotopes,&decayChainParameters);

	background->SetParameter(0,0);
	background->SetParameter(1,500);


	genBateman->SetParameter(2,530);
	genBateman->SetParLimits(2,10,1e9);
	genBateman->FixParameter(3,1.0);
	fitHistogram->GetXaxis()->SetRangeUser(-10,10);

	fitHistogram->Fit(background,"MRL");

	genBateman->FixParameter(0,background->GetParameter(0));
	genBateman->FixParameter(1,background->GetParameter(1));
	//genBateman->SetParLimits(4,0,2);

	fitHistogram->Draw();
	fitHistogram->Fit(genBateman,"MRL");
	background->Draw("SAME");
	genBateman->Draw("SAME");

	std::cout << genBateman->GetChisquare() << std::endl;
	std::cout << genBateman->GetNDF() << std::endl;
	



	return 0;
}