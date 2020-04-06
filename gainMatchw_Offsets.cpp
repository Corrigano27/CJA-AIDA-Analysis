#include "TFile.h"
#include "TTree.h"
#include "TList.h"
#include "TChain.h"
#include "TKey.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TString.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include "TF1.h"
#include "TMath.h"
#include "TMinuit.h"
#include "Math/PdfFuncMathCore.h"
#include <limits>

using namespace std;

TString ToString(int num) {
	ostringstream start;
	start<<num;
	TString start1=start.str();
	return start1;
}

//define first and last p/n strips respectively
const int firstPstrip = 0;
const int lastPstrip = 128;

const int firstNstrip = 0;
const int lastNstrip = 128;

const int pArraysize = lastPstrip - firstPstrip;
const int nArraysize = lastNstrip - firstNstrip;

const int totalStrips = lastPstrip + lastNstrip;

//global array of amplitudes in each single pixel plot produced in singleChannelExEy
vector<Double_t>* AxVecPtr = new vector<Double_t>();
vector<Double_t>* AyVecPtr = new vector<Double_t>();

//Likelihood function to be maximised for each single pixel plot
//chose Lorentzian distribution of natural logs of amplitudes here
Double_t Likelihood(Double_t* AxPtr, Double_t* AyPtr, Double_t par[])
{
	Double_t Ax = *AxPtr;
	Double_t Ay = *AyPtr; 
	Double_t width = par[0];
	Double_t Spnmean = par[1];
	Double_t Opnmean = par[2];
	//cout << Ax << " " << Ay << endl;
	//Double_t logRatio = TMath::Log(Ax/Ay);
	//Double_t logSpn = TMath::Log(Spnmean);
	//Double_t normconst = 1/(width*TMath::Sqrt(2)*TMath::Pi());
	//Double_t R = (Spnmean*Ay-Ax)/width;
	Double_t L = 0;
	if ( Ax > 0 && Ay > 0 && Spnmean > 0 && width > 0){
		L = (1.0/TMath::Pi())*(width/(TMath::Power(TMath::Log((Ax - Opnmean)/Ay)-TMath::Log(Spnmean),2)+TMath::Power(width,2)));
	}
	return L;
}
/*-----------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------*/
//function wrapped in TMinuit to be minimised.
//To maximise the Likelihood defined above, we minimes the negative log of the Likelihood, which is computationally easier.
void minfcn(Int_t &npar, Double_t *deriv, Double_t &f, Double_t par[], Int_t flag)
{
	vector<Double_t> AxVec = *AxVecPtr;
	vector<Double_t> AyVec = *AyVecPtr;
	Int_t arrsize = AxVec.size();

//calculate negative log likelihood
	Double_t LnL = 0.0;
	for (Int_t i=0; i<arrsize; i++){
		Double_t Ax = AxVec[i];
		Double_t Ay = AyVec[i];
		//cout << Ax << " " << Ay << endl;
		Double_t LH = Likelihood(&Ax,&Ay, par);
		if (LH > 0.0){
			LnL += TMath::Log(LH);
		}
		else {
			cout << "Error! Negative Likelihood Function!" <<endl;
		}
	}
	f = -LnL;
}
/*------------------------------------------------------------------------------------------------------------------*/

//global functions for calculating the gain

//global arrays storing parameters returned from first minimisation to be used in second
Double_t widthSpn[pArraysize][nArraysize];

Double_t Spn[pArraysize][nArraysize];

Double_t deltaSpn[pArraysize][nArraysize];

Double_t Opn[pArraysize][nArraysize];

Double_t deltaOpn[pArraysize][nArraysize];
/*
Double_t widthSpn[pArraysize][nArraysize];
*/

//chi squared minimisation, wrapped in TMinuit
void fcn(Int_t &npar, Double_t *deriv, Double_t &f, Double_t par[], Int_t flag)
 {
 
 //calculate chisquare
 	//Double_t nGain = par[n];
 	//Double_t pGain = par[p];

    Double_t chisq = 0;
    Double_t slopeTerms;
    Double_t offsetTerms;

    //cout << "array sizes p: " << pArraysize << " n: " << nArraysize << endl; 
    for (int p=0; p<pArraysize; p++) {
    	for (int n=0; n<nArraysize; n++){
    		//cout << "spn at p = " <<p + firstPstrip<<", n = " <<n + firstNstrip<< ": " << Spn[p][n] << endl;
    		if (Spn[p][n] != 0){
      			slopeTerms  = (Spn[p][n]-(par[n]/par[p+lastNstrip]))/deltaSpn[p][n];

      		}
      		else {
      			slopeTerms = 0;
      		}

      		if (Opn[p][n] != 0){
      			offsetTerms = (Opn[p][n]-((par[n+(totalStrips)]-par[p+(totalStrips+lastNstrip)])/par[p+lastNstrip]))/deltaOpn[p][n]; //need to think about offset parameters
      		}

      		else {
      			offsetTerms = 0;
			}

      	chisq += slopeTerms*slopeTerms + offsetTerms*offsetTerms;
      	}
    }
    f = chisq;
 }


/*------------------------------------------------------------------------------------------------------------------*/
//first part takes data from TTree, then minimises log likelihood function to find slope factor Spn for every p and n pixel 
int gainMatch(std::string iName){


	//----------------------------------------------------------------------------------------------------------------
	//opening multiple files with TChain

	std::string rootFile;
	std::string oName;
	std::ifstream rootFiles;
	rootFiles.open( iName );
	if (!rootFiles.is_open()){
		std::cerr << " Input file not opened exiting program" << std::endl;
		return -1; 
	}
	else {
		std::cout << "File ''" << iName << "'' is open" << std::endl;
	}

	TChain chain("AIDA_hits");

	while ( std::getline ( rootFiles, rootFile )){

		chain.Add( rootFile.c_str() );
		std::cout << "Added " << rootFile.c_str() << " to the chain." << std::endl;

	}

 	size_t lastindex = iName.find_last_of(".");
	oName = iName.substr(0, lastindex);
	oName+="_OffsetsAndGains.txt";
  	//Open the tree and create the branch to write to
  	TFile * ofile = TFile::Open( oName.c_str());



	std::cout << "Input and output files open" << std::endl;

	std::cout << "Output file name " << oName << std::endl;




//old way

	//TFile * inputFile = TFile::Open(fileName.c_str(),"read");

	// *next 3 lines are new* 
	/*
	TTree *inputTTree = chain.GetTree();
	  if( !inputTTree ) {
	  	std::cout << "Problem opening tree" << std::endl;
	  }
	*/

	 /*
	TFile *fileOut = new TFile("140519_TEST.root","RECREATE");
	size_t lastindex = fileName.find_last_of(".");
 	std::string oName = fileName.substr(0, lastindex);
 	//oName+="_Gains.root";
	//TFile * ofile = TFile::Open( oName.c_str(), "recreate");

	oName = fileName.substr(0, lastindex);
	oName+= "_Offsets&OGains.dat";
	std::cout << "Output file name " << oName << std::endl;

	*/
	ofstream outf;
	outf.open(oName);
	//TTree* tree = new TTree("ExEy", "Ex and Ey values for each strip");

	//create histograms
	//vector of histograms(histogram number = (DSSD*128*128)+((y-1)*128)+x-1)
	//loop over(6*128*128)
	std::cout <<" Before histograms" << std::endl;

	std::vector<std::pair<double, double>> amplitudes[6*128*128];
	std::cout << "After histograms" << std::endl;
	//create histogram
	//push histogram back
	struct aidaData{
    aidaData():T(0),Tfast(0),E(0),EX(0),EY(0),x(0),y(0),z(0),nx(0),ny(0),nz(0),ID(8){}
    ULong64_t       T;
    ULong64_t       Tfast; 
    Double_t        E; 
    Double_t        EX; 
    Double_t        EY; 
    Double_t        x; 
    Double_t        y;
    Double_t        z; 
    Int_t           nx; 
    Int_t           ny; 
    Int_t           nz; 
    UChar_t         ID;
  };

  aidaData inputEntry;
  chain.SetBranchAddress("aida_hit", & (inputEntry.T) );
  Long64_t nEntries = chain.GetEntries();
  std::cout << "Nentries" << nEntries << std::endl;
  double dx; 
  double dy;
  std::pair<double, double> ampPair;

  for( Long64_t iEntry = 0; iEntry < nEntries; ++iEntry )
  {
	chain.GetEntry(iEntry);
    //std::cout << inputEntry.x << std::endl;
	//high energy event ID
    if(inputEntry.nx == 1 && inputEntry.ny==1 && inputEntry.ID==5 && inputEntry.EX>2500){ //set selection checks
    	ampPair.first = inputEntry.EX;
    	ampPair.second = inputEntry.EY;
    	amplitudes[(int)((inputEntry.z*128*128)+((inputEntry.y)*128)+inputEntry.x)].push_back(ampPair);
    }

  }

  //---------------------------------------------------------------------------------------------------------------------

	cout << "Finding Spn values" <<endl;

	float DetPGain[6][128];
	float DetNGain[6][128];

	float DetPOffset[6][128];
	float DetNOffset[6][128];

	for (int dssd = 0; dssd < 6; dssd ++){
		cout<<"Processing dssd: "<<dssd<<endl;
		for (int p=firstPstrip; p<lastPstrip; p++){
			TString numstr=ToString(p);
			for (int n=firstNstrip; n<lastNstrip; n++){
				if(AxVecPtr->size()>0){
					AxVecPtr->clear();
					AyVecPtr->clear();
				}
				for(unsigned int i = 0; i < amplitudes[(dssd*128*128)+(n*128)+p].size();i++){
					AxVecPtr->push_back(amplitudes[(dssd*128*128)+(n*128)+p].at(i).first);
					AyVecPtr->push_back(amplitudes[(dssd*128*128)+(n*128)+p].at(i).second);

				}
				//std::cout <<" Finished reading in pixel" << std::endl;

				TString numstr2=ToString(n);
				//TString histoname = "histo" + numstr + numstr2;
				TString sel = "pPixelNo==" + numstr + " && nPixelNo==" + numstr2;

				const int npar = 3;
				TMinuit minuit(npar);
				//minuit.SetPrintLevel(-1);
				minuit.SetFCN(minfcn);

				Double_t par[npar]; //start value
				Double_t stepSize[npar]; //step size
				Double_t minVal[npar]; //minimum bound on parameter
				Double_t maxVal[npar]; //maximum bound on parameter
				string parName[npar];

				//NEED TO THINK ABOUT EMPTY PIXELS FOR OFFSETS

				par[0] = 0.02; //guess
				stepSize[0] = 0.; //eg take 0.01 of start value
				minVal[0] = 0.02; // min and max = 0, parameter unbounded
				maxVal[0] = 0.02;
				parName[0] = "width";

				if (amplitudes[(int)((inputEntry.z*128*128)+((inputEntry.y)*128)+inputEntry.x)].size()==0){

					par[1] = 0.;
					stepSize[1] = 0.; //ed take 0.1 of start value
					minVal[1] = 0.; //min and max = 0, parameter unbounded
					maxVal[1] = 0.;
					parName[1] = "Spn";

					par[2] = 0.;
					stepSize[2] = 0.; //ed take 0.1 of start value
					minVal[2] = 0.; //min and max = 0, parameter unbounded
					maxVal[2] = 0.;
					parName[2] = "Opn";


				}

				else {
					par[1] = 1.0; //guess
					stepSize[1] = 0.001; //ed take 0.1 of start value
					minVal[1] = 0.4; //min and max = 0, parameter unbounded
					maxVal[1] = 1.6;
					parName[1] = "Spn";

					par[2] = 0.; //guess
					stepSize[2] = 0.01; //ed take 0.1 of start value
					minVal[2] = -500.0; //min and max = 0, parameter unbounded
					maxVal[2] = 500.0;
					parName[2] = "Opn";
				}

				//setup Parameters
				for (int k=0; k<npar; k++){
					minuit.DefineParameter(k, parName[k].c_str(), par[k], stepSize[k], minVal[k], maxVal[k]);
				}

				//cout<<"Runnning Migrad()......"<<endl;
				//perform minimization!
				minuit.Migrad();
				//cout<<"Migrad() completed......"<<endl;

				//Get Minuit results
				Double_t outpar[npar], err[npar];
	  			for (int l=0; l<npar; l++){
	    			minuit.GetParameter(l,outpar[l],err[l]);
	    			//cout<<"Fitted parameter: "<<l<<" is: "<<outpar[l]<<" +/- "<<err[l]<<endl;
	  			}

	  			int pcorr = p - firstPstrip;
	  			int ncorr = n - firstNstrip; //these make sure arrays behave

	  			Spn[pcorr][ncorr]=outpar[1];
	  			deltaSpn[pcorr][ncorr]=err[1];
	  			widthSpn[pcorr][ncorr]=outpar[0];
	  			Opn[pcorr][ncorr]=outpar[2];
	  			deltaOpn[pcorr][ncorr]=err[2];

	  			

			}
		}

		for (int q=0; q<pArraysize; q++){
			for (int h=0; h<nArraysize; h++){
				//cout << "x = " <<q + firstPstrip<< ", y = " <<h + firstNstrip<<", Spn = "<<Spn[q][h]<<" +/- "<<deltaSpn[q][h]<< " width of " << widthSpn[q][h] << endl;
			}
		}

	


		cout << "Found Spn values" <<endl;
		/*-----------------------------------------------------------------------------------------------------------------------------------------*/
		//minimisation of chi^2 to find gains

		cout << "Finding gain coefficients and offsets" <<endl;

		const int npar2 = 2*pArraysize + 2*nArraysize;
		TMinuit minuit2(npar2);
		//minuit2.SetPrintLevel(-1);
		minuit2.SetFCN(fcn);

		Double_t par2[npar2]; //start value
		Double_t stepSize2[npar2]; //step size
		Double_t minVal2[npar2]; //minimum bound on parameter
		Double_t maxVal2[npar2]; //maximum bound on parameter
		string parName2[npar2];
		
		for (int i = 0; i<nArraysize; i++){
			if (i != 6){
				par2[i] = 1.0; //guess
				stepSize2[i] = 0.001; //eg take 0.01 of start value
				minVal2[i] = 0.; // min and max = 0, parameter unbounded
				maxVal2[i] = 0.;
				parName2[i] = "sn";
			}
			else {
				par2[i] = 1.0; //guess
				stepSize2[i] = 0.; //eg take 0.01 of start value
				minVal2[i] = 1.0; // min and max = 0, parameter unbounded
				maxVal2[i] = 1.0;
				parName2[i] = "sn";
			}	
		}

		for (int j = nArraysize; j<(nArraysize+pArraysize); j++){
			par2[j] = 1.0; //guess
			stepSize2[j] = 0.001; //ed take 0.1 of start value
			minVal2[j] = 0.; //min and max = 0, parameter unbounded
			maxVal2[j] = 0.;
			parName2[j] = "sp";
		}

		for (int k = (nArraysize+pArraysize); k<(2*nArraysize+pArraysize); k++){
			if (k != 6){
				par2[k] = 0.; //guess
				stepSize2[k] = 0.01;
				minVal2[k] = 0.; //min and max = 0, parameter unbounded
				maxVal2[k] = 0.;
				parName2[k] = "on";
			}
			else {
				par2[k] = 0.;
				stepSize2[k] = 0.;
				minVal2[k] = 0.;
				maxVal2[k] = 0.;
				parName2[k] = "on";
			}
		}

		for (int l = (2*nArraysize+pArraysize); l<npar2; l++){
			par2[l] = 0.; //guess
			stepSize2[l] = 0.01;
			minVal2[l] = 0.; //min and max = 0, parameter unbounded
			maxVal2[l] = 0.;
			parName2[l] = "op";
		}

				//setup Parameters
		for (int k=0; k<npar2; k++){
			minuit2.DefineParameter(k, parName2[k].c_str(), par2[k], stepSize2[k], minVal2[k], maxVal2[k]);
		}

		cout<<"Runnning Migrad()......"<<endl;
		//perform minimization!
		minuit2.Migrad();
		cout<<"Migrad() completed......"<<endl;

		float nGain[nArraysize];
		float pGain[pArraysize];

		float nOffset[nArraysize];
		float pOffset[pArraysize];

		//Get Minuit results
		Double_t outpar2[npar2], err2[npar2];

	  	for (int l=0; l<nArraysize; l++){
	    	minuit2.GetParameter(l,outpar2[l],err2[l]);
	    	nGain[l] = outpar2[l];
	    	//cout<<"nGain "<<l + firstNstrip<<" is "<<nGain[l]<<" +/- "<<err2[l]<<endl;
	  	}

	  	for (int l=nArraysize; l<(nArraysize+pArraysize); l++){
	    	minuit2.GetParameter(l,outpar2[l],err2[l]);
	    	pGain[l-nArraysize] = outpar2[l];
	    	//cout<<"pGain "<<l-nArraysize + firstPstrip<<" is "<<pGain[l-nArraysize]<<" +/- "<<err2[l]<<endl;
	  	}

	  	for (int l = (nArraysize + pArraysize); l<(2*nArraysize+pArraysize); l++){
	    	minuit2.GetParameter(l,outpar2[l],err2[l]);
	    	nOffset[l-(nArraysize + pArraysize)] = outpar2[l];
	    	//cout<<"nGain "<<l + firstNstrip<<" is "<<nGain[l]<<" +/- "<<err2[l]<<endl;
	  	}

	  	for (int l=(2*nArraysize+pArraysize); l<npar2; l++){
	    	minuit2.GetParameter(l,outpar2[l],err2[l]);
	    	pOffset[l-(2*nArraysize+pArraysize)] = outpar2[l];
	    	//cout<<"pGain "<<l-nArraysize + firstPstrip<<" is "<<pGain[l-nArraysize]<<" +/- "<<err2[l]<<endl;
	  	}



	  	for (int p=0; p<pArraysize; p++){
			for (int n=0; n<nArraysize; n++){
				//cout << "p = " << p + firstPstrip << ", n = " << n + firstNstrip << ", Spn = " << Spn[p][n] << " +/- " << deltaSpn[p][n] << ", width = " << widthSpn[p][n] << ", sn/sp = " << nGain[n]/pGain[p] <<", deltaSpn = " << Spn[p][n] - nGain[n]/pGain[p] <<", sp = " << pGain[p] <<  ", sn = " << nGain[n] <<endl;
			}
			
		}

		//this bit could be generalised for any number of n or p strips - job for later

		for (int p=0; p<pArraysize; p++){
			DetPGain[dssd][p] = pGain[p];
			DetNGain[dssd][p] = nGain[p];
			DetPOffset[dssd][p] = pOffset[p];
			DetNOffset[dssd][p] = nOffset[p];
			if (dssd ==4){
			//cout << DetNGain[dssd][p] << " " << DetPGain[dssd][p] << endl;
			}
		}	

	  	cout << "Found gains" <<endl;



	} //end of dssd loop




		/*outputting gains in form that can be used by AIDASort*/
	//int feeDSSDMap[24];						//Which DSSD does a FEE correspond to
	//int feeSideMap[24];						//What side of the detector is a FEE64
	//int feeStripMap[24]; 					//How does the FEE Map to the DSSD (1:Left/Bottom or 2: Right/Top)

	int fee64, dssd, channelID;
	int side;
	double value;
	std::string line;


	const double feeChannelOrder[64]={62., 63., 59., 60., 61., 56., 57., 58., 52., 53., 54., 55., 49., 50., 51., 45.,
											46., 47., 48., 42., 43., 44., 38., 39., 40., 41., 35., 36., 37., 31., 32., 33.,
											34., 28., 29., 30., 24., 25., 26., 27., 21., 22., 23., 17., 18., 19., 20., 14.,
											15., 16., 10., 11., 12.,  7.,  3.,  0.,  8.,  4.,  1.,  9.,  5.,  2., 13.,  6.};

	int feeDSSDMap[24]={5,6,6,5,6,5,5,6,3,4,4,3,4,3,3,4,1,2,2,1,2,1,1,2};

	int feeSideMap[24]={1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0};

	int feeStripMap[24]={1,1,1,1,2,2,2,2,1,1,1,1,2,2,2,2,1,1,1,1,2,2,2,2};

	//absolute gain coefficients determined from alpha background. Each element in array is dssd.
    //old calib using pulser
	//float absGainX[6]={0.7255, 0.7331, 0.721, 0.7254, 0.7429, 0.7388};

	//float absGainY[6]={0.7275, 0.7388, 0.7258, 0.7254, 0.7429, 0.734};

	//new calib
	float absGainX[6] = {0.70698, 0.71615, 0.71250, 0.71188, 0.70516, 0.71617};

	float absGainY[6] = {0.71249, 0.7235, 0.7231, 0.73001, 0.71615, 0.7262}; 

	// uncomment these when performing alpha background absolute calibration
	//float absGainX[6]={1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

	//float absGainY[6]={1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

	float GainChannel[24][64];

	float OffsetChannel[24][64];


	for (int dssd = 0; dssd < 6; dssd++){
		for (int strip = 0; strip<128; strip++){
			//cout<<DetPGain[dssd][i]<<endl;
			if (strip < 64){
				//cout << "Top-Left" << endl;
				for (int feechannel = 0; feechannel<64; feechannel++){
					if (feeChannelOrder[feechannel] == strip){
						//cout << "Strip i = " << i << " corresponds to FEE channel " << feechannel << endl;
						for (int fee=0; fee<24; fee++){
							if (feeSideMap[fee] == 1 && feeStripMap[fee] == 1 && feeDSSDMap[fee] == dssd+1){
								GainChannel[fee][feechannel]=absGainX[dssd]*DetPGain[dssd][strip];
								OffsetChannel[fee][feechannel]=DetPOffset[dssd][strip];
							}
							else if (feeSideMap[fee] == 0 && feeStripMap[fee] == 1 && feeDSSDMap[fee] == dssd+1){
									GainChannel[fee][feechannel]=absGainY[dssd]*DetNGain[dssd][strip];
									OffsetChannel[fee][feechannel]=DetNOffset[dssd][strip];
							}
						}  
					}
				}
			}
			else {
				//cout << "Bottom-Right" << endl;
				for (int feechannel = 0; feechannel<64; feechannel++){
					if (feeChannelOrder[feechannel] == 127-strip){
						//cout << "Strip i = " << i << " corresponds to FEE channel " << feechannel << endl;
						for (int fee=0; fee<24; fee++){
							if (feeSideMap[fee] == 1 && feeStripMap[fee] == 2 && feeDSSDMap[fee] == dssd+1){
								GainChannel[fee][feechannel]=absGainX[dssd]*DetPGain[dssd][strip];
								OffsetChannel[fee][feechannel]=DetPOffset[dssd][strip];
							}
							else if (feeSideMap[fee] == 0 && feeStripMap[fee] == 2 && feeDSSDMap[fee] == dssd+1){
									GainChannel[fee][feechannel]=absGainY[dssd]*DetNGain[dssd][strip];
									OffsetChannel[fee][feechannel]=DetNOffset[dssd][strip];
							}
								
						}

					}
				}
			}
		}
	}

	for (int i = 0; i<24; i++){
		for (int j = 0; j<64; j++){
			outf << "adcGain" << '\t' << i+1 << '\t' << j << '\t' << GainChannel[i][j] << std::endl;
		}
	
	}

	outf << " " << std::endl;

	for (int i = 0; i<24; i++){
		for (int j = 0; j<64; j++){
			outf << "adcOffset" << '\t' << i+1 << '\t' << j << '\t' << OffsetChannel[i][j] << std::endl;
		}
	}


/*-----------------------------------------------------------------------------------------------------------------------------------------*/
  	
  	//script to produce test plot for TESTING PURPOSES



  	cout << "Writing to file" <<endl;
/*
  	float Ax2, Ay2;
	int pPixel2;
	int nPixel2;
	int detNo = 5;


	TH2D * exEy;
	exEy = new TH2D("ExEy","",500,0,20e3,500,0,20e3);
	TH2D * exEyPre = new TH2D("ExEyPre","",500,0,20e3,500,0,20e3);
	TH2D * diffStripX = new TH2D("StripDifX","",128,0,128,1200,-600,600);
	TH2D * diffStripY = new TH2D("StripDifY","",128,0,128,1200,-600,600);
	TH1D * exEyDiff = new TH1D("Ex-Ey","",3000,-1500,1500);
	TH1D * exEyDiffPre= new TH1D("Ex-EyPre","",3000,-1500,1500);

	for (int p=firstPstrip; p<lastPstrip; p++){
		TString numstr=ToString(p);
		for (int n=firstNstrip; n<lastNstrip; n++){
			TString numstr2=ToString(n);
			float yGain = DetNGain[0][n-firstNstrip];
			float xGain = DetPGain[0][p-firstPstrip];
			pPixel2=(int)p;
			nPixel2=(int)n;
			int pixel = (detNo*128*128)+(n*128)+p;
			if (amplitudes[pixel].size() ==0)continue;
			for(unsigned int i = 0; i < amplitudes[(detNo*128*128)+(n*128)+p].size();i++){
				Ax2 = amplitudes[(detNo*128*128)+(n*128)+p].at(i).first*xGain;
				Ay2 = amplitudes[(detNo*128*128)+(n*128)+p].at(i).second*yGain;
				diffStripX->Fill(p,Ax2-Ay2);
				diffStripY->Fill(n,Ax2-Ay2);
				exEy->Fill(Ax2,Ay2);
				exEyPre->Fill(Ax2*xGain,Ay2*yGain);
				exEyDiff->Fill(Ax2-Ay2);
				exEyDiffPre->Fill((Ax2/xGain)-(Ay2/yGain));
			}
			//cout<<"applying gainmatch pPixel: " << numstr << " nPixel: " << numstr2 << endl;
			//cout<<"pPixel number: " << numstr << " nPixel number: " << numstr2 << endl;
		}
	}
	TCanvas * c1 = new TCanvas();
	diffStripX->Draw("colz");
	TCanvas * c2 = new TCanvas();
	diffStripY->Draw("colz");
	TCanvas * c3 = new TCanvas();
	exEyDiff->Draw();
	exEyDiff->GetYaxis()->SetRangeUser(0, 450);
	TCanvas * c4 = new TCanvas();
	exEyDiffPre->Draw();
	exEyDiffPre->GetYaxis()->SetRangeUser(0, 450);
	TCanvas * c5 = new TCanvas();
	exEyPre->Draw();
	TCanvas * c6 = new TCanvas();
	exEy->Draw();
	*/

	//tree2->Write();
	//TH2D* h5 = new TH2D("h5", "h5", 400, 0, 3e3, 400, 0, 5e3);
	//tree2->Draw("AmpYG:AmpXG",,"colz");
	//tree->Draw("AmpY-AmpX", "", "same");
	//tree->SetLineColor(2);
	//h5->Draw("same");

	//AidaHits -> Fill("Ey:Ex>>ExEyhist");
	//outfile->Write();


  	//tree2 -> Draw("AmpYG:AmpXG>>h2");
  	
  	//cout << "Completed writing to file" <<endl;
	cout<<"Fin"<<endl;
	//inputFile->Close();
	//ofile->Close();
	return 0;
	
}
