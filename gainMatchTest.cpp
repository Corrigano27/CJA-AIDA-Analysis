#include "TFile.h"
#include "TTree.h"
#include "TList.h"
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
#include <limits>

using namespace std;

int gainMatchTest(std::string fileName){

    TFile * inputFile = TFile::Open(fileName.c_str(),"read"); 
	TTree * inputTTree = (TTree*)inputFile->Get("AIDA_hits");

	std::string oName;

	oName+="_GainMatchTest.root";
  	//Open the tree and create the branch to write to
  	TFile * ofile = TFile::Open( oName.c_str(), "recreate");

	std::vector<std::pair<double, double>> amplitudes[6*128*128];
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
  inputTTree->SetBranchAddress("aida_hit", & (inputEntry.T) );
  Long64_t nEntries = inputTTree->GetEntries();
  std::cout << "Nentries" << nEntries << std::endl;
  double dx; 
  double dy;

  std::pair<double, double> ampPair;

  std::string histogramName;

  	TH2D * exEy[6];
	//TH2D * exEyPre = new TH2D("ExEyPre","",500,0,20e3,500,0,20e3);
	TH2D * diffStripX[6];
	TH2D * diffStripY[6];
	TH1D * exEyDiff[6];

	for(int j = 0; j < 6 ; j++){	

		char hname[50];
		sprintf(hname,"ExEyDet%i",j);
		histogramName = hname;
		exEy[j] = new TH2D(histogramName.c_str(),"",500,0,4e3,500,0,4e3);
		exEy[j]->SetDirectory(ofile);	

		hname[0] = '\0';
		sprintf(hname,"DiffStripXDet%i",j);
		histogramName = hname;
		diffStripX[j] = new TH2D(histogramName.c_str(),"",128,0,128,1200,-600,600);
		diffStripX[j]->SetDirectory(ofile);	

		hname[0] = '\0';
		sprintf(hname,"DiffStripYDet%i",j);
		histogramName = hname;
		diffStripY[j] = new TH2D(histogramName.c_str(),"",128,0,128,1200,-600,600);
		diffStripY[j]->SetDirectory(ofile);

		hname[0] = '\0';
		sprintf(hname,"exEyDiffDet%i",j);
		histogramName = hname;
		exEyDiff[j] = new TH1D(histogramName.c_str(),"",1200,-600,600);
		exEyDiff[j]->SetDirectory(ofile);

	}


  	for( Long64_t iEntry = 0; iEntry < nEntries; ++iEntry ){
    inputTTree->GetEntry(iEntry);
    //std::cout << inputEntry.x << std::endl;
		if(inputEntry.ID==5 && inputEntry.E>200){
			if(inputEntry.z==0){
				diffStripX[0]->Fill(inputEntry.x,inputEntry.EX-inputEntry.EY);
				diffStripY[0]->Fill(inputEntry.y,inputEntry.EX-inputEntry.EY);
				exEy[0]->Fill(inputEntry.EX,inputEntry.EY);
				//exEyPre->Fill(Ax2,Ay2);
				exEyDiff[0]->Fill(inputEntry.EX-inputEntry.EY);
			}
			else if(inputEntry.z==1){
				diffStripX[1]->Fill(inputEntry.x,inputEntry.EX-inputEntry.EY);
				diffStripY[1]->Fill(inputEntry.y,inputEntry.EX-inputEntry.EY);
				exEy[1]->Fill(inputEntry.EX,inputEntry.EY);
				//exEyPre->Fill(Ax2,Ay2);
				exEyDiff[1]->Fill(inputEntry.EX-inputEntry.EY);
			}
			else if(inputEntry.z==2){
				diffStripX[2]->Fill(inputEntry.x,inputEntry.EX-inputEntry.EY);
				diffStripY[2]->Fill(inputEntry.y,inputEntry.EX-inputEntry.EY);
				exEy[2]->Fill(inputEntry.EX,inputEntry.EY);
				//exEyPre->Fill(Ax2,Ay2);
				exEyDiff[2]->Fill(inputEntry.EX-inputEntry.EY);
			}
			else if(inputEntry.z==3){
				diffStripX[3]->Fill(inputEntry.x,inputEntry.EX-inputEntry.EY);
				diffStripY[3]->Fill(inputEntry.y,inputEntry.EX-inputEntry.EY);
				exEy[3]->Fill(inputEntry.EX,inputEntry.EY);
				//exEyPre->Fill(Ax2,Ay2);
				exEyDiff[3]->Fill(inputEntry.EX-inputEntry.EY);
			}
			else if(inputEntry.z==4){
				diffStripX[4]->Fill(inputEntry.x,inputEntry.EX-inputEntry.EY);
				diffStripY[4]->Fill(inputEntry.y,inputEntry.EX-inputEntry.EY);
				exEy[4]->Fill(inputEntry.EX,inputEntry.EY);
				//exEyPre->Fill(Ax2,Ay2);
				exEyDiff[4]->Fill(inputEntry.EX-inputEntry.EY);
			}
			else if(inputEntry.z==5){
				diffStripX[5]->Fill(inputEntry.x,inputEntry.EX-inputEntry.EY);
				diffStripY[5]->Fill(inputEntry.y,inputEntry.EX-inputEntry.EY);
				exEy[5]->Fill(inputEntry.EX,inputEntry.EY);
				//exEyPre->Fill(Ax2,Ay2);
				exEyDiff[5]->Fill(inputEntry.EX-inputEntry.EY);
			}

    	//ampPair.first = inputEntry.EX;
    	//ampPair.second = inputEntry.EY;
    	//amplitudes[(int)((inputEntry.z*128*128)+((inputEntry.y)*128)+inputEntry.x)].push_back(ampPair);
		}

		
	}

	for (int dssd = 0; dssd<6; dssd++){
			diffStripX[dssd]->Write();		
			diffStripY[dssd]->Write();	
			exEy[dssd]->Write();
			exEyDiff[dssd]->Write();
	}

  



  


	cout<<"Fin"<<endl;
	inputFile->Close();
	ofile->Close();
	//outf.Close();
	return 0;

}