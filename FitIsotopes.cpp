#ifndef FITISOTOPES_CPP_
#define FITISOTOPES_CPP_

struct isotope{
	isotope():Z(0),A(0),halflife(0),halflifeErr(0),bphalflife(0),bphalflifeErr(0),bp(0),bpErr(0),betaEffRel(1.0),protonEff(1.0),decayType(0){}
	std::string name;
	int Z;
	int A;
	double halflife;
	double halflifeErr;
	double bphalflife;
	double bphalflifeErr;
	double bp;
	double bpErr;
	double betaEffRel;
	double protonEff;
	int decayType;
	bool stable;
};
struct isotopeParameters{
	isotopeParameters():halflife(0),halflifeErr(0),bphalflife(0),bphalflifeErr(0),bp(0),bpErr(0),betaEffRel(0),protonEff(0){}
	std::string name;
	int halflife;
	int halflifeErr;
	int bphalflife;
	int bphalflifeErr;
	int bp;
	int bpErr;
	int betaEffRel;
	int protonEff;
};
int CreateIsotopeList(std::vector< isotope > * isotopeVector, std::string isotopeFile){
//reads in isotope data, constructs isotope object
	std::string line;
	std::ifstream fileIn;
	fileIn.open(isotopeFile.c_str());

	isotope reference;

	if (!fileIn.is_open()){
		std::cerr << " Isotope file not opened exiting program" << std::endl;
		return -1; 
	}
	else {
		std::cout << isotopeFile << " is open. About to begin reading in isotopes"	<<std::endl;
	}

	while (fileIn.good()){
		getline(fileIn,line);
		auto commentLine=line.find("#");
		std::string dummyVar;
		auto newLine=line.substr(0,commentLine);
		if(newLine.size()>0){
			std::istringstream iss(line,std::istringstream::in);

			iss >> reference.name;
			iss >> reference.Z;
			iss >> reference.A;
			iss >> reference.halflife;
			iss >> reference.halflifeErr;
			iss >> reference.bphalflife;
			iss >> reference.bphalflifeErr;
			iss >> reference.bp;
			iss >> reference.bpErr;
			//iss >> reference.betaEffRel;
			//iss >> reference.protonEff;

			isotopeVector->push_back(reference);

		}
	}
	return 0;

}
int FindIsotopeLocation(std::vector<isotope> * isotopeList,int A, int Z){

	for(int i = 0; i < isotopeList->size();i++ ){
		if(isotopeList->at(i).A == A && isotopeList->at(i).Z == Z){
			return i;
		}
	}
	return -1;
}
std::vector<isotope>::iterator FindIsotopeLocation(std::vector<isotope> * isotopeList,std::string name){
	std::vector<isotope>::iterator isotopeIt;

	for(isotopeIt = isotopeList->begin(); isotopeIt != isotopeList->end();isotopeIt++ ){
		if(isotopeIt->name == name){
			return isotopeIt;
		}
	}
	return isotopeList->end();
}
int DecayChain(std::vector<std::list<isotope>> * decayChain, std::vector<isotope> * isotopeList, int Z, int A, int callingBranch){

	std::list<isotope> currentChain;
	int isotopeLocation;
	int currentZ = Z;
	int currentA = A;
	int chainNum = decayChain->size();
	bool stability = false;

	std::cout << "chain Num" << chainNum << " Size" << decayChain->size() << std::endl;

	if(callingBranch < 0){
	decayChain->push_back(currentChain);
	chainNum = 0;
	}
	if(callingBranch >= 0){
		//Copy the current part of the chain over to the new one as well
		decayChain->push_back(decayChain->at(callingBranch)); 
	}
	std::cout << "chain Num" << chainNum << " Size" << decayChain->size() << std::endl;


	while(!stability){
		isotopeLocation = FindIsotopeLocation(isotopeList, currentA, currentZ);
		if(isotopeLocation <0 || isotopeLocation >=isotopeList->size()){
			std::cout << "Isotope not located in vector" <<std::endl;
			std::cout << "z:" << currentZ << " a:" << currentA << std::endl;
			return -1;
		}
		decayChain->at(chainNum).push_back(isotopeList->at(isotopeLocation));
		stability = isotopeList->at(isotopeLocation).stable;
	
		if(isotopeList->at(isotopeLocation).bp>0.0){
			decayChain->at(chainNum).back().decayType=1;
			//Start decay chain along new branch
			DecayChain(decayChain, isotopeList, currentZ-1,currentA-1,chainNum);
		}
		
		currentZ = currentZ-1;
		decayChain->at(chainNum).back().decayType=0;
	}
	return 0;

}
void PrintDecayChain(std::vector<std::list<isotope>> * decayChain){
	//Print the isotopes in the decay chain
	std::list<isotope>::iterator isotopeIt;
	for(int i = 0; i < decayChain->size();i++){
		std::cout << "\nPrinting decay path: " << i << std::endl;
		for(isotopeIt = decayChain->at(i).begin(); isotopeIt != decayChain->at(i).end(); isotopeIt++){
			std::cout << isotopeIt->name << " " << isotopeIt->decayType << " -> ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	return;
}
//check if isotope parameters have already been declared or not
bool AlreadyDeclared(std::vector<isotopeParameters> * parametersVector,std::string isotopeName){
	if(parametersVector->size() == 0){
		return false;
	}
	for(int isotope = 0; isotope < parametersVector->size();isotope++){
		if(isotopeName == parametersVector->at(isotope).name){
			return true;
		}
	}
	return false;
}
void DeclareParameters(std::vector<std::list<isotope>> * decayChain, std::vector<isotopeParameters> * parametersVector){
	//Function to set parameters for the isotopes in the decay chain. As isotopes may appear multiple times
	//sets the parameters ahead of time so the parameters are shared
	std::list<isotope>::iterator isotopeIt;
	isotopeParameters isoPar;
	int currentPar = 4;

	for(int path = 0; path < decayChain->size();path++){
		for(isotopeIt = decayChain->at(path).begin();isotopeIt != decayChain->at(path).end(); isotopeIt++){

			if(!AlreadyDeclared(parametersVector,isotopeIt->name)){

				isoPar.name = isotopeIt->name;
				isoPar.halflife = currentPar + 0;
				isoPar.halflifeErr = currentPar + 1;
				isoPar.bphalflife = currentPar + 2;
				isoPar.bphalflifeErr = currentPar + 3;
				isoPar.bp = currentPar + 4;
				isoPar.bpErr = currentPar + 5;
				isoPar.betaEffRel = currentPar + 6;
				isoPar.protonEff = currentPar + 7;

				currentPar = currentPar + 8;
				parametersVector->push_back(isoPar);
			}

		}
	}
	return;
}
void ReadParameters(std::vector<isotopeParameters> * parametersVector){
	//Reads the range of parameters covered by each isotope
	for(int i = 0; i < parametersVector->size(); i++){
		std::cout << parametersVector->at(i).name << " Par " << parametersVector->at(i).halflife << " - " << parametersVector->at(i).protonEff << std::endl;
	}
}
std::vector<isotopeParameters>::iterator GetIsotopeParamPos(std::vector<isotopeParameters> * parametersVector, std::string name){
	//Returns the position of an isotope in the parameters vector
	std::vector<isotopeParameters>::iterator parPos;
	for(parPos = parametersVector->begin(); parPos != parametersVector->end(); parPos++ ){
		if(name == parPos->name){
			return parPos;
		}
	}
	return (parametersVector->end()++);
}
void SetFitParameters(TF1 * function, std::vector<isotope> * isotopeList,std::vector<isotopeParameters> * parametersVector){
	//Set the parameters for the fit function
	std::vector<isotopeParameters>::iterator isotopePar;
	std::vector<isotope>::iterator isotopeValues;
	function->SetParName(0,"background slope");
	function->SetParName(1,"background intercept");
	function->SetParName(2,"N0");
	function->SetParName(3,"beta efficiency");
	function->SetParName(4,"Halflife");

	for(isotopePar = parametersVector->begin(); isotopePar != parametersVector->end();isotopePar++){
		isotopeValues = FindIsotopeLocation(isotopeList, isotopePar->name);

		if(isotopePar == parametersVector->begin()){
			function->SetParameter(isotopePar->halflife,isotopeValues->halflife);
			function->FixParameter(isotopePar->halflifeErr,isotopeValues->halflifeErr);
			function->SetParameter(isotopePar->bphalflife,isotopeValues->bphalflife);
			function->FixParameter(isotopePar->bphalflifeErr,isotopeValues->bphalflifeErr);
			function->FixParameter(isotopePar->bp,isotopeValues->bp);
			function->FixParameter(isotopePar->bpErr,isotopeValues->bpErr);
			function->FixParameter(isotopePar->betaEffRel,isotopeValues->betaEffRel);
			function->FixParameter(isotopePar->protonEff,isotopeValues->protonEff);
		}
		else{
			function->SetParameter(isotopePar->halflife,isotopeValues->halflife);
			function->FixParameter(isotopePar->halflifeErr,isotopeValues->halflifeErr);
			function->SetParameter(isotopePar->bphalflife,isotopeValues->bphalflife);
			function->FixParameter(isotopePar->bphalflifeErr,isotopeValues->bphalflifeErr);
			function->FixParameter(isotopePar->bp,isotopeValues->bp);
			function->FixParameter(isotopePar->bpErr,isotopeValues->bpErr);
			function->FixParameter(isotopePar->betaEffRel,isotopeValues->betaEffRel);
			function->FixParameter(isotopePar->protonEff,isotopeValues->protonEff);
		}

	}
	return;
}


#endif