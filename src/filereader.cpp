#include "filereader.h"
#include<complex>
#include<algorithm>
#include<math.h>
#include<chrono>
#include<random>
#include<fstream>
#include<vector>
#include<string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>

using namespace std;

filereader::filereader(string filename){
	commands = {};
	readFile(filename);

}


void filereader::readFile(string filename){
	ifstream infile(filename);
    string line;
        
    for(std::string line; getline(infile, line); )
    {
        if(line!="") commands.push_back(line);
    }
    infile.close();

}


vector<amplitude> filereader::constructAmps(){
	regex reg_Seed("SetSeed\\(\\s*([0-9]+)\\s*\\)"); //SetSeed(3)
    regex reg_FitRegion("FitRegion\\(\\s*([0-9\\.-]+)\\s*,\\s*([0-9\\.-]+)\\s*\\)"); //FitRegion(0.9975, 2.5)
    regex reg_AddChannel("AddChannel\\(\\s*\"(.*?)\"\\s*,\\s*\\{\\s*([0-9\\.-]+)\\s*,\\s*([0-9\\.-]+)\\s*\\}\\s*\\)");
    regex reg_AddWave("AddWave\\(\\s*\"(.*?)\"\\s*,\\s*\"(.*?)\"\\s*,\\s*\"(.*?)\"\\s*,\\s*([0-9\\.-]+)\\s*,\\s*([0-9\\.-]+)\\s*\\)"); //AddWave("S","kmat","nominal", 0, 1.0)
    regex reg_Cheby("ChebyCoeffs\\(\\s*\"(.*?)\"\\s*,\\s*\"(.*?)\"\\s*,\\s*\"(.*?)\"\\s*,\\s*\\{((\\s*[0-9\\.-]+\\s*\\\\pm\\s*[0-9\\.]+\\s*,?\\s*)+)\\}\\)");
    regex reg_AddPole("AddPole\\(\"(.*?)\"\\s*,\\s*([0-9\\.-]+\\s*\\\\pm\\s*[0-9\\.]+)\\s*,\\s*\\{\\s*(((.*?)\\s*,?\\s*)+)\\}\\s*,\\s*\\{\\s*(([0-9\\.-]+\\s*\\\\pm\\s*[0-9\\.]+\\s*,?\\s*)+)\\}\\s*\\)"); 
    regex reg_Kmat("AddKmatBackground\\(\\s*\"(.*?)\"\\s*,\\s*[0-9\\.-]+\\s*,\\s*\\{((\\s*\\{((\\s*[0-9\\.-]+\\s*\\\\pm\\s*[0-9\\.]+\\s*,?\\s*)+)\\}\\s*,?\\s*)+)\\}\\s*\\)");
	
	return {};
}

void filereader::printCommands(){
	for(string s : commands){

		cout<<s<<endl;
	}

}

string filereader::getCommand(int i){
	if(i>=0 && i < commands.size()) return commands[i];

	return "";
}


int filereader::readSeed(string cmd){
	regex reg_Seed("SetSeed\\(\\s*([0-9]+)\\s*\\)");

	if(regex_search(cmd, match, reg_Seed)){
        return stoi(string(match[1]));
    }

	return 0;
}



vector<double> filereader::readFitReg(string cmd){
	regex reg_FitRegion("FitRegion\\(\\s*([0-9\\.-]+)\\s*,\\s*([0-9\\.-]+)\\s*\\)"); //FitRegion(0.9975, 2.5)
	double smin = 0;
	double smax = 0;
	if(regex_search(cmd, match, reg_FitRegion)){
        smin = stod(string(match[1]));
        smax = stod(string(match[2]));
    }

	return {smin,smax};
}


//only works for channels with two masses (2 particle states)
//Need to edit regex
chanDat filereader::readCh(string cmd){
	regex reg_AddChannel("AddChannel\\(\\s*\"(.*?)\"\\s*,\\s*\\{\\s*([0-9\\.-]+)\\s*,\\s*([0-9\\.-]+)\\s*\\}\\s*\\)");
	string chname ="";
	vector<double> masses={};
	
	if(regex_search(cmd, match, reg_AddChannel)){
        chname = match[1];

		for(int i = 2; i < match.size(); i++){
			masses.push_back(stod(string(match[i])));

		}
    }
	return chanDat(chname,masses);
}



ampDat filereader::readWave(string cmd){
	regex reg_AddWave("AddWave\\(\\s*\"(.*?)\"\\s*,\\s*\"(.*?)\"\\s*,\\s*\"(.*?)\"\\s*,\\s*([0-9\\.-]+)\\s*,\\s*([0-9\\.-]+)\\s*\\)"); //AddWave("S","kmat","nominal", 0, 1.0)
	string wname = "void";
	string kmatflag = "";
	string rhoN = "";
	int J = 0;
	double sL = 0;


	if(regex_search(cmd, match, reg_AddWave)){
			wname = match[1];
			kmatflag = match[2];
			rhoN = match[3];
			J = stod(string(match[4]));
			sL = stod(string(match[5]));
        }
	return ampDat(wname, kmatflag,rhoN,J,sL);
}


chebyDat filereader::readCheby(string cmd){
	regex reg_Cheby("ChebyCoeffs\\(\\s*\"(.*?)\"\\s*,\\s*\"(.*?)\"\\s*,\\s*\"(.*?)\"\\s*,\\s*\\{((\\s*[0-9\\.-]+\\s*\\\\pm\\s*[0-9\\.]+\\s*,?\\s*)+)\\}\\)");
	string wavename = "";
	string chname = "";
	string ptype = "";
	vector<double> coeffs={};
	if(regex_search(cmd, match, reg_Cheby)){
		wavename = match[1];
		chname = match[2];
		ptype = match[3];
		for(int i = 4; i < match.size(); i++){
				coeffs.push_back(stod(string(match[i])));
		}
    }


	return chebyDat(wavename, chname, ptype, coeffs);
}


poleDat filereader::readPole(string cmd){
	regex reg_AddPole("AddPole\\(\"(.*?)\"\\s*,\\s*([0-9\\.-]+\\s*\\\\pm\\s*[0-9\\.]+)\\s*,\\s*\\{\\s*(((.*?)\\s*,?\\s*)+)\\}\\s*,\\s*\\{\\s*(([0-9\\.-]+\\s*\\\\pm\\s*[0-9\\.]+\\s*,?\\s*)+)\\}\\s*\\)"); 

	string wname = "";
	double mass = 0;
	vector<string> chnames = {};
	vector<double> couplings = {};

	if(regex_search(cmd, match, reg_AddPole)){
		wname = match[1];
		mass = stod(string(match[2]));
		int numCoupls = (match.size()-3)/2;

		for(int i = 3; i < numCoupls+3; i++){
				chnames.push_back(string(match[i]));
		}

		for(int i = numCoupls +3; i<match.size(); i++){

			couplings.push_back(stod(string(match[i])));
		}

    }

	if(chnames.size()!=couplings.size()) return poleDat("",0,{},{});

	return poleDat(wname, mass, chnames, couplings);
}



kmatDat filereader::readKmat(string cmd){
	regex reg_Kmat("AddKmatBackground\\(\\s*\"(.*?)\"\\s*,\\s*[0-9\\.-]+\\s*,\\s*\\{((\\s*\\{((\\s*[0-9\\.-]+\\s*\\\\pm\\s*[0-9\\.]+\\s*,?\\s*)+)\\}\\s*,?\\s*)+)\\}\\s*\\)");

	string wname = "";
	int power = 0;
	vector<double> matelems = {};

	if(regex_search(cmd, match, reg_Kmat)){
		wname = match[1];
		power = stoi(string(match[2]));

		for(int i = 3; i < match.size(); i++){
				matelems.push_back(stod(string(match[i])));

		}

    }

	return kmatDat(wname,power,matelems);
}