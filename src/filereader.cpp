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

//generic constructor, reading a file
filereader::filereader(string filename){
	commands = {};
	readFile(filename);
	obsObject = observable();

}

//read a file and save the list of commands in the file
void filereader::readFile(string filename){
	ifstream infile(filename);
    string line;
        
    for(std::string line; getline(infile, line); )
    {
        if(line!="") commands.push_back(line);
    }
    infile.close();

}


//creates all amplitudes specified in the input file, stores in private observable object
void filereader::ConstructBareAmps(){

	vector<channel> chlist = {};
	//makes list of channels
	for(string s : getAddChannelList()){
		chanDat chd = readCh(s);
		channel c = channel(chd.name,chd.ch_masses);
		chlist.push_back(c);
	} 
	
	//get smin and smax
	vector<double> fitreglimits = readFitReg(getFitRegion());
	double smin = pow(fitreglimits[0],2);
	double smax = pow(fitreglimits[1],2);

	vector<ampDat> wvdat_list = {};
	vector<amplitude> amps = {};
	//creates list of amplitude data objects
	for(string s : getAddWaveList()){
		ampDat ampd = readWave(s);
		amplitude amp = amplitude(ampd.wname,ampd.J,ampd.sL,smin,smax,chlist);
		amps.push_back(amp);
	}

	obsObject = observable(amps);
}

void filereader::setChebys(){

	string name1 = "";
	string name2 = "";
	vector<chebyDat> chebydat_list;
	int numAmps = obsObject.getNumAmps();
	//make list of chebyshev data objects
	for(string s : getChebyList()){
		chebyDat chd = readCheby(s);
		for(int i = 0; i < obsObject.amplitudes.size(); i++){
			if(chd.wavename == obsObject.amplitudes[i].getName()){
				int type = 0;
				string my_str = chd.poletype;
				remove(my_str.begin(), my_str.end(), ' ');
				if (my_str == "s") type = 1;
				if (my_str == "p") type = 2;
				if (my_str == "p+s") type = 3;
				obsObject.amplitudes[i].setChebyCoeffs(chd.channame, type, chd.ss0, chd.coeffs);
			}
		}
	}

}

void filereader::setPoles(){

	for(string s : getAddPoleList()){
		poleDat pd = readPole(s);

		for(int i = 0; i < obsObject.amplitudes.size(); i++){
			if(pd.wavename == obsObject.amplitudes[i].getName()){
				obsObject.amplitudes[i].addPole(pd.mass, pd.channames, pd.couplings);
			}
		}
	}

}

void filereader::setKmats(){

	for(string s : getKmatList()){
		kmatDat kmd = readKmat(s);
		
		for(int i = 0; i < obsObject.amplitudes.size(); i++){
			if(kmd.wavename==obsObject.amplitudes[i].getName()){

				vector<vector<double>> kmatelemList = {};
				int counter = 0;
				
				for(int length = obsObject.numChans; length >0; length--){
					vector<double> temp = {};
					for(int k = 0; k<length; k++){
						temp.push_back(kmd.matelems[counter]);
						counter++;
					}
					kmatelemList.push_back(temp);
				}
				obsObject.amplitudes[i].setKParams(kmd.power, kmatelemList);

			}

		}
	}

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

void filereader::SetSeed(){
	regex reg_Seed("SetSeed\\(\\s*([0-9]+)\\s*\\)");
	smatch cmdmatch;
	for(int i = 0; i < commands.size(); i++){
		if(regex_search(commands.at(i), cmdmatch, reg_Seed)){
            Seed = commands[i];
		}
    }		
}

void filereader::SetFitRegion(){
	regex reg_FitRegion("FitRegion\\(\\s*([0-9\\.-]+)\\s*,\\s*([0-9\\.-]+)\\s*\\)"); //FitRegion(0.9975, 2.5)
	smatch cmdmatch;
	for(int i = 0; i < commands.size(); i++){
		if(regex_search(commands.at(i), cmdmatch, reg_FitRegion)){
            FitRegion = commands[i];
		}
    }		
}

void filereader::SetAddChannelList(){
	regex reg_AddChannel("AddChannel\\(\\s*\"(.*?)\"\\s*,\\s*\\{\\s*([0-9\\.-]+)\\s*,\\s*([0-9\\.-]+)\\s*\\}\\s*\\)");
	smatch cmdmatch;
	for(int i = 0; i < commands.size(); i++){
		if(regex_search(commands.at(i), cmdmatch, reg_AddChannel)){
            AddChannel_list.push_back(commands[i]);
		}
    }		
}

void filereader::SetAddWaveList(){
	regex reg_AddWave("AddWave\\(\\s*\"(.*?)\"\\s*,\\s*\"(.*?)\"\\s*,\\s*\"(.*?)\"\\s*,\\s*([0-9\\.-]+)\\s*,\\s*([0-9\\.-]+)\\s*\\)"); //AddWave("S","kmat","nominal", 0, 1.0)
	smatch cmdmatch;
	for(int i = 0; i < commands.size(); i++){
		if(regex_search(commands.at(i), cmdmatch, reg_AddWave)){
            AddWave_list.push_back(commands[i]);
		}
    }		
}

void filereader::SetChebyList(){
	regex reg_Cheby("ChebyCoeffs\\(\\s*\"(.*?)\"\\s*,\\s*\"(.*?)\"\\s*,\\s*\"(.*?)\"\\s*,\\s*\\{((\\s*[0-9\\.-]+\\s*\\\\pm\\s*[0-9\\.]+\\s*,?\\s*)+)\\}\\)");
	smatch cmdmatch;
	for(int i = 0; i < commands.size(); i++){
		if(regex_search(commands.at(i), cmdmatch, reg_Cheby)){
            Cheby_list.push_back(commands[i]);
		}
    }		
}

void filereader::SetAddPoleList(){
	regex reg_AddPole("AddPole\\(\"(.*?)\"\\s*,\\s*([0-9\\.-]+\\s*\\\\pm\\s*[0-9\\.]+)\\s*,\\s*\\{\\s*((\"(.*?)\"\\s*,?\\s*)+)\\}\\s*,\\s*\\{\\s*(([0-9\\.-]+\\s*\\\\pm\\s*[0-9\\.]+\\s*,?\\s*)+)\\}\\s*\\)");
	smatch cmdmatch;
	for(int i = 0; i < commands.size(); i++){
		if(regex_search(commands.at(i), cmdmatch, reg_AddPole)){
            AddPole_list.push_back(commands[i]);
		}
    }		
}

void filereader::SetKmatList(){
	regex reg_Kmat("AddKmatBackground\\(\\s*\"(.*?)\"\\s*,\\s*([0-9\\.-]+)\\s*,\\s*\\{((\\s*\\{((\\s*[0-9\\.-]+\\s*\\\\pm\\s*[0-9\\.]+\\s*,?\\s*)+)\\}\\s*,?\\s*)+)\\}\\s*\\)");
	smatch cmdmatch;
	for(int i = 0; i < commands.size(); i++){
		if(regex_search(commands.at(i), cmdmatch, reg_Kmat)){
            Kmat_list.push_back(commands[i]);
		}
    }		
}

void filereader::SetAllCommandLists(){
	SetSeed();
	SetFitRegion();
	SetAddChannelList();
	SetAddWaveList();
	SetChebyList();
	SetAddPoleList();
	SetKmatList();
}

string filereader::getSeed(){
	return Seed;
}

string filereader::getFitRegion(){
	return FitRegion;
}

vector<string> filereader::getAddChannelList(){
	return AddChannel_list;
}

vector<string> filereader::getAddWaveList(){
	return AddWave_list;
}

vector<string> filereader::getChebyList(){
	return Cheby_list;
}

vector<string> filereader::getAddPoleList(){
	return AddPole_list;
}

vector<string> filereader::getKmatList(){
	return Kmat_list;
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
		remove(chname.begin(), chname.end(), ' ');

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
			remove(wname.begin(), wname.end(), ' ');
			kmatflag = match[2];
			rhoN = match[3];
			J = stod(string(match[4]));
			sL = stod(string(match[5]));
        }
	return ampDat(wname, kmatflag,rhoN,J,sL);
}


chebyDat filereader::readCheby(string cmd){
	regex reg_Cheby("ChebyCoeffs\\(\\s*\"(.*?)\"\\s*,\\s*\"(.*?)\"\\s*,\\s*\"(.*?)\"\\s*,\\s*\\{((\\s*[0-9\\.-]+\\s*\\\\pm\\s*[0-9\\.]+\\s*,?\\s*)+)\\}\\)");
	regex testreg_number("([+-]{0,1}?(?=\\.\\d|\\d)(?:\\d+)?(?:\\.?\\d*))(?:[Ee]([+-]{0,1}?\\d+))?");
	regex testreg_name("[A-Za-z\\+]+");
	string wavename = "";
	string chname = "";
	string ptype = "";
	double ss0 = 0;
	vector<double> coeffs, coeffs_inc={};
	if(regex_search(cmd, match, reg_Cheby)){
		wavename = match[1];
		remove(wavename.begin(), wavename.end(), ' ');
		chname = match[2];
		remove(chname.begin(), chname.end(), ' ');
		
		string use = match[3];
		smatch match2;
        if(regex_search(use, match2, testreg_name)){
            ptype = match2[0];
        }
        if(regex_search(use, match2, testreg_number)){
            ss0 = stod(string(match2[0]));
        }

		string mySuffix, temp;
		mySuffix = match[4];
		int contatore = 0;
        while(regex_search(mySuffix, testmatch, testreg_number))
        {
            temp = testmatch[0];
            mySuffix = testmatch.suffix();
            if (contatore%2 == 0){
                coeffs.push_back(stod(string(temp)));
            }else{
                coeffs_inc.push_back(stod(string(temp)));
            }
            contatore++;
        }

    }

	return chebyDat(wavename, chname, ptype, ss0, coeffs, coeffs_inc);
}


poleDat filereader::readPole(string cmd){
	regex reg_AddPole("AddPole\\(\"(.*?)\"\\s*,\\s*([0-9\\.-]+\\s*\\\\pm\\s*[0-9\\.]+)\\s*,\\s*\\{\\s*((\"(.*?)\"\\s*,?\\s*)+)\\}\\s*,\\s*\\{\\s*(([0-9\\.-]+\\s*\\\\pm\\s*[0-9\\.]+\\s*,?\\s*)+)\\}\\s*\\)"); 
	regex testreg_number("([+-]{0,1}?(?=\\.\\d|\\d)(?:\\d+)?(?:\\.?\\d*))(?:[Ee]([+-]{0,1}?\\d+))?");
	regex testreg_name("[A-Za-z]+");
	
	string wname = "";
	double mass, mass_inc = 0;
	vector<string> chnames = {};
	vector<double> couplings, couplings_inc = {};

	if(regex_search(cmd, match, reg_AddPole)){
		wname = match[1];
		remove(wname.begin(), wname.end(), ' ');
		
		string mySuffix, temp;
		mySuffix = match[2];
		int counter = 0;
        while(regex_search(mySuffix, testmatch, testreg_number))
        {
            temp = testmatch[0];
            mySuffix = testmatch.suffix();
			if (counter%2 == 0){
                mass = stod(string(temp));
            }else{
                mass_inc = stod(string(temp));
            }
			counter++;
        }

		mySuffix = match[3];

        counter = 0;
        while(regex_search(mySuffix, testmatch, testreg_name))
        {
            temp = testmatch[0];
			remove(temp.begin(), temp.end(), ' ');
			chnames.push_back(temp);
            mySuffix = testmatch.suffix();
            counter++;
        }
		
		mySuffix = match[3 + counter];

		counter = 0;
		while(regex_search(mySuffix, testmatch, testreg_number))
		{
			temp = testmatch[0];
			mySuffix = testmatch.suffix();
			if (counter%2 == 0){
				couplings.push_back(stod(string(temp)));
			}else{
				couplings_inc.push_back(stod(string(temp)));
			}
			
			counter++;
		}

	}

	if(chnames.size()!=couplings.size()) return poleDat("",0,0,{},{},{});

	return poleDat(wname, mass, mass_inc, chnames, couplings, couplings_inc);
}



kmatDat filereader::readKmat(string cmd){
	regex reg_Kmat("AddKmatBackground\\(\\s*\"(.*?)\"\\s*,\\s*([0-9\\.-]+)\\s*,\\s*\\{((\\s*\\{((\\s*[0-9\\.-]+\\s*\\\\pm\\s*[0-9\\.]+\\s*,?\\s*)+)\\}\\s*,?\\s*)+)\\}\\s*\\)");
	regex testreg_number("([+-]{0,1}?(?=\\.\\d|\\d)(?:\\d+)?(?:\\.?\\d*))(?:[Ee]([+-]{0,1}?\\d+))?");

	string wname = "";
	int power = 0;
	vector<double> matelems, matelems_inc = {};

	if(regex_search(cmd, match, reg_Kmat)){
		wname = match[1];
		remove(wname.begin(), wname.end(), ' ');
		power = stoi(string(match[2]));

		string mySuffix, temp;
		int cou = 0;
		mySuffix = match[3];

        while(regex_search(mySuffix, testmatch, testreg_number))
        {
            temp = testmatch[0];
            mySuffix = testmatch.suffix();
            if (cou%2 == 0){
                matelems.push_back(stod(string(temp)));
            }else{
                matelems_inc.push_back(stod(string(temp)));
            }
            cou++;
        }

    }

	return kmatDat(wname,power,matelems, matelems_inc);
}


observable filereader::getObs(){

	return obsObject;
}