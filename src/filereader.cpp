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
#include <TRandom3.h>

using namespace std;

//generic constructor, reading a file
filereader::filereader(string filename){
	commands = {};
	readFile(filename);
	obsObject = observable();
	NameOfFile = filename;
	seed = 0;
}

//read a file and save the list of commands in the file
void filereader::readFile(string filename){
	ifstream infile(filename);
	if (!infile) { cout << "Error: file " << filename << " not found!!!" << endl; exit(0); }

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
		amplitude amp = amplitude(ampd.wname,ampd.J,ampd.sL,smin,smax,chlist,ampd.kmatflag,ampd.rhoN);
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
	int counter = 0;
	for(string s : getChebyList()){
		chebyDat chd = readCheby(s);
		for(int i = 0; i < obsObject.amplitudes.size(); i++){
			if(chd.wavename == obsObject.amplitudes[i].getName()){
				int type = 0;
				string my_str = chd.poletype;
				remove(my_str.begin(), my_str.end(), ' ');
				if (my_str == "p") type = 1;
				if (my_str == "s") type = 2;
				if (my_str == "p+s") type = 3;
				obsObject.amplitudes[i].setChebyCoeffs(chd.channame, type, chd.ss0, chd.coeffs,chd.coeffs_inc);
			}
		}
		
	}

}

void filereader::setPoles(){

	for(string s : getAddPoleList()){
		poleDat pd = readPole(s);
		for(int i = 0; i < obsObject.amplitudes.size(); i++){
			if(pd.wavename == obsObject.amplitudes[i].getName()){
				obsObject.amplitudes[i].addPole(pd.mass,pd.mass_inc, pd.channames, pd.couplings,pd.couplings_inc);
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
				obsObject.amplitudes[i].setKParams(kmd.power, kmatelemList,kmd.matelems_inc);

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

void filereader::SetChi2CutOff(){
	regex reg_Chi2CutOff("ReducedChi2CutOff\\(\\s*([0-9\\.-]+)\\s*\\)");
	smatch cmdmatch;
	for(int i = 0; i < commands.size(); i++){
		if(regex_search(commands.at(i), cmdmatch, reg_Chi2CutOff)){
            Chi2CutOffCmd = commands[i];
		}
    }
}

void filereader::SetInclChi2Weight(){
	regex reg_InclChi2Weight("InclChi2Weight\\(\\s*([0-9\\.-]+)\\s*\\)");
	smatch cmdmatch;
	for(int i = 0; i < commands.size(); i++){
		if(regex_search(commands.at(i), cmdmatch, reg_InclChi2Weight)){
            InclChi2WeightCmd = commands[i];
		}
    }
}

void filereader::SetExclChi2Weight(){
	regex reg_ExclChi2Weight("ExclChi2Weight\\(\\s*([0-9\\.-]+)\\s*\\)");
	smatch cmdmatch;
	for(int i = 0; i < commands.size(); i++){
		if(regex_search(commands.at(i), cmdmatch, reg_ExclChi2Weight)){
            ExclChi2WeightCmd = commands[i];
		}
    }
}

void filereader::SetActionCmd(){
	regex reg_action("ChooseAnAction\\(\\s*\"(.*?)\"\\s*\\)");
	smatch cmdmatch;
	for(int i = 0; i < commands.size(); i++){
		if(regex_search(commands.at(i), cmdmatch, reg_action)){
            ActionFlagCmd = commands[i];
		}
    }
}

void filereader::SetRandomizeFlag(){
	regex reg_Randomize("DoRandomize\\(\\s*(.*?)\\s*\\)");
	smatch cmdmatch;
	for(int i = 0; i < commands.size(); i++){
		if(regex_search(commands.at(i), cmdmatch, reg_Randomize)){
            RandomizeFlagCmd = commands[i];
		}
    }
}

void filereader::SetInclCrossSecFlag(){
	regex reg_InclCrossSecFlag("IncludeAlsoInclusiveCrossSection\\(\\s*(.*?)\\s*\\)");
	smatch cmdmatch;
	for(int i = 0; i < commands.size(); i++){
		if(regex_search(commands.at(i), cmdmatch, reg_InclCrossSecFlag)){
            InclCrossSecFlagCmd = commands[i];
		}
    }
}

void filereader::SetSeedCmd(){
	regex reg_Seed("SetSeed\\(\\s*([0-9]+)\\s*\\)");
	smatch cmdmatch;
	for(int i = 0; i < commands.size(); i++){
		if(regex_search(commands.at(i), cmdmatch, reg_Seed)){
            SeedCmd = commands[i];
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

void filereader::SetGrid(){
	regex reg_grid("PolesearchGrid\\(\\s*([0-9\\.-]+)\\s*,\\s*([0-9\\.-]+)\\s*,\\s*([0-9\\.-]+)\\s*,\\s*([0-9\\.-]+)\\s*,\\s*([0-9\\.-]+)\\s*,\\s*([0-9\\.-]+)\\s*\\)"); //PolesearchGrid(113., 121., 31, -1.5, 1.5, 31)
	smatch cmdmatch;
	for(int i = 0; i < commands.size(); i++){
		if(regex_search(commands.at(i), cmdmatch, reg_grid)){
            GridCmd = commands[i]; 
		}
    }
}

void filereader::SetZeroCmd(){
	regex reg_polezero("\\s*PolesearchZero\\(\\s*([0-9\\.-]+)\\s*\\)\\s*"); //PolesearchZero(-7)
	smatch cmdmatch;
	for(int i = 0; i < commands.size(); i++){
		if(regex_search(commands.at(i), cmdmatch, reg_polezero)){
            ZeroCmd = commands[i]; 
		}
    }
}

void filereader::SetFitParamsCmd(){
	regex reg_fitparams("\\s*FittingParameters\\(((\\s*\\{((\\s*[0-9\\.-]+\\s*,?\\s*)+)\\}\\s*,?\\s*)+)\\)\\s*"); //FittingParameters({100000, 10000, 0.001, 1}, {100000, 10000, 0.001, 1})
	smatch cmdmatch;//((\\s*\\{(.*?)\\}\\s*,?\\s*)+)
	for(int i = 0; i < commands.size(); i++){
		if(regex_search(commands.at(i), cmdmatch, reg_fitparams)){
            fitparamsCmd = commands[i]; 
		}
    }
}

void filereader::SetFitSequence(){
	regex reg_fitseq("FittingSequence\\(\\s*\\{\\s*((\"(.*?)\"\\s*,?\\s*)+)\\}\\s*\\)");
	smatch cmdmatch;
	for(int i = 0; i < commands.size(); i++){
		if(regex_search(commands.at(i), cmdmatch, reg_fitseq)){
            FitSequenceCmd = commands[i]; 
		}
    }
}

void filereader::setExpInclCrossSec(){
	regex reg_inclcrosssec("LoadExpInclusiveCrossSection\\(\\s*\"(.*?)\"\\s*\\)");
	smatch cmdmatch;
	for(int i = 0; i < commands.size(); i++){
		if(regex_search(commands.at(i), cmdmatch, reg_inclcrosssec)){
            ExpInclCrossSecCmd = commands[i];
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
			//cout << commands[i] << endl;
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
	SetSeedCmd();
	SetFitRegion();
	SetGrid();
	SetZeroCmd();
	SetFitParamsCmd();
	SetFitSequence();
	SetChi2CutOff();
	SetInclChi2Weight();
	SetExclChi2Weight();
	SetActionCmd();
	setExpInclCrossSec();
	SetInclCrossSecFlag();
	SetRandomizeFlag();
	SetAddChannelList();
	SetAddWaveList();
	SetChebyList();
	SetAddPoleList();
	SetKmatList();
	SetExpDataList();
}

double filereader::getChi2CutOff(){
	return readChi2CutOffCmd(Chi2CutOffCmd);
}

double filereader::GetInclChi2Weight(){
	return ReadInclChi2WeightCmd(InclChi2WeightCmd);
}

double filereader::GetExclChi2Weight(){
	return ReadExclChi2WeightCmd(ExclChi2WeightCmd);
}

bool filereader::getFitFlag(){
	readActionCmd(ActionFlagCmd);
	return FitFlag;
}

bool filereader::getPlotFlag(){
	readActionCmd(ActionFlagCmd);
	return PlotFlag;
}

bool filereader::getPolesearchFlag(){
	readActionCmd(ActionFlagCmd);
	return PolesearchFlag;
}

bool filereader::getInclCrossSecFlag(){
	readInclCrossSecFlag(InclCrossSecFlagCmd);
	return InclCrossSecFlag;
}

bool filereader::getRandomizeFlag(){
	readRandomizeFlag(RandomizeFlagCmd);
	return RandomizeFlag;
}

int filereader::getSeed(){
	return seed;
}

//void filereader::setSeed(int newseed){
//	seed=newseed;
//}

string filereader::getFitRegion(){
	return FitRegion;
}

string filereader::getFitSequence(){
	return FitSequenceCmd;
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

vector<string> filereader::getExpDataList(){
	return ExpData_list;
}

int filereader::readSeed(){
	regex reg_Seed("SetSeed\\(\\s*([0-9]+)\\s*\\)");

	if(regex_search(SeedCmd, match, reg_Seed)){
        seed = stoi(string(match[1]));
    }

	return seed;
}

void filereader::readActionCmd(string cmd){
	regex reg_action("ChooseAnAction\\(\\s*\"(.*?)\"\\s*\\)");

	if(regex_search(cmd, match, reg_action)){
		string my_str = match[1];
		remove(my_str.begin(), my_str.end(), ' '); 
		if (my_str == "Fit"){
			FitFlag = true;
		}
		if (my_str == "Polesearch"){
			PolesearchFlag = true;
		}
		if (my_str == "Plot"){
			PlotFlag = true;
		}
	}

}

void filereader::readRandomizeFlag(string cmd){
	regex reg_Randomize("DoRandomize\\(\\s*(.*?)\\s*\\)");

	if(regex_search(cmd, match, reg_Randomize)){
		string my_str = match[1];
		remove(my_str.begin(), my_str.end(), ' ');
		if (my_str == "No"){
			RandomizeFlag = false;
		}
	}
}

void filereader::readInclCrossSecFlag(string cmd){
	regex reg_InclCrossSecFlag("IncludeAlsoInclusiveCrossSection\\(\\s*(.*?)\\s*\\)");

	if(regex_search(cmd, match, reg_InclCrossSecFlag)){
		string my_str = match[1];
		remove(my_str.begin(), my_str.end(), ' ');
		if (my_str == "No"){
			InclCrossSecFlag = false;
		}
	}
}

double filereader::readChi2CutOffCmd(string cmd){
	regex reg_Chi2CutOff("ReducedChi2CutOff\\(\\s*([0-9\\.-]+)\\s*\\)");

	if(regex_search(cmd, match, reg_Chi2CutOff)){
		return stod(string(match[1]));
	}

	return 0;
}

double filereader::ReadInclChi2WeightCmd(string cmd){
	regex reg_InclChi2Weight("InclChi2Weight\\(\\s*([0-9\\.-]+)\\s*\\)");

	if(regex_search(cmd, match, reg_InclChi2Weight)){
		return stod(string(match[1]));
	}

	return 0;
}

double filereader::ReadExclChi2WeightCmd(string cmd){
	regex reg_ExclChi2Weight("ExclChi2Weight\\(\\s*([0-9\\.-]+)\\s*\\)");

	if(regex_search(cmd, match, reg_ExclChi2Weight)){
		return stod(string(match[1]));
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

gridDat filereader::getGrid(){
	regex reg_grid("PolesearchGrid\\(\\s*([0-9\\.-]+)\\s*,\\s*([0-9\\.-]+)\\s*,\\s*([0-9\\.-]+)\\s*,\\s*([0-9\\.-]+)\\s*,\\s*([0-9\\.-]+)\\s*,\\s*([0-9\\.-]+)\\s*\\)"); //PolesearchGrid(113., 121., 31, -1.5, 1.5, 31)
	double grid_Re_sx = 0;
	double grid_Re_dx = 0;
	int grid_Re_numpts = 0;
	double grid_Im_sx = 0;
	double grid_Im_dx = 0;
	int grid_Im_numpts = 0;
	if(regex_search(GridCmd, match, reg_grid)){
		grid_Re_sx = stod(string(match[1]));
		grid_Re_dx = stod(string(match[2]));
		grid_Re_numpts = stod(string(match[3]));
		grid_Im_sx = stod(string(match[4]));
		grid_Im_dx = stod(string(match[5]));
		grid_Im_numpts = stod(string(match[6]));
    }

	return gridDat(grid_Re_sx, grid_Re_dx, grid_Re_numpts, grid_Im_sx, grid_Im_dx, grid_Im_numpts);
}

double filereader::getZero(){
	regex reg_polezero("\\s*PolesearchZero\\(\\s*([0-9\\.-]+)\\s*\\)\\s*"); //PolesearchZero(-7)
	double zero = 0;
	
	if(regex_search(ZeroCmd, match, reg_polezero)){
		zero = stod(string(match[1]));
    }

	return zero;
}

vector<fitparamsDat> filereader::getFitParams(){

	regex reg_fitparams("\\s*FittingParameters\\(((\\s*\\{((\\s*[0-9\\.-]+\\s*,?\\s*)+)\\}\\s*,?\\s*)+)\\)\\s*"); //FittingParameters({100000, 10000, 0.001, 1}, {100000, 10000, 0.001, 1})
	regex testreg_number("([+-]{0,1}?(?=\\.\\d|\\d)(?:\\d+)?(?:\\.?\\d*))(?:[Ee]([+-]{0,1}?\\d+))?");

	vector<string> test = {};

	vector<fitparamsDat> vec = {};

	int cou = 0;
	
	if(regex_search(fitparamsCmd, match, reg_fitparams)){

		string mySuffix, temp;
		mySuffix = match[1];

		//cout << match[1] << endl;
		//cout << match[2] << endl;
		//cout << match[3] << endl;
		//cout << match[4] << endl;

		fitparamsDat aux(0,0,0.0,0);

        while(regex_search(mySuffix, testmatch, testreg_number))
        {
            temp = testmatch[0];
            mySuffix = testmatch.suffix();
            //if (cou%4 == 0){
            //    aux.maxfuncalls = stoi(string(temp));cout << aux.maxfuncalls << endl;
            //}
			//if (cou%4 == 1){
            //    aux.maxiter = stoi(string(temp));cout << aux.maxiter << endl;
            //}
			//if (cou%4 == 2){
            //    aux.tol = stod(string(temp));cout << aux.tol << endl;
            //}
			//if (cou%4 == 3){
            //    aux.verbose = stoi(string(temp));cout << aux.verbose << endl;
            //}
			//vec.push_back(aux);
            //cou++;
			switch (cou % 4) {
                case 0:
                    aux.maxfuncalls = stoi(temp);
                    //std::cout << "maxfuncalls: " << aux.maxfuncalls << std::endl;
                    break;
                case 1:
                    aux.maxiter = stoi(temp);
                    //std::cout << "maxiter: " << aux.maxiter << std::endl;
                    break;
                case 2:
                    aux.tol = stod(temp);  // usa stod per gestire il tipo double
                    //std::cout << "tol: " << aux.tol << std::endl;
                    break;
                case 3:
                    aux.verbose = stoi(temp);
                    //std::cout << "verbose: " << aux.verbose << std::endl;
                    vec.push_back(aux);  // Aggiungi solo quando tutti i campi sono popolati
                    aux = fitparamsDat(0, 0, 0.0, 0);  // Reinizializza
                    break;
            }
            cou++;
        }

    }

	return vec;
}

vector<string> filereader::readFitSequence(string cmd){
	regex reg_fitseq("FittingSequence\\(\\s*\\{\\s*((\"(.*?)\"\\s*,?\\s*)+)\\}\\s*\\)");
	regex testreg_name("[A-Za-z\\_]+");

	string mySuffix, temp;
	vector<string> vecstr = {};

	if(regex_search(cmd, match, reg_fitseq)){
	
		mySuffix = match[1];

		while(regex_search(mySuffix, testmatch, testreg_name))
		{
			temp = testmatch[0];
			remove(temp.begin(), temp.end(), ' ');
			//cout << temp << endl;
			vecstr.push_back(temp);
			mySuffix = testmatch.suffix();
		}

	}

	return vecstr;
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
	regex testreg_name("[A-Za-z\\_]+");
	
	string wname = "";
	double mass = 0; 
	double mass_inc = 0;
	vector<string> chnames = {};
	vector<double> couplings = {};
	vector<double> couplings_inc = {};

	if(regex_search(cmd, match, reg_AddPole)){
		//for (int i = 0; i<15; i++) cout << i << ": " << match[i] << endl;
		//exit(0);
		wname = match[1];
		remove(wname.begin(), wname.end(), ' ');
		//cout << wname << endl;
		
		string mySuffix, temp;
		mySuffix = match[2];
		int counter = 0;
        while(regex_search(mySuffix, testmatch, testreg_number))
        {
            temp = testmatch[0];
            mySuffix = testmatch.suffix();
			if (counter%2 == 0){
                mass = stod(string(temp));
				//cout << mass << endl;
            }else{
                mass_inc = stod(string(temp));
				//cout << mass_inc << endl;
            }
			counter++;
        }

		mySuffix = match[3];

        counter = 0;
        while(regex_search(mySuffix, testmatch, testreg_name))
        {
            temp = testmatch[0];
			remove(temp.begin(), temp.end(), ' ');
			//cout << temp << endl;
			chnames.push_back(temp);
            mySuffix = testmatch.suffix();
            counter++;
        }
		
		mySuffix = match[6];

		counter = 0;
		while(regex_search(mySuffix, testmatch, testreg_number))
		{
			temp = testmatch[0];
			//cout << temp << endl;
			mySuffix = testmatch.suffix();
			if (counter%2 == 0){
				couplings.push_back(stod(string(temp)));
			}else{
				couplings_inc.push_back(stod(string(temp)));
			}
			
			counter++;
		}

	}

	if(chnames.size()!=couplings.size()) {
		cout << "chnames.size and couplings.size don't match: " << chnames.size() << " vs " << couplings.size() << endl; 
		return poleDat("",0,0,{},{},{});
	}

	//cout << "in filereader::readPole: " << chnames.size() << " vs " << couplings.size() << endl; 

	/*
	cout << wname << " " << mass << " " << mass_inc << endl;
	for(int i = 0; i < chnames.size(); i++){
		cout << chnames[i] << endl;
	}

	for(int i = 0; i < couplings.size(); i++){
		cout << couplings[i] << " " << couplings_inc[i] << endl;
	}
	*/

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
	obsObject.nParams = obsObject.getFitParams().size();
	return obsObject;
}

void filereader::setObs(observable o){
	obsObject = o;
	return;
}

void filereader::loadExpData(){

	ifstream letsread;
	int nWaves = obsObject.getNumAmps();
	int nChans = obsObject.amplitudes[0].getNumOfChans();
	vector<expchan> allDat = {};

	vector<expchan> withDummies = {};

	for(string ampname: obsObject.getAmpNames()){
		for(string chname: obsObject.amplitudes[0].getChanNames()){
			expchan aux = expchan(ampname, chname, {}, {}, {}, {});
			withDummies.push_back(aux);
		}
	}

	for(string s : getExpDataList()){
		//for each command in the input file

		//getting info from the file, chan name, wave name, data name
		expdataDat expd = readExpData(s);
		//read the file
		letsread.open(expd.filename);
		vector<double> a(5);
		vector<double> x = {};
		vector<double> y = {};
		vector<double> y_stat_err = {};
		vector<double> y_sist_err = {};
		while(letsread>>a[0]>>a[1]>>a[2]>>a[3]>>a[4]){
			x.push_back(a[0]);
			y.push_back(a[1]);
			y_stat_err.push_back(a[2]);
			y_sist_err.push_back(a[3]);
		}
		letsread.close();
		
		//make a data object from that file.			
		expchan tempData = expchan(expd.wavename, expd.channame, x, y, y_stat_err, y_sist_err);
		allDat.push_back(tempData);
	}

	//cout << withDummies.size() << " " << allDat.size() << endl;

	for(int i = 0; i < withDummies.size(); i++){
		for(int j = 0; j < allDat.size(); j++){
			//cout << withDummies[i].wavename << " " << allDat[j].wavename << " " << withDummies[i].channame << " " << allDat[j].channame << endl;
			if((withDummies[i].wavename == allDat[j].wavename) && (withDummies[i].channame == allDat[j].channame)){
				withDummies[i] = allDat[j];
			}
		}
	}

	obsObject.setData(withDummies);
	//obsObject.setData(allDat);
	/*for(expchan x: obsObject.getData()){
		if(x.sqrts.size() != 0) cout << x.sqrts[0] << endl;
		else cout << "dummy" << endl;
	}*/
}

void filereader::loadExpInclCrossSec(){

	//cout << ExpInclCrossSecFilename << endl;

	readExpInclCrossSecCmd(ExpInclCrossSecCmd);

	//cout << ExpInclCrossSecFilename << endl;

	ifstream letsread(ExpInclCrossSecFilename);

	vector<double> a(5);
	vector<double> x = {};
	vector<double> y = {};
	vector<double> y_stat_err = {};
	vector<double> y_sist_err = {};
	while(letsread>>a[0]>>a[1]>>a[2]>>a[3]>>a[4]){
		x.push_back(a[0]);
		y.push_back(a[1]);
		y_stat_err.push_back(a[2]);
		y_sist_err.push_back(a[3]);
	}
	letsread.close();

	obsObject.setData_InclCrossSec(expInclCrossSec(x, y, y_stat_err, y_sist_err)); 

}

void filereader::SetExpDataList(){
	regex reg_ExpData("LoadExpData\\(\\s*\"(.*?)\"\\s*,\\s*\"(.*?)\"\\s*,\\s*\"(.*?)\"\\s*\\)");
	smatch cmdmatch;
	for(int i = 0; i < commands.size(); i++){
		if(regex_search(commands.at(i), cmdmatch, reg_ExpData)){
            ExpData_list.push_back(commands[i]);
		}
    }
}


/*void filereader::randomize(){
	vector<double> params = obsObject.getFitParams();
	vector<double> stepsizes = obsObject.getStepSizes();

	TRandom3 gen(getSeed());

	//might need to make sure pole masses are positive somehow 
	for(int i = 0; i < params.size(); i++){
		params[i]= gen.Uniform(params[i] - stepsizes[i], params[i] + stepsizes[i]);
	}

	obsObject.setFitParams(params);

	return;
}*/

void filereader::randomize(int new_seed){
	vector<double> params = obsObject.getFitParams();
	vector<double> stepsizes = obsObject.getStepSizes();
	seed = new_seed;
	TRandom3 gen(new_seed);

	//might need to make sure pole masses are positive somehow 
	for(int i = 0; i < params.size(); i++){
		params[i]= gen.Uniform(params[i] - stepsizes[i], params[i] + stepsizes[i]);
	}

	obsObject.setFitParams(params);

	return;
}

expdataDat filereader::readExpData(string cmd){
	regex reg_ExpData("LoadExpData\\(\\s*\"(.*?)\"\\s*,\\s*\"(.*?)\"\\s*,\\s*\"(.*?)\"\\s*\\)");
	regex testreg_name("[A-Za-z0-9\\.\\-\\_]+");
	string wavename = "";
	string chname = "";
	string fname = "";
	
	if(regex_search(cmd, match, reg_ExpData)){
		//cout << match[0] << endl;
		wavename = match[1];
		remove(wavename.begin(), wavename.end(), ' ');
		chname = match[2];
		remove(chname.begin(), chname.end(), ' ');
		fname = match[3];
		remove(fname.begin(), fname.end(), ' ');
		//cout << fname << endl;
    }

	return expdataDat(wavename, chname, fname);
}

void filereader::readExpInclCrossSecCmd(string cmd){

	regex reg_inclcrosssec("LoadExpInclusiveCrossSection\\(\\s*\"(.*?)\"\\s*\\)");
	
	if(regex_search(cmd, match, reg_inclcrosssec)){
		ExpInclCrossSecFilename = match[1];
    }
}
void filereader::writeOutputFile(string outname){
	
	ofstream letswrite(outname);
	for(string s: getOutputCmds()) letswrite << s << endl;

	letswrite.close();
	 
}


vector<string> filereader::getOutputCmds(){
	vector<string> output_cmds = {};

	output_cmds.push_back("SetSeed("+to_string(getSeed())+")");

	output_cmds.push_back(FitRegion);

	output_cmds.push_back(FitSequenceCmd);

	output_cmds.push_back(fitparamsCmd);

	output_cmds.push_back(InclChi2WeightCmd);

	output_cmds.push_back(ExclChi2WeightCmd);

	output_cmds.push_back(Chi2CutOffCmd);


	output_cmds.push_back(ActionFlagCmd);


	output_cmds.push_back(RandomizeFlagCmd);
	output_cmds.push_back(InclCrossSecFlagCmd);

	output_cmds.push_back(GridCmd);
	output_cmds.push_back(ZeroCmd);
	
	for(string s: AddChannel_list) output_cmds.push_back(s);


	for(string s: AddWave_list) output_cmds.push_back(s);


	vector<string> clist = RewriteChebyList();
	output_cmds.insert(output_cmds.end(),clist.begin(),clist.end());


	vector<string> plist = RewriteAddPoleList();
	output_cmds.insert(output_cmds.end(),plist.begin(),plist.end());


	vector<string> klist = RewriteKmatList();
	output_cmds.insert(output_cmds.end(),klist.begin(),klist.end());


	for(string s: ExpData_list) output_cmds.push_back(s);


	output_cmds.push_back(ExpInclCrossSecCmd);

return output_cmds;
	 
}


void filereader::writeMathematicaOutputFile(string outname){

	ofstream letswrite(outname);

	string temp = "";

	temp = to_string(obsObject.numChans) + " ! number of channels";

	Math_output_cmds.push_back(temp);

	vector<double> masses = {};

	for(int i = 0; i < obsObject.numChans; i++){
		temp = "";
		masses = obsObject.amplitudes[0].getChannels()[i].getMasses();
		for(int j = 0; j < masses.size(); j++){
			temp += to_string(masses[j]) + " ";
		}
		temp += "   ! masses of channel " + to_string(i + 1);
		Math_output_cmds.push_back(temp);
	}

	int num_conform_params = 0;

	for(int i = 0; i < obsObject.numChans; i++){
		temp = "";
		num_conform_params = obsObject.amplitudes[0].getChannels()[i].getChebyCoeffs().size();
		temp += to_string(num_conform_params) + " ! number of conformal parameters in channel " + to_string(i + 1);
		Math_output_cmds.push_back(temp);
	}

	temp = "@@@ ! random conformal (seed)?";//check
	Math_output_cmds.push_back(temp);

	temp = to_string(getKmatList().size()) + " 0 ! Order of K polynomial";
	Math_output_cmds.push_back(temp);

	for(int i = 0; i < obsObject.getNumAmps(); i++){
		temp = to_string(obsObject.amplitudes[i].getResMasses().size()) + " ! number of resonances in wave " + to_string(i + 1);
		Math_output_cmds.push_back(temp);
	}

	vector<double> fitreg = obsObject.amplitudes[0].getFitInterval();

	temp = "";

	for(int i = 0; i < fitreg.size(); i++){
		temp += to_string(sqrt(fitreg[i])) + " ";
	}

	temp += "! Fitted energy region";
	Math_output_cmds.push_back(temp);

	temp = "0 0 0 2 1 !Adler zero";//check
	Math_output_cmds.push_back(temp);

	for(int i = 0; i < obsObject.getNumAmps(); i++){
		temp = to_string(obsObject.amplitudes[i].getSl());
		temp += "                1   0.0000000000000000      ! Scale parameter in CM for wave            " + to_string(i + 1);//check
		Math_output_cmds.push_back(temp);
	}

	double cheb = 0;
	double cheb_step = 0;

	for(int k = 0; k < obsObject.numChans; k++){
		for(int j = 0; j < obsObject.getNumAmps(); j++){
			for(int i = 0; i < num_conform_params; i++){
				//-2173.7290154401403                0   100.00000000000000      ! numerator parameter (channel            1 ,wave            1 ) n.            1
		//cheb_vec = obsObject.amplitudes[amp_index].getChannels()[chan_index].getChebyCoeffs();
		//cheb_vec_steps = obsObject.amplitudes[amp_index].getChannels()[chan_index].getChebySteps();
				cheb = obsObject.amplitudes[j].getChannels()[k].getChebyCoeffs()[i];
				temp = to_string(cheb) + "                ";
				cheb_step = obsObject.amplitudes[j].getChannels()[k].getChebySteps()[i];
				if(cheb_step != 0) temp += "0   ";
				else temp += "1   ";
				temp += to_string(cheb_step) + "      ";
				temp += "! numerator parameter (channel            " + to_string(k + 1) + " ,wave            " + to_string(j + 1) + " ) n.            " + to_string(i + 1);
				Math_output_cmds.push_back(temp);
			}
		}
	}

	// 1.0000000000000000                1   0.0000000000000000      ! scale in conf map (channel            1 ,wave            1 )

	for(int k = 0; k < obsObject.numChans; k++){
		for(int j = 0; j < obsObject.getNumAmps(); j++){
			temp = " 1.0000000000000000                1   0.0000000000000000      ! scale in conf map (channel            " + to_string(k + 1) + " ,wave            " + to_string(j + 1) + " )";
			Math_output_cmds.push_back(temp);
		}
	}

	//-6.0584496703077093                0   1.0000000000000000      ! coupling (channel            1 ,wave            1 , res            1 )

	vector<double> res = obsObject.amplitudes[0].getResMasses();

	double coup = 0;
	double coup_step = 0;	

	for(int j = 0; j < res.size(); j++){
		for(int i = 0; i < obsObject.getNumAmps(); i++){
			for(int k = 0; k < obsObject.numChans; k++){
				coup = obsObject.amplitudes[i].getChannels()[k].getCouplings()[j];
				coup_step = obsObject.amplitudes[i].getChannels()[k].getCouplingSteps()[j];
				temp = to_string(coup) + "                0   " + to_string(coup_step) + "      ! coupling (channel            " + to_string(k + 1) + " ,wave            " + to_string(i + 1) + " , res            " + to_string(j + 1) + " )";
				Math_output_cmds.push_back(temp);
			}
		}
	}

	//6.31602747547477250E-003           0   1.0000000000000000      ! mass (res            1 ,wave            1 )

	double resm = 0;
	double res_step = 0;

	for(int i = 0; i < obsObject.getNumAmps(); i++){
		for(int j = 0; j < res.size(); j++){
			resm = obsObject.amplitudes[i].getResMasses()[j];
			res_step = obsObject.amplitudes[i].getResMassesSteps()[j]; 
			temp = to_string(resm) + "           0   " + to_string(res_step) + "      ! mass (res            " + to_string(j + 1) + " ,wave            " + to_string(i + 1) + " )";
			Math_output_cmds.push_back(temp);
		}
	}

	for(int i = 0; i < 6; i++){
		temp = " 1.0000000000000000                1   0.0000000000000000      ! xn";
		Math_output_cmds.push_back(temp);
	}

	double kcoeff = 0;
	double kcoeff_step = 0;

	//16.190285088168821                0   10.000000000000000      ! c1  1   deg0   iw1
	int count = 0;
	for(int i = 0; i < obsObject.getNumAmps(); i++){
		for(int j = 0; j < getKmatList().size(); j++){

			vector<double> steps = obsObject.amplitudes[i].getKSteps(j);
			MatrixXcd mat_steps = obsObject.amplitudes[i].getkParameters()[j];

			for(int m = 0; m < obsObject.numChans; m++){
				for(int n = 0; n < obsObject.numChans; n++){
					if (m == 0 || n == 0) mat_steps(m,n) = steps[m + n]; 
					else mat_steps(m,n) = steps[m + n + 1];
				}
			}

			for(int l = 0; l < obsObject.numChans; l++){
				for(int k = l; k < obsObject.numChans; k++){
					kcoeff = (obsObject.amplitudes[i].getkParameters()[j](l,k)).real();
					kcoeff_step = (mat_steps(l,k)).real();
					temp = to_string(kcoeff);
					if(kcoeff_step != 0) temp += "                0   ";
					else temp += "                1   ";
					temp += to_string(kcoeff_step) + "      ! c" + to_string(l + 1) + "  1   deg" + to_string(j) + "   iw" + to_string(i + 1);
					Math_output_cmds.push_back(temp);
				}
			}
		}
	}

	temp = "0.0000000000000000                1   0.0000000000000000      ! CM subtraction            1";
	Math_output_cmds.push_back(temp);

	temp = "0.0000000000000000                1   0.0000000000000000      ! CM subtraction            2";	
	Math_output_cmds.push_back(temp);

	temp = "           1 ! CM condition";
	Math_output_cmds.push_back(temp);

	temp = "           1 ! Bootstrap";
	Math_output_cmds.push_back(temp);

	for(string s: Math_output_cmds) letswrite << s << endl;

	letswrite.close();

}


vector<string> filereader::RewriteAddPoleList(){
	vector<string> plist = {};
	regex testreg_number("([+-]{0,1}?(?=\\.\\d|\\d)(?:\\d+)?(?:\\.?\\d*))(?:[Ee]([+-]{0,1}?\\d+))?");
	regex testreg_name("[A-Za-z0-9\\.\\-\\_]+");
	regex reg_AddPole("AddPole\\(\"(.*?)\"\\s*,\\s*([0-9\\.-]+\\s*\\\\pm\\s*[0-9\\.]+)\\s*,\\s*\\{\\s*((\"(.*?)\"\\s*,?\\s*)+)\\}\\s*,\\s*\\{\\s*(([0-9\\.-]+\\s*\\\\pm\\s*[0-9\\.]+\\s*,?\\s*)+)\\}\\s*\\)"); 
	regex reg_AddPole_prefix("AddPole\\(\"(.*?)\"\\s*,\\s*([0-9\\.-]+\\s*\\\\pm\\s*[0-9\\.]+)\\s*,\\s*\\{\\s*((\"(.*?)\"\\s*,?\\s*)+)\\}\\s*,\\s*\\{\\s*");
	regex reg_AddPole_suffix("\\}\\s*\\)");
	regex reg_mass_pm_inc("[0-9\\.-]+\\s*\\\\pm\\s*[0-9\\.]+");

	string wname = "";

	string prefixe = "";
	string suffixe = "";

	int numchans = obsObject.amplitudes[0].getNumOfChans();
	int numpoles = obsObject.amplitudes[0].getResMasses().size();
	double couplings_to_poles[numchans][numpoles];
	double couplings_to_poles_steps[numchans][numpoles];

	int pole_index = 0;

	for(string s : getAddPoleList()){

		if(regex_search(s, match, reg_AddPole)){

			wname = match[1];
			remove(wname.begin(), wname.end(), ' ');
	
		}

		int amp_index = obsObject.getampindex(wname);

		if(regex_search(s, match, reg_AddPole_prefix)) prefixe = match[0];
		if(regex_search(s, match, reg_AddPole_suffix)) suffixe = match[0];

		for(int i = 0; i < numchans; i++){
			for(int j = 0; j < numpoles; j++){
				couplings_to_poles[i][j] = obsObject.amplitudes[amp_index].getChannels()[i].getCouplings()[j];
				couplings_to_poles_steps[i][j] = obsObject.amplitudes[amp_index].getChannels()[i].getCouplingSteps()[j];
			}
		}

		s = prefixe;

		for(int i = 0; i < numchans; i++){
			s += to_string(couplings_to_poles[i][pole_index]) + " \\pm " + to_string(couplings_to_poles_steps[i][pole_index]);
			if(i != numchans - 1) s += ", ";
		}

		s += suffixe;

		if(regex_search(s, testmatch, reg_mass_pm_inc)){

			double mass = obsObject.amplitudes[amp_index].getResMasses()[pole_index];

			double inc_mass = obsObject.amplitudes[amp_index].getResMassesSteps()[pole_index];

			s = string(testmatch.prefix()) + to_string(mass) + " \\pm " + to_string(inc_mass) + string(testmatch.suffix()); 

			plist.push_back(s);

		}

		//cout << s << endl;

		pole_index++;
		
	}
	return plist;
	
}

vector<string> filereader::RewriteChebyList(){
	vector<string> clist = {};
	regex testreg_number("([+-]{0,1}?(?=\\.\\d|\\d)(?:\\d+)?(?:\\.?\\d*))(?:[Ee]([+-]{0,1}?\\d+))?");
	regex testreg_name("[A-Za-z0-9\\.\\-\\_]+");
	regex reg_Cheby("ChebyCoeffs\\(\\s*\"(.*?)\"\\s*,\\s*\"(.*?)\"\\s*,\\s*\"(.*?)\"\\s*,\\s*\\{((\\s*[0-9\\.-]+\\s*\\\\pm\\s*[0-9\\.]+\\s*,?\\s*)+)\\}\\)");
	regex reg_Cheby_prefix("ChebyCoeffs\\(\\s*\"(.*?)\"\\s*,\\s*\"(.*?)\"\\s*,\\s*\"(.*?)\"\\s*,\\s*\\{");
	regex reg_Cheby_suffix("\\}\\)");

	string wname = "";
	string chname = "";

	string prefixe = "";
	string suffixe = "";

	for(string s : getChebyList()){


		if(regex_search(s, match, reg_Cheby)){
			wname = match[1];
			remove(wname.begin(), wname.end(), ' ');
			chname = match[2];
			remove(chname.begin(), chname.end(), ' ');
		}

		if(regex_search(s, match, reg_Cheby_prefix)) prefixe = match[0];
		if(regex_search(s, match, reg_Cheby_suffix)) suffixe = match[0];

		int amp_index = obsObject.getampindex(wname);
		int chan_index = obsObject.getchanindex(wname, chname);

		s = prefixe;

		vector<double> cheb_vec = {};
		vector<double> cheb_vec_steps = {};
		int num_of_chebs = obsObject.amplitudes[0].getChannels()[0].getChebyCoeffs().size();

		cheb_vec = obsObject.amplitudes[amp_index].getChannels()[chan_index].getChebyCoeffs();
		cheb_vec_steps = obsObject.amplitudes[amp_index].getChannels()[chan_index].getChebySteps();

		for(int i = 0; i < num_of_chebs; i++){
			s += to_string(cheb_vec[i]) + " \\pm " + to_string(cheb_vec_steps[i]);
			if(i != num_of_chebs - 1) s += ", ";
		}

		s += suffixe;

		clist.push_back(s);

	}
	return clist;
}

vector<string> filereader::RewriteKmatList(){
	vector<string> klist = {};
	regex reg_Kmat("AddKmatBackground\\(\\s*\"(.*?)\"\\s*,\\s*([0-9\\.-]+)\\s*,\\s*\\{((\\s*\\{((\\s*[0-9\\.-]+\\s*\\\\pm\\s*[0-9\\.]+\\s*,?\\s*)+)\\}\\s*,?\\s*)+)\\}\\s*\\)");
	regex reg_Kmat_core("\\{((\\s*\\{((\\s*[0-9\\.-]+\\s*\\\\pm\\s*[0-9\\.]+\\s*,?\\s*)+)\\}\\s*,?\\s*)+)\\}");

	string wname = {};
	int power = 0;

	for(string s : getKmatList()){

		if(regex_search(s, match, reg_Kmat)){

			wname = match[1];
			remove(wname.begin(), wname.end(), ' ');
			power = stoi(string(match[2]));

		}

		int amp_index = obsObject.getampindex(wname);

		if(regex_search(s, testmatch, reg_Kmat_core)){

			MatrixXcd mat = obsObject.amplitudes[amp_index].getkParameters()[power];
			vector<double> steps = obsObject.amplitudes[amp_index].getKSteps(power);
			MatrixXcd mat_steps = mat;

			int numchan = obsObject.amplitudes[amp_index].getNumOfChans();

			for(int i = 0; i < numchan; i++){
				for(int j = 0; j < numchan; j++){
					if (i == 0 || j == 0) mat_steps(i,j) = steps[i + j]; 
					else mat_steps(i,j) = steps[i + j + 1];
				}
			}

			string prefix = testmatch.prefix();
			string suffix = testmatch.suffix();

			s = prefix + "{";

			for(int i = 0; i < numchan; i++){
				s += "{" ;
				for(int j = i; j < numchan; j++){
					s += to_string(mat(i,j).real()) + " \\pm " + to_string(mat_steps(i,j).real());
					if (j != numchan - 1) s += ", ";
				}
				if (i != numchan - 1) s += "}, ";
				else s += "}";
			}

			s += "}" + suffix;

		}


		klist.push_back(s);

	}
	return klist;
}
