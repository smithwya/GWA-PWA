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
	int counter = 0;
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

void filereader::SetSeed(){
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
	SetSeed();
	SetFitRegion();
	SetAddChannelList();
	SetAddWaveList();
	SetChebyList();
	SetAddPoleList();
	SetKmatList();
	SetExpDataList();
}

int filereader::getSeed(){
	return readSeed(SeedCmd);
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

vector<string> filereader::getExpDataList(){
	return ExpData_list;
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
			//if((withDummies[i].wavename == allDat[j].wavename) && (withDummies[i].channame == allDat[j].channame)) cout << "hello" << endl;
		}
	}//it has to be fixed

	//obsObject.setData(withDummies);
	obsObject.setData(allDat);
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


void filereader::randomize(){
	vector<double> params = obsObject.getFitParams();
	vector<double> stepsizes = obsObject.getStepSizes();

	TRandom3 gen(getSeed());

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
		wavename = match[1];
		remove(wavename.begin(), wavename.end(), ' ');
		chname = match[2];
		remove(chname.begin(), chname.end(), ' ');
		fname = match[3];
		remove(fname.begin(), fname.end(), ' ');

    }

	return expdataDat(wavename, chname, fname);
}

void filereader::writeOutputFile(){

	ofstream letswrite("Data/output_of_" + NameOfFile);

	output_cmds.push_back(SeedCmd);

	output_cmds.push_back("");

	output_cmds.push_back(FitRegion);

	output_cmds.push_back("");

	for(string s: AddWave_list) output_cmds.push_back(s);

	output_cmds.push_back("");

	RewriteChebyList();

	output_cmds.push_back("");

	RewriteAddPoleList();

	output_cmds.push_back("");

	RewriteKmatList();

	output_cmds.push_back("");

	for(string s: ExpData_list) output_cmds.push_back(s);

	for(string s: output_cmds) letswrite << s << endl;

	for(string s: output_cmds) cout << s << endl;
	 
}

void filereader::RewriteAddPoleList(){
	
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

		//cout << s << endl;

		if(regex_search(s, match, reg_AddPole)){

			//for (int i = 0; i<15; i++) cout << i << ": " << match[i] << endl;
			//exit(0);
			wname = match[1];
			remove(wname.begin(), wname.end(), ' ');
			//cout << wname << endl;
	
		}

		int amp_index = obsObject.getampindex(wname);

		if(regex_search(s, match, reg_AddPole_prefix)) prefixe = match[0];
		if(regex_search(s, match, reg_AddPole_suffix)) suffixe = match[0];

		for(int i = 0; i < numchans; i++){
			for(int j = 0; j < numpoles; j++){
				couplings_to_poles[i][j] = obsObject.amplitudes[amp_index].getChannels()[i].getCouplings()[j];
				//cout << couplings_to_poles[i][j] << " ";
				couplings_to_poles_steps[i][j] = obsObject.amplitudes[amp_index].getChannels()[i].getCouplingSteps()[j];
			}
			//cout << endl;
		}

		/*for(string channame: obsObject.amplitudes[amp_index].getChanNames()){
			int chan_index = obsObject.getchanindex(wname, channame);
			for(double coup: obsObject.amplitudes[amp_index].getChannels()[chan_index].getCouplings()){
				cout << coup << " ";
			}
			cout << endl;
		}*/

		//s = regex_replace(s, testreg_number, obsObject.amplitudes[amp_index].)

		s = prefixe;

		for(int i = 0; i < numchans; i++){
			s += to_string(couplings_to_poles[i][pole_index]) + " \\pm " + to_string(couplings_to_poles_steps[i][pole_index]);
			if(i != numchans - 1) s += ", ";
		}

		s += suffixe;

		if(regex_search(s, testmatch, reg_mass_pm_inc)){
			
			//cout << testmatch[0] << endl;

			//cout << testmatch.prefix() << endl << testmatch.suffix() << endl;

			double mass = obsObject.amplitudes[amp_index].getResMasses()[pole_index];

			double inc_mass = obsObject.amplitudes[amp_index].getResMassesSteps()[pole_index];

			s = string(testmatch.prefix()) + to_string(mass) + " \\pm " + to_string(inc_mass) + string(testmatch.suffix()); 

			output_cmds.push_back(s);

		}

		//cout << s << endl;

		pole_index++;
		
	}
	
}

void filereader::RewriteChebyList(){

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

		//cout << s << endl;

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

		//cout << s << endl;

		output_cmds.push_back(s);

	}

}

void filereader::RewriteKmatList(){
	
	regex reg_Kmat("AddKmatBackground\\(\\s*\"(.*?)\"\\s*,\\s*([0-9\\.-]+)\\s*,\\s*\\{((\\s*\\{((\\s*[0-9\\.-]+\\s*\\\\pm\\s*[0-9\\.]+\\s*,?\\s*)+)\\}\\s*,?\\s*)+)\\}\\s*\\)");
	regex reg_Kmat_core("\\{((\\s*\\{((\\s*[0-9\\.-]+\\s*\\\\pm\\s*[0-9\\.]+\\s*,?\\s*)+)\\}\\s*,?\\s*)+)\\}");

	string wname = {};
	int power = 0;

	for(string s : getKmatList()){

		//cout << s << endl;

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

			/*for(int i = 0; i < numchan; i++){
				for(int j = 0; j < numchan; j++){
					cout << mat_steps(i,j) << " ";
				}
				cout << endl;
			}

			for(double x: steps) cout << x << endl;*/
			

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

		//cout << s << endl;

		output_cmds.push_back(s);

	}

}