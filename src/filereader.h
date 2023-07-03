#pragma once
#include <iostream>
#include<vector>
#include<string>
#include "amplitude.h"
#include "channel.h"
#include "observable.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>


using namespace std;

struct chanDat{
	string name;
	vector<double> ch_masses;

	chanDat(string n, vector<double> m){
		name = n;
		ch_masses = m;
	}

};


struct ampDat{

	string wname;
	string kmatflag;
	string rhoN;
	int J;
	double sL;

	ampDat(string wn, string km, string rn, int Jj, double ssl){
		wname = wn;
		kmatflag = km;
		rhoN = rn;
		J = Jj;
		sL = ssl; 
	}


};

struct chebyDat{
	string wavename;
	string channame;
	string poletype;
	double ss0;
	vector<double> coeffs;
	vector<double> coeffs_inc;

	chebyDat(string wn,string cn,string p, double S0, vector<double> c, vector<double> c_inc){
		wavename = wn;
		channame = cn;
		poletype = p;
		ss0 = S0;
		coeffs = c;
		coeffs_inc = c_inc;
	}


};

struct poleDat{
	string wavename;
	double mass, mass_inc;
	vector<string> channames;
	vector<double> couplings, couplings_inc;

	poleDat(string wn, double m, double m_inc, vector<string> cn, vector<double> cs, vector<double> cs_inc){

		wavename = wn;
		mass = m;
		mass_inc = m_inc;
		channames = cn;
		couplings =cs;
		couplings_inc = cs_inc;
	}


};

struct kmatDat{
	string wavename;
	int power;
	vector<double> matelems, matelems_inc;

	kmatDat(string wn, int p, vector<double> m, vector<double> m_inc){
		wavename = wn;
		power = p;
		matelems = m;
		matelems_inc = m_inc;
	}


};

struct expdataDat{

	string wavename;
	string channame;
	string filename;

	expdataDat(string wn, string cn, string fn){
		wavename = wn;
		channame = cn;
		filename = fn;
	}


};

class filereader {

public:
	filereader(string fname);
	void readFile(string fname);
	void ConstructBareAmps();
	//void CompleteBareAmps(vector<amplitude> amplist);
	void printCommands();
	string getCommand(int i);
	void SetChi2CutOff();
	void SetFitFlag();
	void SetRandomizeFlag();
	void SetSeed();
	void SetFitRegion();
	void SetInclCrossSecFlag();
	void SetAddChannelList();
	void SetAddWaveList();
	void SetChebyList();
	void SetAddPoleList();
	void SetKmatList();
	void SetAllCommandLists();
	int getSeed();
	bool getFitFlag();
	bool getInclCrossSecFlag();
	bool getRandomizeFlag();
	double getChi2CutOff();
	string getFitRegion();
	vector<string> getAddChannelList();
	vector<string> getAddWaveList();
	vector<string> getChebyList();
	vector<string> getAddPoleList();
	vector<string> getKmatList();
	double readChi2CutOffCmd(string cmd);
	void readFitFlag(string cmd);
	void readInclCrossSecFlag(string cmd);
	void readRandomizeFlag(string cmd);
	int readSeed(string cmd);
	vector<double> readFitReg(string cmd);
	chanDat readCh(string cmd);
	ampDat readWave(string cmd);
	chebyDat readCheby(string cmd);
	poleDat readPole(string cmd);
	kmatDat readKmat(string cmd);
	observable getObs();
	void setObs(observable o);
	void setChebys();
	void setKmats();
	void setPoles();
	void loadExpData();
	void loadExpInclCrossSec();
	expdataDat readExpData(string cmd);
	void readExpInclCrossSecCmd(string cmd);
	void SetExpDataList();
	vector<string> getExpDataList();
	void randomize();
	void randomize(int newseed);
	void writeOutputFile(string outfile);
	void RewriteChebyList();
	void RewriteAddPoleList();
	void RewriteKmatList();
	void setSeed(int newseed);
	void setExpInclCrossSec();

private:
	bool RandomizeFlag = true;
	bool FitFlag = true;
	bool InclCrossSecFlag = true; 
	string NameOfFile;
	vector<std::string> commands = {}, output_cmds = {};
	smatch match, testmatch;
	string FitFlagCmd, SeedCmd, FitRegion, Chi2CutOffCmd, RandomizeFlagCmd, ExpInclCrossSecCmd, ExpInclCrossSecFilename, InclCrossSecFlagCmd;
	vector<string> AddChannel_list, AddWave_list, Cheby_list, AddPole_list, Kmat_list,ExpData_list;
	observable obsObject;
	int seed;
};