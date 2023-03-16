#pragma once
#include <iostream>
#include<vector>
#include<string>
#include "amplitude.h"
#include "channel.h"
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
	vector<double> coeffs;

	chebyDat(string wn,string cn,string p, vector<double> c){
		wavename = wn;
		channame = cn;
		poletype = p;
		coeffs = c;
	}


};

struct poleDat{
	string wavename;
	double mass;
	vector<string> channames;
	vector<double> couplings;

	poleDat(string wn, double m, vector<string> cn, vector<double> cs){

		wavename = wn;
		mass = m;
		channames = cn;
		couplings =cs;
	}


};

struct kmatDat{
	string wavename;
	int power;
	vector<double> matelems;

	kmatDat(string wn, int p, vector<double> m){
		wavename = wn;
		power = p;
		matelems = m;
	}


};

class filereader {

public:
	filereader(string fname);
	void readFile(string fname);
	vector<amplitude> constructAmps();
	void printCommands();
	string getCommand(int i);
	int readSeed(string cmd);
	vector<double> readFitReg(string cmd);
	chanDat readCh(string cmd);
	ampDat readWave(string cmd);
	chebyDat readCheby(string cmd);
	poleDat readPole(string cmd);
	kmatDat readKmat(string cmd);




private:
	vector<std::string> commands;
	smatch match;
};