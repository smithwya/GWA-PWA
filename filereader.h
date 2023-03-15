#pragma once
#include <iostream>
#include<vector>
#include<string>
#include "amplitude.h"
#include "channel.h"



using namespace std;

class filereader {

public:
	filereader(string fname);
	void readFile(string fname);
	vector<amplitude> constructAmps();
	void printCommands();

private:
	vector<std::string> commands;

};