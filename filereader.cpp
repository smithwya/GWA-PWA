#include "filereader.h"
#include<complex>
#include<algorithm>
#include<math.h>
#include<chrono>
#include<random>
#include<fstream>
#include<vector>
#include<string>

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


	return {};
}

void filereader::printCommands(){
	for(string s : commands){

		cout<<s<<endl;
	}

}