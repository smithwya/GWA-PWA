#pragma once
#include <iostream>
#include<vector>
#include<string>
#include "amplitude.h"
#include "channel.h"
using namespace std;


class observable {


private:
	vector<channel> channels;
	vector<amplitude> amplitudes;

public:
	observable();
	observable(vector<channel> c, vector<amplitude> a);
	double getIntensity();
	double getPhase();

	void addChannel(channel c);
	void addAmp(amplitude a);
	void getAmps();

};