#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <deque>
#include <cmath>
#include <algorithm>
#include <getopt.h>

#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"

#include "TrackManager.h"
#include "TrackFunctions.h"

void Usage(char* Name);
int main(int argc, char** argv)
{
	const struct option longopts[] = 
	{
		{"intput", 	required_argument, 	0, 'i'},
		{"spectrum", 	required_argument, 	0, 'o'},
		{"peakvall", 	required_argument, 	0, 'l'},
		{"rootout", 	required_argument, 	0, 'r'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	std::ifstream InFile;
	std::ofstream OutFile;
	TFile *RootFile;
	unsigned int Integ = 0;
	std::string ListFile, BaseOut;
	unsigned int A = 0, B = 1000;
	
	while((iarg = getopt_long(argc,argv, "I:f:s:o:l:A:B:h", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'I':	//smoothening window
				Integ = strtol(optarg, NULL, 10);
				break;
			case 'f':
				ListFile.assign(optarg);
				break;
			case 's':
				BaseOut.assign(optarg);
				break;
			case 'o':
				OutFile.open(optarg);
				break;
			case 'r':
				RootFile = new TFile(optarg, "RECREATE");
				break;
			case 'A':
				A = std::strtol(optarg, NULL, 10);
				break;
			case 'B':
				B = std::strtol(optarg, NULL, 10);
				break;
			case 'h':
				Usage(argv[0]);
				return 1;
			default:
				break;
		}
	
	}

	std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;

	TrackManager   *Manager   = new TrackManager(ListFile);
	if (Manager->GetEntries() < 2)
	{
		std::cerr << "Not enough entires" << std::endl;
		return 1;
	}

	TrackFunctions *Functions = new TrackFunctions();

	std::vector<double> v0;
	Manager->AverageTrack(v0, 0, 1);		//v0 contains the error
	bool Err = (2*Manager->GetX() == v0.size());	//means error is in there, so must be propagated
	//Functions->Smoothen(v0, 1);
	Functions->Baseline(v0, 1, 1);
	Functions->NormaliseArea(v0, A, A+40, 1);

	for (unsigned int i = 1; i < Manager->GetEntries(); ++i)
	{
		std::vector<double> vY, vAbs;

		if (Manager->AverageTrack(vY, i, 1))
		{
			//Functions->Smoothen(vY, 1);
			Functions->Baseline(vY, 1, 1);
			Functions->NormaliseArea(vY, A, A+40, 1);

			vAbs.insert(vAbs.end(), vY.begin(), vY.end());
			Functions->AbsorptionVar(v0, vAbs, 1);

			std::stringstream ssL;
			ssL << BaseOut << Manager->GetPercentage(i) << ".dat";

			std::ofstream TrackOut(ssL.str().c_str());
			for (unsigned int j = 0; j < vY.size(); ++j)
			{
				TrackOut << Manager->GetX(j) << "\t" << v0.at(j) << "\t" << sqrt(v0Var.at(j)) << "\t";
				TrackOut << vY.at(j) << "\t" << sqrt(vYVar.at(j)) << "\t";
				TrackOut << vAbs.at(j) << "\t" << sqrt(vAbs.at(j)) << std::endl;
			}

			std::vector<unsigned int> iP, iV;
			Functions->FindPeakValley(vAbs, iP, iV, A, B, 1);

			/*
			std::cout << "Peaks" << std::endl;
			for (unsigned int j = 0; j < iP.size(); ++j)
				std::cout << Manager->GetX(iP.at(j)) << "\t" << vAbs.at(iP.at(j)) << std::endl;
				//std::cout << iP.at(j) << "\t" << vAbs.at(iP.at(j)) << std::endl;
			std::cout << "Deeps" << std::endl;
			for (unsigned int j = 0; j < iV.size(); ++j)
				std::cout << Manager->GetX(iV.at(j)) << "\t" << vAbs.at(iV.at(j)) << std::endl;
			*/

			double AbsDiff    = Functions->AbsorptionDiff(vY, iP.at(0), iP.at(1));
			double AbsDiffVar = Functions->AbsorptionDiffVar(vY, vYVar, iP.at(0), iP.at(1));
			Out << Manager->GetPercentage(i) << "\t";
			Out << vAbs.at(iP.at(0)) << "\t" << sqrt(vAbsVar.at(iP.at(0))) << "\t";
			Out << vAbs.at(iP.at(1)) << "\t" << sqrt(vAbsVar.at(iP.at(1))) << "\t";
			Out << AbsDiff << "\t" << sqrt(AbsDiffVar) << std::endl;

		}
	}

	return 0;
}

void Usage(char *Name)
{
	std::cout << "Description" << std::endl;
	std::cout << "Usage : " << Name << " [OPTIONS]" << std::endl;
	std::cout << std::endl;
	std::cout <<"\n  -i,  --input\t\t";
	std::cout << "Input test file containing the spectrum (wavelength, intensity)" << std::endl;
	std::cout <<"\n  -o,  --output\t\t";
	std::cout << "Output file for testing now" << std::endl;
	std::cout <<"\n  -r,  --rootout\t\t";
	std::cout << "Output ROOT file with the absorption spectrum and other stuff" << std::endl;
	std::cout <<"\n  -I,  --intwindow\t\t";
	std::cout << "Define width of integration window" << std::endl;
	std::cout <<"\n  -h,  --help\t\t";
	std::cout << "Print this message and exit" << std::endl;
}
