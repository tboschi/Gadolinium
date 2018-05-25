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
	Functions->Baseline(v0, 1, Err);
	Functions->NormaliseArea(v0, A, A+40, Err);

	for (unsigned int j = 1; j < Manager->GetEntries(); ++j)
	{
		std::vector<double> vY, vAbs;

		if (Manager->AverageTrack(vY, j, Err))
		{
			//Functions->Smoothen(vY, 1);
			Functions->Baseline(vY, 1, Err);
			Functions->NormaliseArea(vY, A, A+40, Err);

			//vAbs.insert(vAbs.end(), vY.begin(), vY.end());
			vAbs = Functions->Absorption(v0, vY, Err);

			std::stringstream ssL;
			ssL << BaseOut << Manager->GetPercentage(j) << ".dat";

			std::ofstream TrackOut(ssL.str().c_str());
			for (unsigned int i = 0, e = v0.size()/2; i < v0.size()/(1+Err); ++i, ++e)
			{
				TrackOut << Manager->GetX(i) << "\t" << v0.at(i) << "\t" << sqrt(v0.at(e)) << "\t";
				TrackOut << vY.at(i) << "\t" << sqrt(vY.at(e)) << "\t";
				TrackOut << vAbs.at(i) << "\t" << sqrt(vAbs.at(e)) << std::endl;
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
			double AbsDiffVar = -1.0;
			if (Err) 
				AbsDiffVar = Functions->AbsorptionDiffErr(vY, iP.at(0), iP.at(1));

			Out << Manager->GetPercentage(j) << "\t";
			Out << vAbs.at(iP.at(0)) << "\t";
			if (Err)
				Out << sqrt(vAbs.at(iP.at(0)+vAbs.size()/2)) << "\t";
			Out << vAbs.at(iP.at(1)) << "\t";
			if (Err)
				Out << sqrt(vAbs.at(iP.at(1)+vAbs.size()/2)) << "\t";
			Out << AbsDiff;
			if (Err)
				Out << "\t" << sqrt(AbsDiffVar);
			Out << std::endl;
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
