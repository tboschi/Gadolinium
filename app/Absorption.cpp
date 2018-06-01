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
	double A = 0, B = 1000;
	double Volume = 30, Length = 5;	//cm
	
	while((iarg = getopt_long(argc,argv, "I:f:s:o:l:A:B:V:L:h", longopts, &index)) != -1)
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
				A = std::strtod(optarg, NULL);
				break;
			case 'B':
				B = std::strtod(optarg, NULL);
				break;
			case 'V':
				Volume = std::strtod(optarg, NULL);
				break;
			case 'L':
				Length = std::strtod(optarg, NULL);
				break;
			case 'h':
				Usage(argv[0]);
				return 1;
			default:
				break;
		}
	
	}

	std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;

	TrackManager   *Manager   = new TrackManager(ListFile, Volume);
	if (Manager->GetEntries() < 2)
	{
		std::cerr << "Not enough entires" << std::endl;
		return 1;
	}

	TrackFunctions *Functions = new TrackFunctions(0.05);
	Functions->SetInterval(Manager, A, B);

	double Err = true, Subtract = true;
	std::vector<double> v0;
	if (!Manager->AverageTrack(v0, 0, Err))		//v0 contains the error
		std::cout << "I0 is corrupted" << std::endl;

	//Functions->Smoothen(v0, 1);
	//std::cout << "baseline " << Functions->Baseline(v0, Subtract, Err) << std::endl;
	//Functions->NormaliseArea(v0, Err);	//normalise using terzile

	unsigned int R = 0;
	for (unsigned int j = 1; j < Manager->GetEntries(); ++j)
	{
		std::vector<double> vY, vAbs;

		bool Set = Manager->AverageTrack(vY, j, Err);	//if True, one data set
								//if False, more data set, must reloop over folder j
		//Functions->Smoothen(vY, 1);
	//std::cout << "baseline " << Functions->Baseline(v0, Subtract, Err) << std::endl;
		//Functions->NormaliseArea(vY, Err);

		//vAbs.insert(vAbs.end(), vY.begin(), vY.end());
		vAbs = Functions->Absorption(v0, vY, Err);
		//Functions->Baseline

		std::stringstream ssL;
		ssL << BaseOut << "_track_" << std::fixed << std::setprecision(4) << Manager->GetPercentage(j) << "_" << R << ".dat";
		std::ofstream TrackOut(ssL.str().c_str());
		for (unsigned int i = 0, e = v0.size()/2; i < v0.size()/(1+Err); ++i, ++e)
		{
			if (Manager->GetX(i) >= A && Manager->GetX(i) <= B)
			{
				TrackOut << Manager->GetX(i) << "\t" << v0.at(i) << "\t" << sqrt(v0.at(e)) << "\t";
				TrackOut << vY.at(i) << "\t" << sqrt(vY.at(e)) << "\t";
				TrackOut << vAbs.at(i)/Length << "\t" << sqrt(vAbs.at(e))/Length << std::endl;
			}
		}
		TrackOut.close();

		std::vector<unsigned int> iP, iV;
		Functions->FindPeakValley(vAbs, iP, iV, 1);

		ssL.clear();
		ssL.str("");
		ssL << BaseOut << "_peaks_" << std::fixed << std::setprecision(4) << Manager->GetPercentage(j) << "_" << R << ".dat";
		std::ofstream PeaksOut(ssL.str().c_str());

		PeaksOut << "#Peaks" << std::endl;
		for (unsigned int j = 0; j < iP.size(); ++j)
			PeaksOut << Manager->GetX(iP.at(j)) << "\t" << vY.at(iP.at(j)) << "\t" << vAbs.at(iP.at(j))/Length << std::endl;
		PeaksOut << "\n\n#Valleys" << std::endl;
		for (unsigned int j = 0; j < iV.size(); ++j)
			PeaksOut << Manager->GetX(iV.at(j)) << "\t" << vY.at(iV.at(j)) << "\t" << vAbs.at(iV.at(j))/Length << std::endl;
		PeaksOut.close();

		if (iP.size() < 2)
			continue;
		double Abs0    = Functions->AbsorptionDiff(v0, iP.at(0), iP.at(1));
		double Abs0Var = Functions->AbsorptionDiffErr(v0, iP.at(0), iP.at(1));
		double AbsY    = Functions->AbsorptionDiff(vY, iP.at(0), iP.at(1));
		double AbsYVar = Functions->AbsorptionDiffErr(vY, iP.at(0), iP.at(1));

		Out << Manager->GetPercentage(j) << "\t";
		Out << vAbs.at(iP.at(0))/Length << "\t" << sqrt(vAbs.at(iP.at(0)+vAbs.size()/2))/Length << "\t";
		Out << vAbs.at(iP.at(1))/Length << "\t" << sqrt(vAbs.at(iP.at(1)+vAbs.size()/2))/Length << "\t";
		Out << (Abs0-AbsY)/Length << "\t" << sqrt(Abs0Var + AbsYVar)/Length << std::endl;

		if (!Set) 		//repeat, for multiple data set in same folder
		{
			--j;
			++R;
		}
		else
			R = 0;

		/*
		std::cout << "Peaks" << std::endl;
		for (unsigned int j = 0; j < iP.size(); ++j)
			std::cout << Manager->GetX(iP.at(j)) << "\t" << vAbs.at(iP.at(j)) << std::endl;
			//std::cout << iP.at(j) << "\t" << vAbs.at(iP.at(j)) << std::endl;
		std::cout << "Deeps" << std::endl;
		for (unsigned int j = 0; j < iV.size(); ++j)
			std::cout << Manager->GetX(iV.at(j)) << "\t" << vAbs.at(iV.at(j)) << std::endl;
		*/

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
