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

class Sort
{
	private:
		std::vector<double> vX;
	public:
		Sort(std::vector<double> &vv) : vX(vv) {}
		bool operator()(int i, int j) const { return vX.at(i) < vX.at(j); }
};

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
	std::ofstream SpecFile, PeakFile;
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
			case 'l':
				PeakFile.open(optarg);
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

	std::ostream &OutS = (SpecFile.is_open()) ? SpecFile : std::cout;
	std::ostream &OutP = (PeakFile.is_open()) ? PeakFile : std::cout;

	TrackManager   *Manager   = new TrackManager(ListFile);
	std::cout << "H0" << std::endl;
	if (Manager->GetEntries() < 2)
	{
		std::cerr << "Not enough entires" << std::endl;
		return 1;
	}

	std::cout << "H1" << std::endl;
	TrackFunctions *Functions = new TrackFunctions();

	std::cout << "H2" << std::endl;
	std::vector<double> v0;
	std::cout << "H3" << std::endl;
	Manager->AverageTrack(v0, 0);
	std::cout << "H4" << std::endl;
	Functions->Normalise(v0, A, B);
	std::cout << "H5" << std::endl;

	std::cout << "H6 " << Manager->GetEntries() <<  std::endl;
	for (unsigned int i = 1; i < Manager->GetEntries(); ++i)
	{
		std::vector<double> vY, vAbs;

		std::cout << "F0" << std::endl;
		if (Manager->AverageTrack(vY, i))
		{
			std::cout << "I0" << std::endl;
			vAbs.insert(vAbs.end(), vY.begin(), vY.end());
			std::cout << "I1" << std::endl;
			Functions->Normalise(vAbs, A, B);
			std::cout << "I2" << std::endl;
			Functions->Absorption(v0, vAbs);

			std::cout << "I3" << std::endl;
			std::stringstream ssL(BaseOut);
			std::cout << "I4" << std::endl;
			ssL << Manager->GetPercentage(i) << ".dat";

			std::cout << "I5" << std::endl;
			std::ofstream Out(ssL.str().c_str());
			std::cout << "I6" << std::endl;
			for (unsigned int j = 0; j < vY.size(); ++j)
				Out << Manager->GetX(j) << "\t" << vY.at(j) << "\t" << vAbs.at(j) << std::endl;
			std::cout << "I0" << std::endl;

		}
		std::cout << "F1" << std::endl;
	}
	std::cout << "H7 " << Manager->GetEntries() <<  std::endl;

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
