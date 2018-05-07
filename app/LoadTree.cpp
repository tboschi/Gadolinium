#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <getopt.h>

#include "TFile.h"
#include "TH1D.h"

void Usage(char* Name);
std::vector<double> Markov(const std::vector<double> &Spectrum, int Window);
int main(int argc, char** argv)
{
	const struct option longopts[] = 
	{
		{"intput", 	required_argument, 	0, 'i'},
		{"output", 	required_argument, 	0, 'o'},
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
	
	while((iarg = getopt_long(argc,argv, "I:i:r:o:h", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'I':
				Integ = strtol(optarg, NULL, 10);
				break;
			case 'i':
				InFile.open(optarg);
				break;
			case 'o':
				OutFile.open(optarg);
				break;
			case 'r':
				RootFile = new TFile(optarg, "RECREATE");
				break;
			case 'h':
				Usage(argv[0]);
				return 1;
			default:
				break;
		}
	
	}

	std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;

	double Wavelength, Intensity;
	std::vector<double> vWave, vInty;

	std::string Line;
	std::stringstream ssL;
	while (std::getline(InFile, Line))
	{
		if (Line[0] == '#') continue;

		ssL.clear();
		ssL.str("");

		ssL << Line;
		ssL >> Wavelength >> Intensity;

		vWave.push_back(Wavelength);
		vInty.push_back(Intensity);
	}
	//					4 points per 100nm = 400 bins
	TH1D* hSpectrum = new TH1D("hspectrum", "Spectrum", 400, 300, 400);
	TH1D* hDerivate = new TH1D("hderivate", "Derivate", 400, 300, 400);

	
	double mean = 0, mean2 = 0, rms = 0;
	int sample = 200;
	int start = 0;
	for (unsigned int i = start; i < start+sample; ++i)
	{
		mean += vInty.at(i)/sample;
		rms += pow(vInty.at(i), 2)/sample;
	}
	rms = sqrt(rms);
	for (unsigned int i = start; i < start+sample; ++i)
		mean2 += pow(vInty.at(i)-mean, 2)/(sample-1);
	mean2 = sqrt(mean2);

	unsigned int Window = Integ;

	std::vector<double> vMarkov = Markov(vInty, Integ/2);

	std::vector<double> vCopy(vInty);
	for (int rep = 1; rep < Window; ++rep)
	{
		for (int i = rep; i < vCopy.size()-rep; ++i)
		{
			vCopy.at(i) = (vCopy.at(i+rep) + vCopy.at(i-rep))/2.0;
		}
	}

					
	/*
	std::vector<double> vSmooth = Markov(vCopy, Integ/2);
	for (unsigned int i = 0; i < vInty.size(); ++i)
		Out << vWave.at(i) << "\t" << vInty.at(i) << "\t" << vMarkov.at(i) << "\t" << vCopy.at(i) << "\t" << vSmooth.at(i) << std::endl;
	*/

	//vCopy is almost backgroundless
	//must find a way to subtract baseline autmatically
	//find the lowest and highest point of Copy
	double Min = 999, Max = -999;
	unsigned int iMin, iMax;
	for (unsigned int i = 0; i < vCopy.size(); ++i)
	{
		if (Min > vCopy.at(i))
		{
			Min = vCopy.at(i);
			iMin = i;
		}
		if (Max < vCopy.at(i))
		{
			Max = vCopy.at(i);
			iMax = i;
		}
	}

	double Perc = 0.02;	//2% of the peak is start
	double Thr = Perc*(Max-Min);
	std::cout << "Min " << Min << "\tMax " << Max << "\tDiff " << Max - Min << "\tPerc " << Thr << std::endl;
	double Baseline = 0, Rms = 0;
	int SampleLength = 0;
	int iA, iB;
	for (unsigned int i = 0; i < vCopy.size(); ++i)
	{
		double s0 = vCopy.at(i);
		double sum = 0, sum2 = 0;
		unsigned int j;
		for (j = i; j < vCopy.size(); ++j)
		{
			sum  += vCopy.at(j);
			if (fabs(vCopy.at(j) - s0) > Perc*(Max-Min))
				break;
		}
		if (j-i > SampleLength)
		{
			Baseline = sum/(j-i);
			SampleLength = j-i;
			iA = i, iB = j;
		}
		i = j;		//skip forward
	}
	double Std = 0;
	for (unsigned int i = iA; i < iB; ++i)
		Std += pow(vCopy.at(i) - Baseline, 2)/(iB-iA-1);
	Std = sqrt(Std);
	std::cout << "Sample " << SampleLength << "\tBaseline " << Baseline << std::endl;

	for (unsigned int i = 0; i < vCopy.size(); ++i)
		vCopy.at(i) -= Baseline;

	std::vector<double> vPeak, vVall;
	std::vector<int> iPeak, iVall;
	vPeak.push_back(vCopy.at(iMax));
	iPeak.push_back(iMax);

	//going left first, and then right
	for (int Dir = -1; Dir < 2; Dir += 2)
	{
		double PoV = false;		//F looking for Valley, T looking for P
		int iD = iPeak.front(), iS = iPeak.front();
		double fS = vPeak.front();
		std::cout << "dir " << Dir << "\t" << vWave.at(iS) << "\t" << fS << std::endl;
		while (iD > -1 && iD < vCopy.size())
		{
			int Sign = 2*PoV - 1;	//-1 looking for Valley, +1 looking for Peak
			if ( Sign*(vCopy.at(iD) - fS) > Thr/2.0)
			{
				iS = iD ;
				fS = vCopy.at(iD);
			}
			if (-Sign*(vCopy.at(iD) - fS) > Thr/2.0)
			{
				double fX = -Sign*vPeak.front();
				int    iX = -1;
				//for (unsigned j = iD; Dir*(iS-j) < 0 && j > -1 && j < vCopy.size(); j -= Dir)
				for (int j = iD; Dir*(iS-j) < 0; j -= Dir)
				{
					if (Sign*(vCopy.at(j)-fX) > 0)
					{
						fX = vCopy.at(j);
						iX = j;
					}
				}
				if (iX > -1)
				{
					std::vector<double> &vRef = PoV ? vPeak : vVall;
					std::vector<int> &iRef = PoV ? iPeak : iVall;
					vRef.push_back(fX);
					iRef.push_back(iX);
					fS = fX;
					iD = iS = iX;
					PoV = !PoV;
				}
			}
			iD += Dir;
		}
	}

	std::cout << "Peaks\t" << vPeak.size() << "\tValleys " << vVall.size() << std::endl;
	for (unsigned int j = 0; j < vPeak.size(); ++j)
		std::cout << "1\t" << vWave.at(iPeak.at(j)) << "\t" << vPeak.at(j) << std::endl;
	for (unsigned int j = 0; j < vVall.size(); ++j)
		std::cout << "0\t" << vWave.at(iVall.at(j)) << "\t" << vVall.at(j) << std::endl;

	for (unsigned int i = 0; i < vInty.size(); ++i)
		Out << vWave.at(i) << "\t" << vInty.at(i) << "\t" << vCopy.at(i) << std::endl;


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
	std::cout <<"\n  -h,  --help\t\t";
	std::cout << "Print this message and exit" << std::endl;
}

std::vector<double> Markov(const std::vector<double> &Spectrum, int Window)
{
	double Area = 0;
	double MaxI = 0;
	for (unsigned int i = 0; i < Spectrum.size(); ++i)
	{
		if (MaxI < Spectrum.at(i))
			MaxI = Spectrum.at(i);
		Area += Spectrum.at(i);
	}

	std::vector<double> vMarkov(Spectrum.size());
	vMarkov.front() = 1.0;
	double Norm = 1.0;
	for (int i = 0; i < Spectrum.size()-1; ++i)
	{
		double sumF = 0, sumB = 0;
		double p0 = Spectrum.at(i) / MaxI;
		double p1 = Spectrum.at(i+1) / MaxI;
		for (int j = 1; j < Window; ++j)
		{
			double pF, ComptF;	//forward
			if (i+j > Spectrum.size()-1)
				pF = Spectrum.at(Spectrum.size()-1);
			else
				 pF = Spectrum.at(i+j)/MaxI;
			if(pF + p0 <= 0)
				ComptF = exp(pF - p0);
			else
				ComptF = exp((pF-p0)/sqrt(pF+p0));
			sumF += ComptF;

			double pB, ComptB;	//backward
			if (i+1-j < 0)
				pB = Spectrum.at(0);
			else
				 pB = Spectrum.at(i+1-j)/MaxI;
			if(pB + p1 <= 0)
				ComptB = exp(pB - p1);
			else
				ComptB = exp((pB - p1)/sqrt(pB+p1));
			sumB += ComptB;
		}

		vMarkov.at(i+1) = vMarkov.at(i) * sumF/sumB;
		Norm += vMarkov.at(i+1);
	}
	for (unsigned int i = 0; i < Spectrum.size(); ++i)
		vMarkov.at(i) *= Area/Norm;

	return vMarkov;
}
