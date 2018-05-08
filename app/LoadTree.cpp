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

class Sort
{
	private:
		std::vector<double> vX;
	public:
		Sort(std::vector<double> &vv) : vX(vv) {}
		bool operator()(int i, int j) const { return vX.at(i) < vX.at(j); }
};

void Usage(char* Name);
std::vector<double> Markov(const std::vector<double> &Spectrum, int Window);
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
	
	while((iarg = getopt_long(argc,argv, "I:i:r:o:l:h", longopts, &index)) != -1)
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
				SpecFile.open(optarg);
				break;
			case 'l':
				PeakFile.open(optarg);
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

	std::ostream &OutS = (SpecFile.is_open()) ? SpecFile : std::cout;
	std::ostream &OutP = (PeakFile.is_open()) ? PeakFile : std::cout;

	double Wavelength, Intensity;
	std::vector<double> xWave, yInty;

	std::string Line;
	std::stringstream ssL;
	while (std::getline(InFile, Line))
	{
		if (Line[0] == '#') continue;

		ssL.clear();
		ssL.str("");

		ssL << Line;
		ssL >> Wavelength >> Intensity;

		xWave.push_back(Wavelength);
		yInty.push_back(Intensity);
	}
	//					4 points per 100nm = 400 bins
	TH1D* hIntensity = new TH1D("hintensity", "Intensity", 4000, 100, 1100);
	TH1D* hAbsorbanc = new TH1D("habsorbanc", "Absorbanc", 4000, 100, 1100);

	TTree *tPeak = new TTree("Peak", "Peak");
	TTree *tVall = new TTree("Vall", "Valley");

	double tx, ty;
	tPeak->Branch("x", &tx, "tx/D");
	tPeak->Branch("y", &ty, "ty/D");
	tVall->Branch("x", &tx, "tx/D");
	tVall->Branch("y", &ty, "ty/D");
	
	//smoothen spectrum
	std::vector<double> yWave(yInty);
	for (int rep = 1; rep < Integ; ++rep)
		for (int i = rep; i < yWave.size()-rep; ++i)
			yWave.at(i) = (yWave.at(i+rep) + yWave.at(i-rep))/2.0;

					
	//yWave is almost backgroundless now
	//find the lowest and highest point of it
	//to subtract baseline autmatically
	//
	double Min = 999, Max = -999;
	unsigned int iMin, iMax;
	for (unsigned int i = 0; i < yWave.size(); ++i)
	{
		if (Min > yWave.at(i))
		{
			Min = yWave.at(i);
			iMin = i;
		}
		if (Max < yWave.at(i))
		{
			Max = yWave.at(i);
			iMax = i;
		}
	}

	double Perc = 0.02;	//2% of the peak is start
	double Thr = Perc*(Max-Min);

	/*
	double Baseline = 0, Rms = 0;
	int SampleLength = 0;
	int iA, iB;
	for (unsigned int i = 0; i < xWave.size(); ++i)
	{
		double s0 = yInty.at(i);
		double sum = 0;
		unsigned int j = i;
		while (fabs(yWave.at(j) - s0) < Thr)
		{
			sum += yWave.at(j);
			++j;
		}

		if (j-i > SampleLength)
		{
			Baseline = sum/(j-i);
			SampleLength = j-i;
			iA = i, iB = j;
		}
		i = j;		//skip forward
	}
	*/

	//for (unsigned int i = 0; i < xWave.size(); ++i)
	//	yWave.at(i) -= Baseline;

	//find all peak and valleys
	//
	std::deque<double> yPeak, yVall;
	std::deque<double> xPeak, xVall;

	double xMax = xWave.at(iMax);
	double yMax = yWave.at(iMax);
	xPeak.push_back(xMax);
	yPeak.push_back(yMax);

	//going left first, and then right
	for (int Dir = -1; Dir < 2; Dir += 2)
	{
		double PoV = false;		//F looking for Valley, T looking for P
		int iD = iMax, iS = iMax;
		double fS = yMax;
		while (iD > -1 && iD < xWave.size())
		{
			int Sign = 2*PoV - 1;	//-1 looking for Valley, +1 looking for Peak
			if ( Sign*(yWave.at(iD) - fS) > Thr)
			{
				iS = iD ;
				fS = yWave.at(iD);
			}
			if (-Sign*(yWave.at(iD) - fS) > Thr)
			{
				double fZ = -Sign*yPeak.front();
				int    iZ = -1;
				for (int j = iD; Dir*(iS-j) < 1 && j > -1 && j < xWave.size(); j -= Dir)
				{
					if (Sign*(yWave.at(j)-fZ) > 0)
					{
						fZ = yWave.at(j);
						iZ = j;
					}
				}
				if (iZ > -1)
				{
					std::deque<double> &yRef = PoV ? yPeak : yVall;
					std::deque<double> &xRef = PoV ? xPeak : xVall;

					if (Dir < 0)
					{
						yRef.push_front(fZ);
						xRef.push_front(xWave.at(iZ));
					}
					else
					{
						yRef.push_back(fZ);
						xRef.push_back(xWave.at(iZ));
					}
					fS = fZ;
					iD = iS = iZ;
					PoV = !PoV;
				}
			}
			iD += Dir;
		}
	}

	OutP << "#Peaks\t" << xPeak.size() << "\tValleys " << xVall.size() << std::endl;
	for (unsigned int j = 0; j < xPeak.size(); ++j)
	{
		tx = xPeak.at(j);
		ty = yPeak.at(j);
		tPeak->Fill();
		OutP << "1\t" << xPeak.at(j) << "\t" << yPeak.at(j)/Max << std::endl;
	}
	for (unsigned int j = 0; j < xVall.size(); ++j)
	{
		tx = xVall.at(j);
		ty = yVall.at(j);
		tVall->Fill();
		OutP << "0\t" << xVall.at(j) << "\t" << yVall.at(j)/Max << std::endl;
	}

	//Creating absorbance spectrum
	std::vector<double> vAmp, vAbs, vLine;
	for (unsigned int i = 0, j = 1; i < xWave.size(); ++i)
	{
		double y0;
		if (xWave.at(i) < xPeak.back())
		{
			double m = (yPeak.at(j) - yPeak.at(j-1))/(xPeak.at(j) - xPeak.at(j-1));
			y0 = m * (xWave.at(i) - xPeak.at(j)) + yPeak.at(j);
			if (j < xPeak.size() - 1 && xWave.at(i) > xPeak.at(j))
				++j;
		}
		else 
			y0 = yWave.at(i);
		vAmp.push_back( y0-yWave.at(i) );
		vAbs.push_back( log10(y0/yWave.at(i)) );
	}

	for (unsigned int i = 0; i < xWave.size(); ++i)
	{
		hIntensity->Fill(xWave.at(i), yWave.at(i));
		hAbsorbanc->Fill(xWave.at(i), vAbs.at(i));
		OutS << xWave.at(i) << "\t" << yInty.at(i)/Max << "\t" << yWave.at(i)/Max << "\t" << vAmp.at(i)/Max << "\t" << vAbs.at(i) << std::endl;
	}

	if (RootFile->IsOpen())
	{
		hIntensity->Write();
		hAbsorbanc->Write();
		tPeak->Write();
		tVall->Write();
		RootFile->Close();
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

//std::vector<double> vMarkov = Markov(vInty, Integ/2);
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
