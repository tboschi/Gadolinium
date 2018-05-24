#include "TrackFunctions.h"

TrackFunctions::TrackFunctions(double Percentage) :
	fPerc(Percentage)
{
}

void TrackFunctions::SetPercentage(double Percentage)
{
	fPerc = Percentage;
}

/*
void TrackFunctions::SetTrack(std::vector<double> &vSpectrum)
{
	vTrack.clear();
	vTrack.insert(vTrack.end(), vSpectrum.begin(), vSpectrum.end());
}

void TrackFunctions::SetReferenceTrack(std::vector<double> &vSpectrum)
{
	vTRef.clear();
	vTRef.insert(vTRef.end(), vTRef.begin(), vTRef.end());
}
*/

void TrackFunctions::Analyse(const std::vector<double> &vTrack)
{
	fArea = 0;
	fMax = 0;
	iMax = 0;
	fMin = 0;
	iMin = 0;

	for (unsigned int i = 0; i < vTrack.size(); ++i)
	{
		if (fMax < vTrack.at(i))
		{
			fMax = vTrack.at(i);
			iMax = i;
		}
		if (fMin > vTrack.at(i))
		{
			fMin = vTrack.at(i);
			iMin = i;
		}

		fArea += vTrack.at(i);
	}

	fThr = fPerc*(fMax-fMin);
}

double TrackFunctions::Baseline(std::vector<double> &vTrack, bool Subtract)
{
	fBaseline = 0;
	int SampleLength = 0;
	for (unsigned int i = 0; i < vTrack.size(); ++i)
	{
		double s0 = vTrack.at(i);
		double sum = 0;
		unsigned int j = i;
		while (fabs(vTrack.at(j) - s0) < GetThr())
		{
			sum += vTrack.at(j);
			++j;
		}

		if (j-i > SampleLength)
		{
			fBaseline = sum/(j-i);
			SampleLength = j-i;
		}
		i = j;		//skip forward
	}

	if (Subtract)
	{
		for (unsigned int i = 0; i < vTrack.size(); ++i)
			vTrack.at(i) -= fBaseline;
	}

	return fBaseline;
}


//normalise to reference track, using highest peak in region [iA:iB]
bool TrackFunctions::Normalise(std::vector<double> &vTrack, unsigned int iA, unsigned int iB)
{
	if (iA > iB)
	{
		unsigned int tmp = iB;
		iB = iA;
		iA = tmp;
	}

	if (iA > vTrack.size())
	       return false;
	if (iB > vTrack.size())
		iB = vTrack.size();

	double fNorm = 0;
	for (unsigned int j = iA; j < iB; ++j)
		if (fNorm < vTrack.at(j))
			fNorm = vTrack.at(j);

	for (unsigned int j = 0; j < vTrack.size(); ++j)
		vTrack.at(j) / fNorm;
}


bool TrackFunctions::Absorption(const std::vector<double> &vTrack0, std::vector<double> &vTrack)
{
	for (unsigned int i = 0; i < vTrack0.size(); ++i)
		vTrack.at(i) = log10(vTrack0.at(i) / vTrack.at(i)); 
}


//valleys and peaks of a smoothened track
//gived the position in the vector
/*
void TrackFunctions::FindPeakValley(std::vector<unsigned int> &iPeak, std::vector<unsigned int> &iVall)
{
	FindPeakValley(vTrack, iPeak, iVall)
}
*/

void TrackFunctions::FindPeakValley(const std::vector<double> &vTrack, std::vector<unsigned int> &iPeak, std::vector<unsigned int> &iVall)
{
	Analyse(vTrack);	//needed
	//std::deque<double> yPeak, yVall;
	//std::deque<double> xPeak, xVall;
	std::deque<unsigned int> dPeak, dVall;

	iPeak.push_back(GetMax_i());

	//going left first, and then right
	for (int Dir = -1; Dir < 2; Dir += 2)
	{
		double PoV = false;		//F looking for Valley, T looking for P
		int iD = GetMax_i(), iS = GetMax_i();
		double fS = GetMax();
		while (iD > -1 && iD < vTrack.size())
		{
			int Sign = 2*PoV - 1;	//-1 looking for Valley, +1 looking for Peak
			if ( Sign*(vTrack.at(iD) - fS) > GetThr())
			{
				iS = iD ;
				fS = vTrack.at(iD);
			}
			if (-Sign*(vTrack.at(iD) - fS) > GetThr())
			{
				double fZ = -Sign*GetMax();
				int    iZ = -1;
				for (int j = iD; Dir*(iS-j) < 1 && j > -1 && j < vTrack.size(); j -= Dir)
				{
					if (Sign*(vTrack.at(j)-fZ) > 0)
					{
						fZ = vTrack.at(j);
						iZ = j;
					}
				}
				if (iZ > -1)
				{
					//std::deque<double> &yRef = PoV ? yPeak : yVall;
					//std::deque<double> &xRef = PoV ? xPeak : xVall;
					std::deque<unsigned int> &dRef = PoV ? dPeak : dVall;

					if (Dir < 0)
					{
						//yRef.push_front(fZ);
						//xRef.push_front(xWave.at(iZ));
						dRef.push_front(iZ);
					}
					else
					{
						//yRef.push_back(fZ);
						//xRef.push_back(xWave.at(iZ));
						dRef.push_back(iZ);
					}
					fS = fZ;
					iD = iS = iZ;
					PoV = !PoV;
				}
			}
			iD += Dir;
		}
	}

	iPeak.clear();
	iPeak.insert(iPeak.end(), dPeak.begin(), dPeak.end());
	iVall.clear();
	iVall.insert(iVall.end(), dVall.begin(), dVall.end());
}

/*
void TrackFunctions::Smoothen(std::vector<double> &vTrack, unsigned int Integrate)
{
	SetTrack(vTrack);
	vTrack.swap(Smoothen(Integrate));
}
*/

void TrackFunctions::Smoothen(std::vector<double> &vTrack, unsigned int Integrate)
{
	for (int rep = 1; rep < Integrate; ++rep)
		for (int i = rep; i < vTrack.size()-rep; ++i)
			vTrack.at(i) = (vTrack.at(i+rep) + vTrack.at(i-rep))/2.0;
}

/*
void TrackFunctions::Markov(std::vector<double> &vTrack, int Window)
{
	SetTrack(vTrack);
	vTrack.swap(Markov(Window));
}
*/

void TrackFunctions::Markov(std::vector<double> &vTrack, int Window)
{
	std::vector<double> vMarkov(vTrack.size());
	vMarkov.front() = 1.0;
	double Norm = 1.0;
	for (int i = 0; i < vTrack.size()-1; ++i)
	{
		double sumF = 0, sumB = 0;
		double p0 = vTrack.at(i) / GetMax();
		double p1 = vTrack.at(i+1) / GetMax();
		for (int j = 1; j < Window; ++j)
		{
			double pF, ComptF;	//forward
			if (i+j > vTrack.size()-1)
				pF = vTrack.at(vTrack.size()-1);
			else
				 pF = vTrack.at(i+j)/GetMax();
			if(pF + p0 <= 0)
				ComptF = exp(pF - p0);
			else
				ComptF = exp((pF-p0)/sqrt(pF+p0));
			sumF += ComptF;

			double pB, ComptB;	//backward
			if (i+1-j < 0)
				pB = vTrack.at(0);
			else
				 pB = vTrack.at(i+1-j)/GetMax();
			if(pB + p1 <= 0)
				ComptB = exp(pB - p1);
			else
				ComptB = exp((pB - p1)/sqrt(pB+p1));
			sumB += ComptB;
		}

		vMarkov.at(i+1) = vMarkov.at(i) * sumF/sumB;
		Norm += vMarkov.at(i+1);
	}
	for (unsigned int i = 0; i < vTrack.size(); ++i)
		vMarkov.at(i) *= GetArea()/Norm;

	vTrack.swap(vMarkov);
}

double TrackFunctions::GetPercentage()
{
	return fPerc;
}

double TrackFunctions::GetArea()
{
	return fArea;
}

double TrackFunctions::GetThr()
{
	return fThr;
}

double TrackFunctions::GetBaseline()
{
	return fBaseline;
}

double TrackFunctions::GetMax()
{
	return fMax;
}

unsigned int TrackFunctions::GetMax_i()
{
	return iMax;
}

double TrackFunctions::GetMin()
{
	return fMin;
}

unsigned int TrackFunctions::GetMin_i()
{
	return iMin;
}
