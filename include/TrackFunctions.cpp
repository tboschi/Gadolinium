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

void TrackFunctions::Analyse(const std::vector<double> &vTrack, int iA, int iB, bool Err)
{
	if (iA < -1)
		iA = 0;
	if (iB < -1)
		iB = vTrack.size()/(Err+1);

	fArea = 0.0;
	fAreaVar = -1.0+Err;
	fMax = -100.0;
	fMaxVar = -1.0;
	iMax = 0.0;
	fMin = 100.0;
	fMinVar = -1.0;
	iMin = 0.0;

	for (int i = iA; i < iB; ++i)
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
		if (Err)
			fAreaVar += vTrack.at(i+vTrack.size()/2);
	}
	if (Err)
	{
		fMaxVar = vTrack.at(iMax);
		fMinVar = vTrack.at(iMin);
	}

	fThr = fPerc*(fMax-fMin);
}

double TrackFunctions::Baseline(std::vector<double> &vTrack, bool Subtract, bool Err)
{
	Analyse(vTrack);

	fBaseline = 0, fBaselVar = 0;
	int SampleLength = 0;
	for (unsigned int i = 0, e = vTrack.size()/2; i < vTrack.size()/(1+Err); ++i, ++e)
	{
		double s0 = vTrack.at(i);
		double sum = 0, var = 0;
		unsigned int j = i;
		while (j < vTrack.size() && fabs(vTrack.at(j) - s0) < GetThr())
		{
			if (Err)
			{
				sum += vTrack.at(j)/vTrack.at(e);
				var += 1/vTrack.at(e);
			}
			else
				sum    += vTrack.at(j);
			++j;
		}

		if (j-i > SampleLength)
		{
			if (Err)
			{
				fBaseline = sum/var;
				fBaselVar =   1/var;
			}
			else
				fBaseline = sum/(j-1);

			SampleLength = j-i;
		}
		i = j;		//skip forward
	}

	if (Subtract)
	{
		for (unsigned int i = 0, e = vTrack.size()/2; i < vTrack.size()/(1+Err); ++i, ++e)
		{
			vTrack.at(i) -= fBaseline;
			if (Err)
				vTrack.at(e) += fBaselVar;
		}
	}

	return fBaseline;
}


//normalise to reference track, using highest peak in region [iA:iB]
void TrackFunctions::Normalise(std::vector<double> &vTrack, double Norm, double NormVar)
{
	bool Err = NormVar >= 0.0;

	for (unsigned int i = 0, e = vTrack.size()/2; i < vTrack.size()/(1+Err); ++i, ++e)
	{
		double t = vTrack.at(i);	//o ~ x^2/y^2 (o_y / y + 0_x / x)
		vTrack.at(e) = pow(t/Norm, 2.0) * (vTrack.at(e)/t/t + NormVar/Norm/Norm);
		vTrack.at(i) /= Norm;
	}
}

void TrackFunctions::NormaliseLine(std::vector<double> &vTrack, int Line, bool Err)
{
	double LineVar = -1.0;
	if (Err)
		LineVar = vTrack.at(Line+vTrack.size()/2);

	if (Line > -1 && Line < vTrack.size()/(Err+1))
		return Normalise(vTrack, vTrack.at(Line), LineVar);
}

void TrackFunctions::NormalisePeak(std::vector<double> &vTrack, int iA, int iB, bool Err)
{
	Analyse(vTrack, iA, iB, Err);
	Normalise(vTrack, GetMax(), GetMaxVar());
}

void TrackFunctions::NormaliseArea(std::vector<double> &vTrack, int iA, int iB, bool Err)
{
	Analyse(vTrack, iA, iB, Err);
	Normalise(vTrack, GetArea(), GetAreaVar());
}

std::vector<double> TrackFunctions::Absorption(const std::vector<double> &vRef, 
					       const std::vector<double> &vTrack, 
					       bool Err)
{
	std::vector<double> vAbs(vRef.size());
	for (unsigned int i = 0, e = vRef.size()/2; i < vRef.size()/(1+Err); ++i, ++e)
	{
		if (Err)
			vAbs.at(e) = vRef.at(e)/pow(log(10.0)*vRef.at(i), 2) + 
				     vTrack.at(e)/pow(log(10.0)*vTrack.at(i), 2);
		vAbs.at(i) = log10(vRef.at(i) / vTrack.at(i)); 
	}

	return vAbs;
}

double TrackFunctions::AbsorptionDiff(const std::vector<double> &vTrack, unsigned int A, unsigned int B)
{
	return log10(vTrack.at(B)/vTrack.at(A));
}

double TrackFunctions::AbsorptionDiffErr(const std::vector<double> &vTrack, unsigned int A, unsigned int B)
{
	double Abs = AbsorptionDiff(vTrack, A, B);
	return vTrack.at(A+vTrack.size()/2) / pow(log(10.0)*vTrack.at(A), 2) + 
	       vTrack.at(B+vTrack.size()/2) / pow(log(10.0)*vTrack.at(B), 2);
}

//valleys and peaks of a smoothened track
//gived the position in the vector
/*
void TrackFunctions::FindPeakValley(std::vector<unsigned int> &iPeak, std::vector<unsigned int> &iVall)
{
	FindPeakValley(vTrack, iPeak, iVall)
}
*/

void TrackFunctions::FindPeakValley(const std::vector<double> &vTrack, std::vector<unsigned int> &iPeak, std::vector<unsigned int> &iVall, int iA, int iB, bool Sorting)
{
	if (iA < -1)
		iA = -1;
	if (iB < 0)
		iB = vTrack.size();

	Analyse(vTrack, iA, iB);	//needed
	std::deque<unsigned int> dPeak, dVall;

	dPeak.push_back(GetMax_i());

	//going left first, and then right
	for (int Dir = -1; Dir < 2; Dir += 2)
	{
		double PoV = false;		//F looking for Valley, T looking for P
		int iD = GetMax_i(), iS = GetMax_i();
		double fS = GetMax();
		while (iD > iA && iD < iB)
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
					std::deque<unsigned int> &dRef = PoV ? dPeak : dVall;

					if (Dir < 0)
						dRef.push_front(iZ);
					else
						dRef.push_back(iZ);

					fS = fZ;
					iD = iS = iZ;
					PoV = !PoV;
				}
			}
			iD += Dir;
		}
	}
	if (Sorting)
	{
		std::sort(dPeak.begin(), dPeak.end(), Sort(vTrack,  1));
		std::sort(dVall.begin(), dVall.end(), Sort(vTrack, -1));
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
	Analyse(vTrack);
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

double TrackFunctions::GetAreaVar()
{
	return fAreaVar;
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

double TrackFunctions::GetMaxVar()
{
	return fMaxVar;
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
