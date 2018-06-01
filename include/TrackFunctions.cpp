#include "TrackFunctions.h"

TrackFunctions::TrackFunctions(double Percentage) :
	fPerc(Percentage)
{
}

void TrackFunctions::SetPercentage(double Percentage)
{
	fPerc = Percentage;
}

void TrackFunctions::SetInterval(TrackManager* TM, double A, double B)
{
	iA = 0;
	iB = 0;

	if (A > B)
	{
		double tmp = A;
		A = B;
		B = tmp;
	}

	for (unsigned int i = 0; i < TM->GetX(); ++i)
	{
		if (TM->GetX(i) < A && A < TM->GetX(i+1))
			iA = i;
		if (TM->GetX(i) < B && B < TM->GetX(i+1))
			iB = i;
	}
	std::cout << "from iA " << iA << "\t" << " to iB " << iB << std::endl;
}

void TrackFunctions::Analyse(const std::vector<double> &vTrack, bool Err)
{
	if (iA < 0)
		iA = 0;
	if (iB > vTrack.size()/(1+Err))
		iB = vTrack.size()/(Err+1);

	fArea = 0.0;
	fA1  = 0.0;
	fA2  = 0.0;
	fAreaVar = -(!Err);
	fA1Var   = -(!Err);	//-1 if no error
	fA2Var   = -(!Err);	//0 if error
	fT1 = (2*iA + iB)/3;
	fT2 = (iA + 2*iB)/3;

	fMax    = -100.0;
	fMin    =  100.0;
	fMaxVar = -1.0;
	fMinVar = -1.0;
	iMax    =  0.0;
	iMin    =  0.0;

	unsigned int Nf1 = 0, Nf2 = 2;
	for (unsigned int j = iA, e = iA + vTrack.size()/(1+Err); j < iB; ++j, ++e)
	{
		if (fMax < vTrack.at(j))
		{
			fMax = vTrack.at(j);
			iMax = j;
		}
		if (fMin > vTrack.at(j))
		{
			fMin = vTrack.at(j);
			iMin = j;
		}


		fArea += vTrack.at(j);
		if (Err)
			fAreaVar += vTrack.at(e);
		if (j < fT1)
		{
			fA1    += vTrack.at(j);
			if (Err)
				fA1Var += vTrack.at(e);
			++Nf1;
		}
		if (j > fT2)
		{
			fA2    += vTrack.at(j);
			if (Err)
				fA2Var += vTrack.at(e);
			++Nf2;
		}
	}

	//min and max error
	if (Err)
	{
		fMaxVar = vTrack.at(iMax+vTrack.size()/2);
		fMinVar = vTrack.at(iMin+vTrack.size()/2);
	}

	fThr = fPerc*(fMax-fMin);
}

double TrackFunctions::Baseline(std::vector<double> &vTrack, bool Subtract, bool Err)
{
	Analyse(vTrack, Err);

	fBaseline = 0, fBaselVar = 0;
	int SampleLength = 0;
	std::cout << "the " << GetThr() << std::endl;
	unsigned int A, B;
	for (unsigned int i = 0, e = vTrack.size()/2; i < vTrack.size()/(1+Err); ++i, ++e)
	{
		double s0 = vTrack.at(i);
		double sum = 0, var = 0;
		unsigned int j = i;
		while (j < vTrack.size()/(1+Err) && fabs(vTrack.at(j) - s0) < GetThr())
		{
			if (Err)
			{
				sum += vTrack.at(j)/vTrack.at(e);
				var += 1/vTrack.at(e);
			}
			else
				sum += vTrack.at(j);
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
			A = i;
			B = j;
		}
		i = j;		//skip forward
	}

	std::cout << "baseline from " << A << "\tto\t" << B << std::endl;
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

	for (unsigned int j = iA, e = iA + vTrack.size()/(1+Err); j < iB; ++j, ++e)
	{
		double t = vTrack.at(j);	//o ~ x^2/y^2 (o_y / y + 0_x / x)
		vTrack.at(e) = pow(t/Norm, 2.0) * (vTrack.at(e)/t/t + NormVar/Norm/Norm);
		vTrack.at(j) /= Norm;
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

void TrackFunctions::NormalisePeak(std::vector<double> &vTrack, bool Err)
{
	Analyse(vTrack, Err);
	Normalise(vTrack, GetMax(), GetMaxVar());
}

//using T1 and Q3 to define two wings and normalise with median value
void TrackFunctions::NormaliseArea(std::vector<double> &vTrack, bool Err)
{
	Analyse(vTrack, Err);
	double VarA = Err ? 1 / (1.0/GetAreaT1Var() + 1.0/GetAreaT2Var()) : -1.0;
	double AvgA = (GetAreaT1()/GetAreaT1Var() + GetAreaT2()/GetAreaT2Var()) * VarA;

	Normalise(vTrack, AvgA, VarA);
}

std::vector<double> TrackFunctions::Absorption(const std::vector<double> &vRef, 
					       const std::vector<double> &vTrack, 
					       bool Err)
{
	std::vector<double> vAbs(vRef.size());
	for (unsigned int j = iA, e = iA + vRef.size()/(1+Err); j < iB; ++j, ++e)
	{
		if (Err)
			vAbs.at(e) = vRef.at(e)/pow(vRef.at(j), 2) + vTrack.at(e)/pow(vTrack.at(j), 2);
		vAbs.at(j) = log(vRef.at(j) / vTrack.at(j)); 
	}

	return vAbs;
}

double TrackFunctions::AbsorptionDiff(const std::vector<double> &vTrack, unsigned int A, unsigned int B)
{
	return log(vTrack.at(A)/vTrack.at(B));
}

double TrackFunctions::AbsorptionDiffErr(const std::vector<double> &vTrack, unsigned int A, unsigned int B)
{
	return vTrack.at(A+vTrack.size()/2)/pow(vTrack.at(A), 2) + 
	       vTrack.at(B+vTrack.size()/2)/pow(vTrack.at(B), 2);
}

//valleys and peaks of a smoothened track
//gived the position in the vector
/*
void TrackFunctions::FindPeakValley(std::vector<unsigned int> &iPeak, std::vector<unsigned int> &iVall)
{
	FindPeakValley(vTrack, iPeak, iVall)
}
*/

void TrackFunctions::FindPeakValley(const std::vector<double> &vTrack, std::vector<unsigned int> &iPeak, std::vector<unsigned int> &iVall, bool Sorting)
{
	Analyse(vTrack);	//needed
	std::deque<unsigned int> dPeak, dVall;

	dPeak.push_back(GetMax_i());

	//going left first, and then right
	for (int Dir = -1; Dir < 2; Dir += 2)
	{
		double PoV = false;		//F looking for Valley, T looking for P
		int iD = GetMax_i(), iS = GetMax_i();
		double fS = GetMax();
		while (iD > GetT1() && iD < GetT2())
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

double TrackFunctions::GetT1()
{
	return fT1;
}

double TrackFunctions::GetT2()
{
	return fT2;
}

double TrackFunctions::GetAreaT1()
{
	return fA1;
}

double TrackFunctions::GetAreaT2()
{
	return fA2;
}

double TrackFunctions::GetAreaT1Var()
{
	return fA1Var;
}

double TrackFunctions::GetAreaT2Var()
{
	return fA2Var;
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
