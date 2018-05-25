#ifndef TRACKANALYSIS_H
#define TRACKANALYSIS_H

#include <iostream>
#include <vector>
#include <deque>
#include <cstring>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <algorithm>

struct Sort
{
	std::vector<double> vX;
	bool bPoV;
	bool operator()(int i, int j) const { return (2*bPoV - 1)*(vX.at(i)-vX.at(j)) > 0; }

	Sort(const std::vector<double> &vv, bool PoV = true) : vX(vv), bPoV(PoV) {}
};

class TrackFunctions
{
	public:
		TrackFunctions(double Percentage = 0.02);

		void SetPercentage(double Percentage);
		void Analyse(const std::vector<double> &vTrack, int iA = -100, int iB = -100, bool Err = false);
		double Baseline(std::vector<double> &vTrack, bool Subtract = false, bool Err = false);
		void Normalise(std::vector<double> &vTrack, double Norm = 1, double NormVar = -1.0);
		void NormaliseLine(std::vector<double> &vTrack, int Line = 0, bool Err = false);
		void NormalisePeak(std::vector<double> &vTrack, int iA = -100, int iB = -100, bool Err = false);
		void NormaliseArea(std::vector<double> &vTrack, int iA = -100, int iB = -100, bool Err = false);
		std::vector<double> Absorption(const std::vector<double> &vRef, 
					       const std::vector<double> &vTrack, bool Err = false);
		double AbsorptionDiff(const std::vector<double> &vTrack, unsigned int A, unsigned int B);
		double AbsorptionDiffErr(const std::vector<double> &vTrack, unsigned int A, unsigned int B);
		void FindPeakValley(const std::vector<double> &vTrack, 
					  std::vector<unsigned int> &iPeak, 
					  std::vector<unsigned int> &iVall, 
				    int iA = -100, int iB = -100, bool Sorting = false);
		void Smoothen(std::vector<double> &vTrack, unsigned int Integrate);
		void Markov(std::vector<double> &vTrack, int Window);

		double GetPercentage();
		double GetArea();
		double GetAreaVar();
		double GetThr();
		double GetBaseline();
		double GetMax();
		double GetMaxVar();
		double GetMin();
		double GetMinVar();
		unsigned int GetMax_i();
		unsigned int GetMin_i();

	private:
		double fPerc;
		double fArea;
		double fAreaVar;
		double fThr;
		double fBaseline;
		double fBaselVar;
		double fMax;
		double fMaxVar;
		double fMin;
		double fMinVar;
		int iMax;
		int iMin;

};

#endif

