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
		void Analyse(const std::vector<double> &vTrack, int iA = -100, int iB = -100);
		double Baseline(std::vector<double> &vTrack, bool Subtract = false);
		void Normalise(std::vector<double> &vTrack, double Norm = 1);
		void NormaliseLine(std::vector<double> &vTrack, int Line = 0);
		void NormalisePeak(std::vector<double> &vTrack, int iA = -100, int iB = -100);
		void NormaliseArea(std::vector<double> &vTrack, int iA = -100, int iB = -100);
		bool Absorption(const std::vector<double> &vRef, std::vector<double> &vTrack);
		bool AbsorptionVar(const std::vector<double> &vRef, const std::vector<double> &vRefVar, 
				   std::vector<double> &vTrack, std::vector<double> &vTrackRef);
		double AbsorptionDiff(const std::vector<double> &vTrack, unsigned int A, unsigned int B);
		double AbsorptionDiffVar(const std::vector<double> &vTrack, const std::vector<double> &vTrackVar, 
					 unsigned int A, unsigned int B);
		void FindPeakValley(const std::vector<double> &vTrack, std::vector<unsigned int> &iPeak, std::vector<unsigned int> &iVall, 
				    int iA = -100, int iB = -100, bool Sorting = false);
		void Smoothen(std::vector<double> &vTrack, unsigned int Integrate);
		void Markov(std::vector<double> &vTrack, int Window);

		double GetPercentage();
		double GetArea();
		double GetThr();
		double GetBaseline();
		double GetMax();
		unsigned int GetMax_i();
		double GetMin();
		unsigned int GetMin_i();

	private:
		double fPerc;
		double fArea;
		double fThr;
		double fBaseline;
		double fMax;
		int iMax;
		double fMin;
		int iMin;

};

#endif

