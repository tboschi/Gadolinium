#ifndef TRACKANALYSIS_H
#define TRACKANALYSIS_H

#include <iostream>
#include <vector>
#include <deque>
#include <cstring>
#include <sstream>
#include <cstdlib>
#include <cmath>

class TrackFunctions
{
	public:
		TrackFunctions(double Percentage = 0.02);

		void SetPercentage(double Percentage);
		void Analyse(const std::vector<double> &vTrack);
		double Baseline(std::vector<double> &vTrack, bool Subtract = false);
		bool Normalise(std::vector<double> &vTrack, unsigned int iA, unsigned int iB);
		bool Absorption(const std::vector<double> &vTrack0, std::vector<double> &vTrack);
		void FindPeakValley(const std::vector<double> &vTrack, std::vector<unsigned int> &iPeak, std::vector<unsigned int> &iVall);
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

