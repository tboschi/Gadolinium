#ifndef TRACKMANAGER_H
#define TRACKMANAGER_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cstdlib>
#include <cmath>

class TrackManager
{
	public:
		TrackManager(std::string ListFile, double V);

		bool LoadFileList(std::vector<std::string> &vFile, unsigned int n);
		bool LoadFileList(std::vector<std::string> &vFile, std::string Folder);
		double LoadTrack(std::vector<double> &vTrack, std::string FileName);
		void AverageBack();
		bool AverageTrack(std::vector<double> &vAvg, unsigned int n, bool Err);
		unsigned int GetEntries();
		double GetPercentage(unsigned int i);
		std::string GetFolderName(unsigned int i);
		unsigned int GetX();
		double GetX(unsigned int i);

	private:
		std::vector<double> vPerc;
		std::vector<std::string> vFold;

		//backgrund
		std::vector<double> vBack;
		std::string BackFold;

		std::vector<double> vXaxis;

		unsigned int _n, _i;
};

#endif
