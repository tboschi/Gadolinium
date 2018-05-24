#ifndef TRACKMANAGER_H
#define TRACKMANAGER_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cstdlib>

class TrackManager
{
	public:
		TrackManager(std::string ListFile);

		bool LoadFileList(std::vector<std::string> &vFile, unsigned int n);
		bool LoadTrack(std::vector<std::string> &vFile, std::vector<double> &vTrack, unsigned int nFile);
		bool LoadTrack(std::vector<double> &vTrack, std::string FileName);
		bool AverageTrack(std::vector<double> &vTrack, unsigned int n);
		unsigned int GetEntries();
		double GetPercentage(unsigned int i);
		std::string GetFolderName(unsigned int i);
		unsigned int GetX();
		double GetX(unsigned int i);

	private:
		std::vector<double> vPerc;
		std::vector<std::string> vFold;
		std::vector<double> vXaxis;
};

#endif
