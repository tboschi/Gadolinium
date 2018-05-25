#include "TrackManager.h"

//constructor, requires a configuration file in which each line is a couple (double, string) as (percentage, folder)
//the couples are loeaded in synced vector, vPerc, vFold
TrackManager::TrackManager(std::string ListFile)
{
	std::string Line, Folder;
	double Percentage;
	std::stringstream ssL;

	std::ifstream InFile(ListFile.c_str());
	while (std::getline(InFile, Line))
	{
		if (Line[0] == '#') continue;

		ssL.str("");
		ssL.clear();
		ssL << Line;

		ssL >> Percentage >> Folder;

		vPerc.push_back(Percentage);
		vFold.push_back(Folder);
	}
}

//load the files in a folder inside the vector vFile
bool TrackManager::LoadFileList(std::vector<std::string> &vFile, unsigned int n)
{
	vFile.clear();

	std::string Cmd = std::string("ls ") + vFold.at(n) + "/ledj* > .tmp_listfile";
	system(Cmd.c_str());

	std::string Line;
	std::stringstream ssL;

	std::ifstream InFile(".tmp_listfile");
	while (std::getline(InFile, Line))
	{
		if (Line[0] == '#') continue;

		ssL.str("");
		ssL.clear();
		ssL << Line;
		//ssL << vFold.at(n) << "/" << Line;

		vFile.push_back(ssL.str());
	}

	return vFile.size();
}

//load the spectrum contained in a single file
bool TrackManager::LoadTrack(std::vector<double> &vTrack, std::string FileName)
{
	vTrack.clear();
	bool Test = (GetX() == 0);

	std::string Line;
	std::stringstream ssL;

	std::ifstream InFile(FileName);
	while (std::getline(InFile, Line))
	{
		if (Line[0] == '#') continue;

		ssL.str("");
		ssL.clear();
		ssL << Line;

		std::string Word;
		std::vector<std::string> vWord;
		while (std::getline(ssL, Word, ';'))
			vWord.push_back(Word);

		if (vWord.size() == 2)
		{
			if (Test)
				vXaxis.push_back( atof(vWord.at(0).c_str()) );

			vTrack.push_back( atof(vWord.at(1).c_str()) );
		}
	}

	return vTrack.size();
}

//averages all the spectrum inside folder n
//the variance is appendended at the end of the spectrum, so the track is actually twice as long
bool TrackManager::AverageTrack(std::vector<double> &vAvg, unsigned int n, bool Err)
{
	vAvg.clear();

	std::vector<double> vTrack;
	std::vector<std::string> vList;
	if (LoadFileList(vList, n))		//vList is a list of the files in the folder n
	{					//looping over the list and averaging the spectra
		for (unsigned int i = 0; i < vList.size(); ++i)
		{
			LoadTrack(vTrack, vList.at(i));

			if (vAvg.size() == 0)
				vAvg.resize((Err+1)*vTrack.size());	//twice the size if error wanted

			for (unsigned int i = 0; i < vTrack.size(); ++i)
				vAvg.at(i) += vTrack.at(i)/vList.size();
		}
		if (Err)
			for (unsigned int i = 0, j = vTrack.size(); i < vTrack.size(); ++i, ++j)
				vAvg.at(j) += pow(vAvg.at(i)-vTrack.at(i), 2)/(vList.size()-1);
	}

	return vAvg.size();
}

unsigned int TrackManager::GetEntries()
{
	return vPerc.size();
}

double TrackManager::GetPercentage(unsigned int i)
{
	if (i < GetEntries())
		return vPerc.at(i);
	else return -1.0;
}

std::string TrackManager::GetFolderName(unsigned int i)
{
	if (i < GetEntries())
		return vFold.at(i);
	else return std::string("");
}

unsigned int TrackManager::GetX()
{
	return vXaxis.size();
}

double TrackManager::GetX(unsigned int i)
{
	if (i < GetX())
		return vXaxis.at(i);
	else return -1.0;
}
