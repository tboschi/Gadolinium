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

	std::string Cmd = std::string("ls ") + vFold.at(n) + " > .tmp_listfile";
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

		vFile.push_back(Line);
	}

	return vFile.size();
}

//load the spectrum contained in a single file
//requires the list of files
bool TrackManager::LoadTrack(std::vector<std::string> &vFile, std::vector<double> &vTrack, unsigned int nFile)
{
	LoadTrack(vTrack, vFile.at(nFile));
}

bool TrackManager::LoadTrack(std::vector<double> &vTrack, std::string FileName)
{
	std::string Line;
	std::stringstream ssL;

	std::ifstream InFile(FileName);
	std::cout << "
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
			if (vXaxis.size() == 0)
				vXaxis.push_back( atof(vWord.at(0).c_str()) );

			vTrack.push_back( atof(vWord.at(1).c_str()) );
		}
	}

	return vTrack.size();
}

//averages all the spectrum inside a folder
bool TrackManager::AverageTrack(std::vector<double> &vTrack, unsigned int n)
{
	vTrack.clear();

	std::vector<std::string> vList;
	std::cout << "TM 0" << std::endl;
	if (LoadFileList(vList, n))		//vList is a list of the files in the folder n
	{					//looping over the list and averaging the spectra
		std::cout << "if TM 0 " << vList.size() << std::endl;
		for (unsigned int i = 0; i < vList.size(); ++i)
		{
			vTrack.clear();
			std::cout << "TM L " << vList.at(i) << std::endl;
			LoadTrack(vTrack, vList.at(i));

			if (vTrack.size() == 0)
				vTrack.resize(vTrack.size());

			for (unsigned int i = 0; i < vTrack.size(); ++i)
				vTrack.at(i) += vTrack.at(i)/vList.size();
		}
		std::cout << "if TM 1" << std::endl;
	}

	std::cout << "TM 1" << std::endl;
	return vTrack.size();
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
