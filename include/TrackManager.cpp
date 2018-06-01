#include "TrackManager.h"

//constructor, requires a configuration file in which each line is a couple (double, string) as (percentage, folder)
//the couples are loeaded in synced vector, vPerc, vFold
TrackManager::TrackManager(std::string ListFile, double V0)
{
	std::string Line, Folder;
	double Percentage, Volume;
	std::stringstream ssL;

	std::ifstream InFile(ListFile.c_str());
	while (std::getline(InFile, Line))
	{
		if (Line[0] == '#') continue;

		ssL.str("");
		ssL.clear();
		ssL << Line;

		ssL >> Percentage >> Volume >> Folder;

		if (Volume < 0)
			BackFold = Folder;
		else
		{
			vPerc.push_back(Percentage*V0/Volume);	//correct concentration
			vFold.push_back(Folder);
		}
	}
	_n = -1;	//unsigned int, so these are INT_MAX
	_i = -1;

	if (BackFold != "")
		AverageBack();
}

//load the files in a folder inside the vector vFile
bool TrackManager::LoadFileList(std::vector<std::string> &vList, unsigned int n)
{
	vList.clear();
	return LoadFileList(vList, vFold.at(n));
}

bool TrackManager::LoadFileList(std::vector<std::string> &vList, std::string Folder)
{
	std::cout << "Reading " << Folder << "\t";
	std::string Cmd = std::string("ls ") + Folder + "/ledj* > .tmp_listfile";
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

		vList.push_back(ssL.str());
	}
	std::cout << " -- " << vList.size() << " files" << std::endl;

	return vList.size();
}

//load the spectrum contained in a single file
//returns the average
double TrackManager::LoadTrack(std::vector<double> &vTrack, std::string FileName)
{
	vTrack.clear();
	bool Zero = (GetX() == 0);

	std::string Line;
	std::stringstream ssL;

	std::ifstream InFile(FileName);
	double yAvg = 0.0;
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
			if (Zero)
				vXaxis.push_back( atof(vWord.at(0).c_str()) );

			double yVal = atof(vWord.at(1).c_str());
			yAvg += yVal;
			vTrack.push_back( yVal );
		}
	}
	
	return yAvg/vTrack.size();
}

void TrackManager::AverageBack()
{
	std::cout << "averaging background " << BackFold << std::endl;
	vBack.clear();
	std::vector<double> vTrack;
	std::vector<std::string> vList;
	if (LoadFileList(vList, BackFold))		//vList is a list of the background files
	{						//looping over the list and averaging the spectra
		for (unsigned int i = 0; i < vList.size(); ++i)		//if n is different from _n starts from 0
		{							//else start from previous status
			double yAvg = LoadTrack(vTrack, vList.at(i));

			if (vBack.size() == 0)
				vBack.resize(2*vTrack.size());	//twice the size if error wanted

			for (unsigned int j = 0; j < vTrack.size(); ++j)
				vBack.at(j) += vTrack.at(j)/vList.size();
		}

		for (unsigned int i = 0; i < vList.size(); ++i)
		{
			double yAvg = LoadTrack(vTrack, vList.at(i));

			for (unsigned int j = 0, e = vTrack.size(); j < vTrack.size(); ++j, ++e)
				vBack.at(e) += pow(vBack.at(j) - vTrack.at(j), 2)/(vList.size()-1);
		}
	}
}

//averages all the spectrum inside folder n
//the variance is appendended at the end of the spectrum, so the track is actually twice as long
bool TrackManager::AverageTrack(std::vector<double> &vAvg, unsigned int n, bool Err)
{
	vAvg.clear();

	std::vector<double> vTrack;
	std::vector<std::string> vList;
	unsigned int Nf = 0;
	double yMean = 0, yRms2 = 0;
	bool Exit = true;
	unsigned int iStart = _n != n ? 0 : _i;
	if (LoadFileList(vList, n))		//vList is a list of the files in the folder n
	{					//looping over the list and averaging the spectra
		unsigned int iEnd = vList.size();
		for (unsigned int i = iStart; i < iEnd; ++i)		//if n is different from _n starts from 0
		{							//else start from previous status
			double yAvg = LoadTrack(vTrack, vList.at(i));
			if (yAvg < 0.04)	//points away from the average, like the led was off
				continue;	//there are discarded
			else
				++Nf;

			yMean = ((Nf-1)*yMean + yAvg)/Nf;
			yRms2 = ((Nf-1)*yRms2 + pow(yMean-yAvg, 2))/Nf;

			if (sqrt(yRms2) > 0.75e-3)	//different set of points basically
			{
				_n = n;		//save current status for later
				//_i = i+1;
				_i = i;
				iEnd = i+1;
				Exit = false;
				break;
			}

			if (vAvg.size() == 0)
				vAvg.resize((Err+1)*vTrack.size());	//twice the size if error wanted

			for (unsigned int j = 0; j < vTrack.size(); ++j)
				vAvg.at(j) += vTrack.at(j);
		}

		if (--Nf > 0)
			for (unsigned int j = 0; j < vAvg.size(); ++j)
				vAvg.at(j) /= Nf;

		for (unsigned int i = iStart; i < iEnd*Err; ++i)
		{
			double yAvg = LoadTrack(vTrack, vList.at(i));
			if (yAvg < 0.04)	//points away from the average, like the led was off
				continue;	//there are discarded
			for (unsigned int j = 0, e = vTrack.size(); j < vTrack.size(); ++j, ++e)
				vAvg.at(e) += pow(vAvg.at(j) - vTrack.at(j), 2)/(Nf-1);
		}
	}

	if (vBack.size() > 0)
	{
		for (unsigned int j = 0, e = vAvg.size()/(1+Err); j < vAvg.size()/(1+Err); ++j, ++e)
		{
			vAvg.at(j) -= vBack.at(j);
			vAvg.at(e) += vBack.at(e);
		}
	}

	std::cout << "... done! Analysing now." << std::endl;

	return Exit;		//full data set ok!
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
