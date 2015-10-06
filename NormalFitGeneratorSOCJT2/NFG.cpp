#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <cstdio>
//#include <tgmath.h>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <ctime>
#include <locale>
#include <windows.h>
#include "omp.h"
#include <sys/stat.h>

using namespace std;

/* This program takes a SOCJT2 fit file and outputs a fit file using the same assignments
where the levels are varied using a normal distribution with standard deviation hard
coded below. This was made to be used to create varying "experiment" data sets from
a pure calculation and run fitting tests against such.

USAGE:
- Make sure to have SOCJT 2.exe in the same folder.
- The fit file is formatted exactly as for SOCJT 2, and each iteration is saved as
(fit output name).(iteration index).nfgfit
- The output files are regular output files and are named
(output name).(iteration index).nfgout
- In the input file, set the fit file to "USENFG"
- The total output file is a file with all the fit values listed in columns left for data
  analysis. */

/* This is the Box Muller Method for generating random numbers following a normal distribution
where the two arguments are the mean or expectation, and the standard deviation 
In particular, this is the POLAR FORM of the Box Muller transformation which is computationally
more efficient (not that it really matters...) and a bit safer when x and y are closer to 0. */
double rand_normal(double mean, double stddev)
{
	static double n2 = 0.0;
	static int n2_cached = 0;
	if (!n2_cached)
	{
		double x, y, r;
		do
		{
			x = 2.0*rand() / RAND_MAX - 1;
			y = 2.0*rand() / RAND_MAX - 1;

			r = x*x + y*y;
		} while (r == 0.0 || r > 1.0);
		{
			double d = sqrt(-2.0*log(r) / r);
			double n1 = x*d;
			n2 = y*d;
			double result = n1*stddev + mean;
			n2_cached = 1;
			return result;
		}
	}
	else
	{
		n2_cached = 0;
		return n2*stddev + mean;
	}
}

/* Gets the CWD */
string GetCWD()
{
	char buffer[MAX_PATH];
	GetModuleFileName(NULL, buffer, MAX_PATH);
	string::size_type pos = string(buffer).find_last_of("\\/");
	return string(buffer).substr(0, pos);
}

/* This function executes SOCJT2 with input and output names. */
void RunSOCJT(string inName, string outName)
{
	string CWD = GetCWD();
	string Path = CWD + "\\SOCJT 2.exe"; // Full path to SOCJT2, assuming SOCJT2 is in the same folder.
	//ShellExecute(NULL, "open", ("\"" + Path + "\"").c_str(), (inName + " " + outName).c_str(), NULL, SW_SHOWDEFAULT);
	system(("\"" + Path + "\" " + inName + " " + outName).c_str());
}
//{ // No idea why any of the below didn't work
//	string CWD = GetCWD();
//	string Path = "\\SOCJT 2.exe"; // Full path to SOCJT2, assuming SOCJT2 is in the same folder.
//	STARTUPINFO si;
//	PROCESS_INFORMATION pi;
//
//	string PathWithArguments = "\"\"" + Path = "\" " + inName + " " + outName + "\"";
//	size_t LengthPWA = PathWithArguments.length();
//	LPSTR CPPathWithArguments = new char[LengthPWA + 1];
//	PathWithArguments._Copy_s(CPPathWithArguments, LengthPWA, LengthPWA);
//	CPPathWithArguments[LengthPWA] = '\0';
//
//	ZeroMemory(&si, sizeof(si));
//	si.cb = sizeof(si);
//	ZeroMemory(&pi, sizeof(pi));
//
//	CreateProcess(Path.c_str(), CPPathWithArguments, NULL, NULL, FALSE, 0, NULL, NULL, &si, &pi);
//
//	CloseHandle(pi.hProcess);
//	CloseHandle(pi.hThread);
//}
//{
//	string CWD = GetCWD();
//	string Path = CWD + "\\SOCJT 2.exe"; // Full path to SOCJT2, assuming SOCJT2 is in the same folder.
//	LPCSTR Parameters = (inName + " " + outName).c_str();
//
//	string PathWithArguments = "\"\"" + Path = "\" " + inName + " " + outName + "\"";
//	size_t LengthPWA = PathWithArguments.length();
//	LPSTR CPPathWithArguments = new char[LengthPWA + 1];
//	PathWithArguments._Copy_s(CPPathWithArguments, LengthPWA, LengthPWA);
//	CPPathWithArguments[LengthPWA] = '\0';
//
//	SHELLEXECUTEINFO ShExecInfo = { 0 };
//	ShExecInfo.cbSize = sizeof(SHELLEXECUTEINFO);
//	ShExecInfo.fMask = SEE_MASK_NOCLOSEPROCESS;
//	ShExecInfo.hwnd = NULL;
//	ShExecInfo.lpVerb = "open";
//	ShExecInfo.lpFile = ("\"" + Path + "\"").c_str(); // "\"SOCJT 2.exe\"";
//	ShExecInfo.lpParameters = NULL;
//	ShExecInfo.lpDirectory = NULL;
//	ShExecInfo.nShow = SW_SHOW;
//	ShExecInfo.hInstApp = NULL;
//	ShellExecuteEx(&ShExecInfo);
//	WaitForSingleObject(ShExecInfo.hProcess, INFINITE);
//	CloseHandle(ShExecInfo.hProcess);
//}

/* This generates an input file using the new fit file. Note that the input file is saved and accessed
from the hard drive so don't run two instances of NFG in the same folder without different output names. */
void GenerateInput(string inName, string tmpinName, string IterationFitFile)
{
	ifstream InputFile(inName.c_str()); // Reads the base input file.
	ofstream NewInputFile(tmpinName.c_str()); // Generates the temporary input file.
	string FitFileLine = "FITFILE = " + IterationFitFile + "\n" + "NFG = True"; // String which has the iteration's fit file to fit, then sets NFG flag in SOCJT to true.

	string SetFitFile = "FITFILE = USENFG"; // Line to be replaced.

	size_t len = SetFitFile.length();
	string strTemp;

	bool useNFG = false; // Used to check if flag is off. There is no reason to not use the flag.

	while (getline(InputFile, strTemp))
	{
		while (true) // Loops through document to replace FITFILE line.
		{
			size_t pos = strTemp.find(SetFitFile);
			if (pos != string::npos)
			{
				strTemp.replace(pos, len, FitFileLine);
				useNFG = true; // Turns off error flag.
			}
			else
			{
				break;
			}
		}
		NewInputFile << strTemp << "\n"; // Outputs temporary input file line by line.
	}
	NewInputFile.close();
	if (useNFG == false) // Incase you forget to set the switch.
	{
		throw 0;
	}
}

void GenerateInputPlus(string inName, string tmpinName, vector<string> OldLine, vector<string> NewLine) // Used for scan input files.
{
	ifstream InputFile(inName.c_str()); // Reads the base input file.
	ofstream NewInputFile(tmpinName.c_str()); // Generates the temporary input file.

	vector<size_t> len; // Length of each string to be searched for.
	for (int i = 0; i < OldLine.size(); i++)
	{
		len.push_back(OldLine[i].length());
	}
	string strTemp;

	while (getline(InputFile, strTemp)) // Reads each line of input file.
	{
		while (true) 
		{
			for (int i = 0; i < OldLine.size(); i++) // Loop over all lines to be replaced
			{
				size_t pos = strTemp.find(OldLine[i]);
				if (pos != string::npos) // If a line to be replaced is found, replace it.
				{
					strTemp.replace(pos, len[i], NewLine[i]); // Replace with new line.
				}
			}
			break;
		}
		NewInputFile << strTemp << "\n"; // Outputs temporary input file line by line.
	}
	NewInputFile.close();
}

inline bool FileExists(const std::string& name) { // Checks if file exists.
	struct stat buffer;
	return (stat(name.c_str(), &buffer) == 0);
}

void NormalFitGenerator(string inFit, string Input, string Output, double StdDev, int N, string strThread)
{
	srand(time(0));

	//string inFit, outFit, Input, Output; // Fit file, Basename for the generated fit files, Input filename, Basename for the generated output files
	//cout << "Enter fit file name:" << endl;
	//getline(cin, inFit);
	//cout << "Enter input file name:" << endl;
	//getline(cin, Input);
	//cout << "Enter output file name or press enter to use " << Input << ".out:" << endl;
	//getline(cin, Output);

	//int N = 0; // Number of iterations.
	//string strThread; // Number of parallel threads to run.
	//double StdDev; // Standard Deviation of fit levels.
	//cout << "Enter Standard Deviation (Normal Distribution):" << endl;
	//cin >> StdDev;
	//cout << "Enter number of iterations:" << endl;
	//cin >> N;
	//cin.ignore();
	//cout << "Enter degree of parallelization or press enter to use default:" << endl;
	//getline(cin, strThread);

	if (Output.empty()) // If nothing is entered, use Input.out name.
	{
		Output = Input + ".out";
	}

	if (!strThread.empty()) // If something was entered, then set that as the number of threads.
	{
		int intThread = atoi(strThread.c_str());
		omp_set_dynamic(0);
		omp_set_num_threads(intThread);
	}

	string outFit = inFit + "_" + Output; // The base name for the fit files generated.

	//string tmpInput = "tmpInput_" + Output + ".tmp"; // Name for the temporary input which SOCJT uses to fit.

	std::ifstream infile(inFit.c_str()); // Checks for fit file.
	if (!infile)
	{
		cout << endl << "Failed to open fit file: " << inFit;
		return;
	}

	/* Definitions to keep track of the speed */
	clock_t start;
	double duration;
	start = clock();

	/* The following stores all numbers in the fit file into ExactLevels */
	std::vector<double> ExactLevels; // Stores the original fit file.
	double val;
	while (infile >> val)
	{
		ExactLevels.push_back(val);
	}

	ofstream TotalOutput(Output + ".total.nfgout"); // This is the output file with each fit value running down a column.
	string Marker = "NFG_OUTPUT"; // The program searches for the line with this string, which is the line I've made to have all the values of interest.

	/* Now we check that the USENFG switch is present */
	string tmp1, tmp2;
	try
	{
		GenerateInput(Input, tmp1, tmp2); // Check to see if "USENFG" switch is present. Only the first string matters.
	}
	catch (int k)
	{
		if (k == 0) // Means NFG switch is off. There is no reason for this to be off.
		{
			TotalOutput << "\n" << "No switch (USENFG) detected in input file." << endl;
			return;
		}
	}
	std::remove(tmp1.c_str());

	int j; // Iteration index.
	int i; // For fit file generation

	/* First we generate all the fit files outside of the parallel loop */
	for (j = 0; j < N; j++)
	{
		/* This converts the index (iteration) to a string for naming purposes*/
		ostringstream jToString;
		jToString << j + 1;

		string IterationFitFile = outFit + "." + jToString.str() + ".nfgfit"; // These are the names of the generated fit files for each iteration

		if (FileExists(IterationFitFile.c_str())) // If fit file already exists, it means we are resuming and don't need to generate fit files.
		{
			break;
		}

		std::ofstream FitFile(IterationFitFile); // This holds the fit file each iteration.

		FitFile << ExactLevels[0] << "\n"; // This generates the fit files. First the number of lines is defined and then the levels are changed by a normal distribution. The format is exactly that of a SOCJT2 fit file.
		for (i = 1; i < ExactLevels.size(); i += 4)
		{
			FitFile << setprecision(10) << rand_normal(ExactLevels[i], StdDev) << "\t" << ExactLevels[i + 1] << "\t" << ExactLevels[i + 2] << "\t" << ExactLevels[i + 3] << "\n"; // Level, j, nj, Sigma
		}
		FitFile.close();

	}

#pragma omp parallel for
	for (j = 0; j < N; j++) // Begin iterations.
	{
		/* This converts the index (iteration) to a string for naming purposes*/
		ostringstream jToString;
		jToString << j + 1;

		string IterationFitFile = outFit + "." + jToString.str() + ".nfgfit"; // These are the names of the generated fit files for each iteration
		string IterationOutFile = Output + "." + jToString.str() + ".nfgout"; // These are the names of the SOCJT2 outputs for each iteration.
		string tmpInput = "tmpInput_" + Output + "." + jToString.str() + ".tmpinput"; // Temporary input file which is read by SOCJT2.

		if (FileExists(IterationOutFile.c_str())) // Checks if the file is already made, moves on if so. Made to continue calculations.
		{
			string strTemp; // Temporary string used to read the file.
			ifstream OutputToRead(IterationOutFile); // Reads the output of the current iteration.
			/* This loop checks each line until Marker is found and puts that line into the total output file */
			while (getline(OutputToRead, strTemp))
			{
				size_t pos = strTemp.find(Marker);
				if (pos != string::npos)
				{
					TotalOutput << jToString.str() << "\t" << strTemp << endl; // Records exactly the line which contains Marker and ends the search.
					break;
				}
			}
			continue;
		}

		//std::ofstream FitFile(IterationFitFile); // This holds the fit file each iteration.

		//FitFile << ExactLevels[0] << "\n"; // This generates the fit files. First the number of lines is defined and then the levels are changed by a normal distribution. The format is exactly that of a SOCJT2 fit file.
		//for (i = 1; i < ExactLevels.size(); i += 4)
		//{
		//	FitFile << setprecision(10) << rand_normal(ExactLevels[i], StdDev) << "\t" << ExactLevels[i + 1] << "\t" << ExactLevels[i + 2] << "\t" << ExactLevels[i + 3] << "\n"; // Level, j, nj, Sigma
		//}
		//FitFile.close();

		GenerateInput(Input, tmpInput, IterationFitFile); // Creates an input file with the name in tmpInput which is the same as the input file but with the fit file modified.

		RunSOCJT(tmpInput, IterationOutFile); // Runs SOCJT2 with the temporary input file and outputs the iteration output file.
		std::remove(tmpInput.c_str()); // Removes temporary input file after SOCJT2 finishes.

		string strTemp; // Temporary string used to read the file.
		ifstream OutputToRead(IterationOutFile); // Reads the output of the current iteration.
		/* This loop checks each line until Marker is found and puts that line into the total output file */
		while (getline(OutputToRead, strTemp))
		{
			size_t pos = strTemp.find(Marker);
			if (pos != string::npos)
			{
				TotalOutput << jToString.str() << "\t" << strTemp << endl; // Records exactly the line which contains Marker and ends the search.
				break;
			}
		}
	} // End Iterations Loop
	//remove(tmpInput.c_str()); // Deletes the temporary input file.

	duration = (clock() - start) / (double)CLOCKS_PER_SEC; // Calculates the time of the process.
	TotalOutput << "\n" << "NFG took " << duration << " seconds.";
	TotalOutput.close();

	return;
}

void RMSGridScan(vector<string> ParameterName, vector<double> ParameterStart, vector<double> ParameterStep, vector<int> ParameterStepNum, string InputName, string OutputName)
{
	int TotalSteps = 1;
	for (int i = 0; i < ParameterStepNum.size(); i++)
	{
		TotalSteps *= ParameterStepNum[i]; // Product of all numbers of steps gives us total number of steps to be taken.
	}

	string* RMSArray = new string[TotalSteps]; // Holds each RMS. I use this to keep values in order.

	vector<string> Replace; // String to replace in the input file
	for (int i = 0; i < ParameterName.size(); i++)
	{
		Replace.push_back(ParameterName[i] + " = SCANP" + to_string(i + 1)); // i.e "MODED = SCANP2" will be replaced. No check is done on fit boolean.
	}

	ofstream tmpTotal(OutputName + "_tmp.total.scan"); // Stores iterations as the process moves, incase of unexpected interuptions.
	tmpTotal << "Input information:" << endl;
	for (int i = 0; i < ParameterName.size(); i++)
	{
		tmpTotal << ParameterName[i] << "\t" << ParameterStart[i] << "\t" << ParameterStep[i] << "\t" << ParameterStepNum[i] << endl;
	}
	tmpTotal << "\n\n";

#pragma omp parallel for
	for (int i = 0; i < TotalSteps; i++) // Loop through all grid coordinates.
	{	
		int* GridIndex = new int[ParameterName.size()]; // This tells us what step we are on for each parameter. An "abacus" if you will.

		/* This generates what each grid index should be based on total index */
		int tmpIndex = i;
		for (int k = ParameterStep.size() - 1; k >= 0; k--)
		{
			int mod = 1;
			for (int kk = 0; kk < k; kk++)
			{
				mod *= ParameterStepNum[kk]; // Size of grid before the k'th parameter, i.e how many steps we have to take to add one more to the k'th index.
			}

			GridIndex[k] = tmpIndex / mod; // Floor division tells us how many times we've passed the k'th dimension, and thus how much we've ++'d.
			tmpIndex = tmpIndex - GridIndex[k] * mod; // Lowers dimension of index so that we can do the same for k-1'th dimension. Essential takes remainder.
		}

		vector<string> NewLine; // Line with new parameter value.
		RMSArray[i] += to_string(i + 1) + "\t"; // Puts step index onto each line.
		for (int j = 0; j < ParameterName.size(); j++)
		{
			ostringstream osst;
			osst << ParameterStart[j] + ParameterStep[j] * (double)GridIndex[j];
			NewLine.push_back(ParameterName[j] + " = " + osst.str()); // String that contains the parameter value.
			RMSArray[i] += osst.str() + "\t"; // Records parameter value into each line.
		}

		string tmpScan = "tmpScan_" + OutputName + "." + to_string(i + 1) + ".in"; // temporary input file name.
		string IterationOutput = OutputName + "." + to_string(i + 1) + ".scan"; // Iteration output file name.
		if (FileExists(IterationOutput.c_str())) // Checks if the file is already made, moves on if so. Made to continue calculations.
		{
			string Marker = "RMS Error ="; // Will search for this string.
			string tmpString;
			ifstream OutputToRead(IterationOutput); // Reads the output of the current iteration.
			/* This loop checks each line until Marker is found and puts that line into the total output file */
			while (getline(OutputToRead, tmpString))
			{
				size_t pos = tmpString.find(Marker);
				if (pos != string::npos)
				{
					tmpString.erase(0, 11); // Deletes the "RMS Error ="
					RMSArray[i] += tmpString; // Records exactly the line which contains Marker and ends the search.
					break;
				}
			}
			tmpTotal << RMSArray[i] << endl;
			continue;
		}

		GenerateInputPlus(InputName, tmpScan, Replace, NewLine); // Generates input with all SCANPi replaced with actual values.

		RunSOCJT(tmpScan, IterationOutput); // Runs SOCJT2 with temp input and iteration output name.
		std::remove(tmpScan.c_str()); // Deletes temp input

		string Marker = "RMS Error ="; // Will search for this string.
		string tmpString;
		ifstream OutputToRead(IterationOutput); // Reads the output of the current iteration.
		/* This loop checks each line until Marker is found and puts that line into the total output file */
		while (getline(OutputToRead, tmpString))
		{
			size_t pos = tmpString.find(Marker);
			if (pos != string::npos)
			{
				tmpString.erase(0, 11); // Deletes the "RMS Error ="
				RMSArray[i] += tmpString; // Records exactly the line which contains Marker and ends the search.
				break;
			}
		}

		tmpTotal << RMSArray[i] << endl;
		delete[] GridIndex;
	} // End iterations

	string ScanName = OutputName + ".total.scan";
	ofstream ScanTotal(ScanName);

	//ScanTotal << "Input information:" << endl;
	//for (int i = 0; i < ParameterName.size(); i++)
	//{
	//	ScanTotal << ParameterName[i] << "\t" << ParameterStart[i] << "\t" << ParameterStep[i] << "\t" << ParameterStepNum[i] << endl;
	//}
	//ScanTotal << "\n\n";

	ScanTotal << "Step" << "\t";
	for (int i = 0; i < ParameterName.size(); i++)
	{
		ScanTotal << ParameterName[i] << "\t";
	}
	ScanTotal << "RMS Error" << endl;
	for (int i = 0; i < TotalSteps; i++)
	{
		ScanTotal << RMSArray[i] << endl;
	}
	std::remove((OutputName + "_tmp.total.scan").c_str());
	delete[] RMSArray;
}

int main()
{
	int Option;
	cout << "_______________________________________________" << endl;
	cout << "[                                             ]" << endl;
	cout << "[           SOCJT_2 Utility Program           ]" << endl;
	cout << "[                                             ]" << endl;
	cout << "[                                             ]" << endl;
	cout << "[ 1. Normal Fit Generator                     ]" << endl;
	cout << "[ 2. RMS Grid Scan                            ]" << endl;
	cout << "[_____________________________________________]" << endl;
	cin >> Option;
	cin.ignore();

	if (Option == 1)
	{
		string inFit, outFit, Input, Output; // Fit file, Basename for the generated fit files, Input filename, Basename for the generated output files
		cout << "Enter fit file name:" << endl;
		getline(cin, inFit);
		cout << "Enter input file name:" << endl;
		getline(cin, Input);
		cout << "Enter output file name or press enter to use " << Input << ".out:" << endl;
		getline(cin, Output);

		int N = 0; // Number of iterations.
		string strThread; // Number of parallel threads to run.
		double StdDev; // Standard Deviation of fit levels.
		cout << "Enter Standard Deviation (Normal Distribution):" << endl;
		cin >> StdDev;
		cout << "Enter number of iterations:" << endl;
		cin >> N;
		cin.ignore();
		cout << "Enter degree of parallelization or press enter to use default:" << endl;
		getline(cin, strThread);

		NormalFitGenerator(inFit, Input, Output, StdDev, N, strThread);
	}
	if (Option == 2)
	{
		string InputName, OutputName, strThread;
		vector<string> ParameterName;
		vector<double> ParameterStart, ParameterStep;
		vector<int> ParameterStepNum;
		
		cout << "Description: This utility will scan desired parameters and create a grid of RMS fit errors at each parameter step." << "\n" << endl;
		cout << "Instructions: Before running this program, create an input file where each parameter to be scanned is set equal to the string \"SCANPi\" "
			<< "where i denotes the order of the parameter to be scanned. Be sure to set the fit boolean to false. When prompted to enter parameter names, enter the string before "
			<< "the \"=\", i.e. \"MODED\"." << "\n" << endl;
		cout << "Enter input file name:" << endl;
		getline(cin, InputName);
		cout << "Enter output file name or press enter to use " << InputName << ".out" << endl;
		getline(cin, OutputName);
		cout << "Enter degree of parallelization or press enter to use default:" << endl;
		getline(cin, strThread);
		if (OutputName.empty()) // If nothing is entered, use Input.out name.
		{
			OutputName = InputName + ".out";
		}
		if (!strThread.empty()) // If something was entered, then set that as the number of threads.
		{
			int intThread = atoi(strThread.c_str());
			omp_set_dynamic(0);
			omp_set_num_threads(intThread);
		}

		string tmpString;
		cout << "Enter the names of the parameters in order. Enter an empty line to finish:" << endl;
		while (getline(cin, tmpString) && !tmpString.empty())
		{
			ParameterName.push_back(tmpString);
		}

		cout << "Enter the starting values of your parameters in order:" << endl;
		while (getline(cin, tmpString))
		{
			ParameterStart.push_back(stod(tmpString));
			if (ParameterStart.size() == ParameterName.size()) break;
		}

		cout << "Enter the step size of each parameter in order:" << endl;
		while (getline(cin, tmpString))
		{
			ParameterStep.push_back(stod(tmpString));
			if (ParameterStep.size() == ParameterName.size()) break;
		}
		cout << "Enter the number of steps in each parameter in order:" << endl;
		while (getline(cin, tmpString))
		{
			ParameterStepNum.push_back(stoi(tmpString));
			if (ParameterStepNum.size() == ParameterName.size()) break;
		}
		
		RMSGridScan(ParameterName, ParameterStart, ParameterStep, ParameterStepNum, InputName, OutputName);
	}
	else
	{
		cout << "No option selected: Terminating Program." << endl;
	}

	return 0;
}
