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

int main()
{
	srand(time(0));

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

	outFit = inFit + "_" + Output; // The base name for the fit files generated.

	//string tmpInput = "tmpInput_" + Output + ".tmp"; // Name for the temporary input which SOCJT uses to fit.

	std::ifstream infile(inFit.c_str()); // Checks for fit file.
	if (!infile)
	{
		cout << endl << "Failed to open fit file: " << inFit;
		return 0;
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
			return 0;
		}
	}
	remove(tmp1.c_str());

	int j; // Iteration index.
	int i; // For fit file generation

	/* First we generate all the fit files outside of the parallel loop */
	for (j = 0; j < N; j++)
	{
		/* This converts the index (iteration) to a string for naming purposes*/
		ostringstream jToString;
		jToString << j + 1;

		string IterationFitFile = outFit + "." + jToString.str() + ".nfgfit"; // These are the names of the generated fit files for each iteration

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

		//std::ofstream FitFile(IterationFitFile); // This holds the fit file each iteration.

		//FitFile << ExactLevels[0] << "\n"; // This generates the fit files. First the number of lines is defined and then the levels are changed by a normal distribution. The format is exactly that of a SOCJT2 fit file.
		//for (i = 1; i < ExactLevels.size(); i += 4)
		//{
		//	FitFile << setprecision(10) << rand_normal(ExactLevels[i], StdDev) << "\t" << ExactLevels[i + 1] << "\t" << ExactLevels[i + 2] << "\t" << ExactLevels[i + 3] << "\n"; // Level, j, nj, Sigma
		//}
		//FitFile.close();

		GenerateInput(Input, tmpInput, IterationFitFile); // Creates an input file with the name in tmpInput which is the same as the input file but with the fit file modified.

		RunSOCJT(tmpInput, IterationOutFile); // Runs SOCJT2 with the temporary input file and outputs the iteration output file.
		remove(tmpInput.c_str()); // Removes temporary input file after SOCJT2 finishes.

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

	return 0;
}
