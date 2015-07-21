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
// #include <ctime>
#include <locale>

using namespace std;

/* This program takes a SOCJT2 fit file and outputs a fit file using the same assignments
where the levels are varied using a normal distribution with standard deviation hard
coded below. This was made to be used to create varying "experiment" data sets from
a pure calculation and run fitting tests against such.

USAGE:
- Make sure to have the modified SOCJT 2 executable named "NFGSOCJT2.exe" in the same folder
- The fit file is formatted exactly as for SOCJT 2, and each iteration is saved as
(fit output name).(iteration index).nfgfit
- The output files are regular output files and are named
(output name).(iteration index).nfgout
- In the input file, set the fit file to "USENFG"									*/

/* This is the Box Muller Method for generating random numbers following a normal distribution
where the two arguments are the mean or expectation, and the standard deviation */
double rand_normal(double mean, double stddev) //Box muller method
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

/* This function calls SOCJT2. Note that Terrance doesn't leave the input and output file
names as command line inputs, so I had to modify the SOCJT2 source code to do so. Therefore,
the typical "SOCJT_2.exe" will not work. I renamed my modification "NFGSOCJT2.exe" so make
sure to have this modified SOCJT executable when running NFG */
void RunSOCJT(string inName, string outName, string Iteration)
{
	string NFGName = outName + "." + Iteration + ".nfgout";
	// system(("./NFGSOCJT2 " + inName + " " + NFGName).c_str()); // For bash
	system(("NFGSOCJT2.exe " + inName + " " + NFGName).c_str()); // For CMD
}

/* This generates an input file using the new fit file. Note that the input file is saved and accessed
from the hard drive so don't run two instances of NFG in the same folder. */
void GenerateInput(string inName, string IterationFitFile)
{
	ifstream InputFile(inName.c_str());
	ofstream NewInputFile("tmpInput");
	string FitFileLine = "FITFILE = " + IterationFitFile;

	string SetFitFile = "FITFILE = USENFG";

	size_t len = SetFitFile.length();
	string strTemp;

	while (getline(InputFile, strTemp))
	{
		while (true)
		{
			size_t pos = strTemp.find(SetFitFile);
			if (pos != string::npos)
			{
				strTemp.replace(pos, len, FitFileLine);
			}
			else
			{
				break;
			}
		}
		NewInputFile << strTemp << "\n";
	}
	NewInputFile.close();
}

int main()
{
	srand(time(0));

	string inFit, outFit, Input, Output; // Original fit file, randomized fit file, SOCJT2 input, name for output
	cout << "Name of original fit file?" << endl;
	getline(cin, inFit);
	cout << "Name of randomized fit files? Enter to use the same name as the input." << endl;
	getline(cin, outFit);
	if (outFit == " ");
	{
		outFit = inFit;
	}
	cout << "Name of input file?" << endl;
	getline(cin, Input);
	cout << "Name of output file?" << endl;
	getline(cin, Output);
	if (Output == " ")
	{
		Output = Input + ".out";
	}

	int N = 0;
	double StdDev;
	cout << "How many iterations?" << endl;
	cin >> N;
	cout << "Standard Deviation? (Normal Distribution)" << endl;
	cin >> StdDev;

	std::ifstream infile(inFit.c_str());

	if (!infile)
	{
		cout << endl << "Failed to open fit file: " << inFit;
		return 0;
	}

	//	clock_t start;
	//	double duration;

	//	start = clock();

	std::vector<double> ExactLevels;
	double val;
	while (infile >> val)
	{
		ExactLevels.push_back(val);
	}

	for (int j = 0; j < N; j++)
	{
		ostringstream jToString;
		jToString << j;
		string IterationFitFile = outFit;
		IterationFitFile.append(".");
		IterationFitFile.append(jToString.str());
		IterationFitFile.append(".nfgfit");

		//  string IterationDifference = output;
		//  IterationDifference.append(".dif.");
		//  IterationDifference.append(jToString.str());
		//  IterationDifference.append(".dif");

		string strItIndex = jToString.str();

		std::ofstream FitFile(IterationFitFile.c_str());
		//  std::ofstream DifFile(IterationDifference.c_str());

		FitFile << ExactLevels[0] << "\n"; // This generates the fit files. First the number of lines is defined and then the levels are changed by a normal distribution.
		for (int i = 1; i < ExactLevels.size(); i += 4)
		{
			FitFile << setprecision(10) << rand_normal(ExactLevels[i], StdDev) << "\t" << ExactLevels[i + 1] << "\t" << ExactLevels[i + 2] << "\t" << ExactLevels[i + 3] << "\n";
		}
		FitFile.close();

		GenerateInput(Input, IterationFitFile);

		RunSOCJT("tmpInput", Output, strItIndex);
	} // End Iterations Loop
	remove("tmpInput");

	// duration = (clock() - start) / (double)CLOCKS_PER_SEC;

	//	cout << "This took " << duration << " seconds."
	return 0;
}