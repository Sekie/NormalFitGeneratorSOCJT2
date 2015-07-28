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

/* This function calls SOCJT2. Note that Terrance doesn't leave the input and output file
names as command line inputs, so I had to modify the SOCJT2 source code to do so. Therefore,
the typical "SOCJT_2.exe" will not work. I renamed my modification "NFGSOCJT2.exe" so make
sure to have this modified SOCJT executable when running NFG. This function executes SOCJT2
with input and output names. SOCJT also needs to be modified to display all the fit parameters
in one line after the Marker string. */
void RunSOCJT(string inName, string outName)
{
	// system(("./NFGSOCJT2 " + inName + " " + outName).c_str()); // For bash
	system(("\"SOCJT 2\".exe " + inName + " " + outName).c_str()); // For CMD // I wish I could figure out a way to avoid system and to make this compatible for both Windows and Linux.
}

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

	int N = 0;
	double StdDev;
	cout << "Enter Standard Deviation (Normal Distribution):" << endl;
	cin >> StdDev;
	cout << "Enter number of iterations:" << endl;
	cin >> N;

	if (Output.empty() == 1) // If nothing is entered, use Input.out name.
	{
		Output = Input + ".out";
	}

	outFit = inFit + "_" + Output; // The base name for the fit files generated.

	string tmpInput = "tmpInput_" + Output + ".tmp"; // Name for the temporary input which SOCJT uses to fit.

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

	for (int j = 0; j < N; j++) // Begin iterations.
	{
		/* This converts the index (iteration) to a string for naming purposes*/
		ostringstream jToString;
		jToString << j + 1;

		string IterationFitFile = outFit + "." + jToString.str() + ".nfgfit"; // These are the names of the generated fit files for each iteration
		string IterationOutFile = Output + "." + jToString.str() + ".nfgout"; // These are the names of the SOCJT2 outputs for each iteration.

		std::ofstream FitFile(IterationFitFile); // This holds the fit file each iteration.

		FitFile << ExactLevels[0] << "\n"; // This generates the fit files. First the number of lines is defined and then the levels are changed by a normal distribution. The format is exactly that of a SOCJT2 fit file.
		for (int i = 1; i < ExactLevels.size(); i += 4)
		{
			FitFile << setprecision(10) << rand_normal(ExactLevels[i], StdDev) << "\t" << ExactLevels[i + 1] << "\t" << ExactLevels[i + 2] << "\t" << ExactLevels[i + 3] << "\n"; // Level, j, nj, Sigma
		}
		FitFile.close();

		try
		{
			GenerateInput(Input, tmpInput, IterationFitFile); // Creates an input file with the name in tmpInput which is the same as the input file but with the fit file modified.
		}
		catch (int k)
		{
			if (k == 0) // Means NFG switch is off. There is no reason for this to be off.
			{
				TotalOutput << "\n" << "No switch (USENFG) detected in input file." << endl;
				return 0;
			}
		}

		RunSOCJT(tmpInput, IterationOutFile); // Runs SOCJT2 with the temporary input file and outputs the iteration output file.

		string strTemp; // Temporary string used to read the file.
		ifstream OutputToRead(IterationOutFile); // Reads the output of the current iteration.
		/* This loop checks each line until Marker is found and puts that line into the total output file */
		while (getline(OutputToRead, strTemp))
		{
			size_t pos = strTemp.find(Marker);
			if (pos != string::npos)
			{
				TotalOutput << strTemp << endl; // Records exactly the line which contains Marker and ends the search.
				break;
			}
		}
	} // End Iterations Loop
	remove(tmpInput.c_str()); // Deletes the temporary input file.

	duration = (clock() - start) / (double)CLOCKS_PER_SEC; // Calculates the time of the process.
	TotalOutput << "\n" << "NFG took " << duration << " seconds.";
	TotalOutput.close();

	return 0;
}
