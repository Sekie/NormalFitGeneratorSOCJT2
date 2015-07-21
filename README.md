NormalFitGeneratorSOCJT2 is a script written to utilize SOCJT2's vibronic fitting capabilities. The SOCJT2 program was developed by Terrance J. Codd (https://github.com/tcodd86) for the computation and fitting of vibronic spectra. Pure calculations could be fit using known parameters which left an open question of how experimental error could affect the parameters of the fit, if we knew the right parameters but let the levels vary by an experimental error. NormalFitGeneratorSOCJT2 was written for that purpose. 

NormalFitGeneratorSOCJT2 takes a set of vibronic levels (in the fit file format) and varies the levels by a normal distribution, outputing a fit file using the exact same assignments but with varied levels. Then SOCJT2 is executed and fits these levels giving an output. This is repeated as many times as the user wishes.

The input files required are the fit file, same as described by SOCJT2 and intended to be a set of pure calculations, and the input file, same as described by SOCJT2 but with the addition that the FITFILE is set to USENFG: "FITFILE = USENFG".

The user must also input the number of iterations and the standard deviation on the normal distribution.

A modified console application for SOCJT2 will be necessary to run this program and should be named "NFGSOCJT2.exe." The modified application should allow for command line inputs, which the original did not. 

