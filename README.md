# GrowthRates

C implementation of the codes described in "Challenges and pitfalls of inferring microbial growth rates from lab cultures", by Ana-Hermina Ghenu, Loïc Marrec and Claudia Bank, as well as additional figures and tables that report the fits and estimates obtained from the analysis of four data sets.

Briefly, we perform stochastic simulations of population dynamics using different growth patterns (Baranyi, Gompertz, Logistic and Richards).

First define your parameters directly in the .c file. Then, to compile C code, open a terminal, change the working directory to where the code is, compile the code using the command g++ followed by the file name of the code (e.g., g++ Baranyi.c) and run the code using the executable (e.g., ./ExecutableFileName.out). A .txt file with the results will be generated where each column corresponds to a stochastic realization and each row to the population size recorded at the time defined by the user.

The source code is freely available under the GNU GPLv3 license.

If you find this code useful for your research, please cite the associated reference, "Challenges and pitfalls of inferring microbial growth rates from lab cultures", by Ana-Hermina Ghenu, Loïc Marrec and Claudia Bank.
