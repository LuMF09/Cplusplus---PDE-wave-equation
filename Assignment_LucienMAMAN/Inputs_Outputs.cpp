//Inputs_Outputs.cpp
/*
Author: Lucien MAMAN
Date: 05/12/2016

Description: Implementation of
*/
#include "Inputs_Outputs.h"
#include "Errors.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

//FUNCTION which ask inputs to the user, check if they are conforms and return a vector with correct parameters

std::vector<double> Inputs_Outputs::askImputs()
{
	std::vector<double>result(0);
	std::string choiceFCT(""), gridUser(""), choiceScheme(""), quickChoiceScheme("");
	double t(0.00), cflUser(0.00);
	bool isCorrect(false);

	//Choose the function to work with
	std::cout << "Choose the function to work with:" << std::endl <<
				 "Enter 0 to work with f0=1/2 * (sign(x)+1)" << std::endl << 
				 "Enter 1 to work with f1=1/2*(exp(-x^2))" << std::endl;
	std::cin >> choiceFCT;

	//Force the user to enter a correct value
	while (choiceFCT != "0" && choiceFCT != "1")
	{
		std::cerr << std::endl << "Please choose a correct function:" << std::endl << 
								  "Enter 0 to work with f0=1/2 * (sign(x)+1)" << std::endl << 
								  "Enter 1 to work with f1=1/2*(exp(-x^2))" << std::endl;
		std::cin >> choiceFCT;
	}
	result.push_back(stoi(choiceFCT)); // convert from string to int

	//Choose the time to run
	std::cout << std::endl << "Choose the time t to run:" << std::endl;
	std::cin >> t;

	//Force the user to enter a correct value
	while (t <= 0)
	{
		std::cerr << std::endl << "Please choose a time step >0" << std::endl;
		std::cin >> t;
	}
	result.push_back(t);

	//Choose the number of point in the grid
	std::cout << std::endl << "Choose the number of points in the grid:" << std::endl;
	std::cin >> gridUser;
	if (stoi(gridUser) == stod(gridUser) && stoi(gridUser) > 1) // test if it is an integer and if it is superior at 1
	{
		isCorrect = true;
		result.push_back(stoi(gridUser));
	}

		while (!isCorrect) //test if it is and integer or not
		{
			std::cerr << std::endl << "Please choose a correct number of point in your grid" << std::endl;
			std::cin >> gridUser;
			if (stoi(gridUser) == stod(gridUser) && stoi(gridUser)>0)
			{
				isCorrect = true;
				result.push_back(stoi(gridUser));
			}
		}

	//Force the user to enter a correct value
	std::cout << std::endl << "Choose the scheme:" << std::endl << 
							  "1 for the UPWIND EXPLICIT scheme" << std::endl << 
							  "2 for the UPWIND IMPLICIT scheme" << std::endl << 
							  "3 for the LAX WENDROFF scheme" << std::endl << 
							  "4 for the RICHARD-MAYER MULTISTEP scheme" << std::endl;
	std::cin >> choiceScheme;
	while (choiceScheme != "1" && choiceScheme != "2" && choiceScheme != "3" && choiceScheme != "4")
	{
		std::cerr << std::endl << "Please choose a correct number to have the appropriate scheme" << std::endl;
		std::cin >> choiceScheme;
	}
	result.push_back(stoi(choiceScheme)); //convert from string to integer

	//Choose the cfl
	std::cout << std::endl << "Choose the CFL number:" << std::endl;
	std::cin >> cflUser;

	//Force the user to enter a correct value
	while (cflUser <= 0)
	{
		std::cerr << std::endl << "Please choose a CFL >0" << std::endl;
		std::cin >> cflUser;
	}
	result.push_back(cflUser);

	//Choose the way to compute imlplicit scheme
	if (result[3] == 2)
	{
		std::cout << std::endl << "Two different way to compute:" << std::endl <<
			"1 for a fast computation" << std::endl <<
			"2 for a long computation using the Thomas decomposition" << std::endl;

		std::cin >> quickChoiceScheme;
		while (quickChoiceScheme != "1" && quickChoiceScheme != "2")
		{
			std::cerr << std::endl << "Please choose a correct number to have the correct way to compute" << std::endl;
			std::cin >> quickChoiceScheme;
		}
		result.push_back(stoi(quickChoiceScheme));
	}
	if (result.size()<4)
	{
		throw TooLittleTab("Bad initialisation in askInputs function");
	}
	return result;
}

void Inputs_Outputs::writeVectorInCSV(std::vector<double> &v, int nbPoint, double step, std::string a)
{
	std::ofstream fichier("Results.csv", std::ios::out | std::ios::app);  // write at the end of the file
	std::string vector(a);
	if (fichier)
	{
		if (vector == "POSITION")
		{
			fichier << ";";
		}
		else fichier << "t=" << step <<" "<< vector << " ;"; 

		for (int i = 0; i < nbPoint; i++)
		{
			fichier << v[i] << ";"; // write vector in the Results CSV file
		}

		fichier << std::endl;
		fichier.close();
	}
	else throw DocOpen("Results.csv");
}

void Inputs_Outputs::writeNumericalInCSV(std::vector<std::vector<double>> &v, std::string a, double step, int t, int size)
{
	std::string scheme(a);
	std::ofstream fichier("Results.csv", std::ios::out | std::ios::app);  // write at the end of the file

	if (fichier)
	{
		fichier << scheme << std::endl;
		//Write in the CSV file the last line of the numerical solution
		fichier << "t= " << step << " NUMERICAL "<< scheme << ";" << v[t][0] << "; ";//first element of the last line

		for (int i = 1; i < size - 1; i++)
		{
			fichier << v[t][i] << ";";
		}
		fichier << v[t][size - 1]; // last element of the last line
		fichier << std::endl << std::endl;

		fichier.close();
	}
	else throw DocOpen("Results.csv");
}

void Inputs_Outputs::writeNormInCSV(double step, double norm0, double norm2Normalized)
{
	std::ofstream fichierNorm("Norms.csv", std::ios::out | std::ios::app);  // write at the end of the file
	if (fichierNorm)
	{
		fichierNorm << ";" << "Norm 0" << "; " << "Norm 2" << "; " << std::endl;
		fichierNorm << "t=" << step << ";";
		fichierNorm << norm0 << ";";
		fichierNorm << norm2Normalized << ";";
		
		fichierNorm.close();
	}
	else throw DocOpen("Norms.csv");
}