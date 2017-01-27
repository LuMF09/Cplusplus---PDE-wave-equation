//main.cpp
/*
Author: Lucien MAMAN
Date: 05/12/2016

Description: Main of the project
*/
#include "Schemes.h"
#include "Inputs_Outputs.h"
#include "Errors.h"
#include <iostream>

int main()
{
	int choice(0), grid(0), scheme(0), quickChoice(0);
	double step(0.00), cflU(0.00);
	std::vector<double> results;
	try
	{
		results = Inputs_Outputs::askImputs(); //vector results composed by corrects type of arguments

		//Recuperation of data
		choice = (int)results[0];
		step = results[1];
		grid = (int)results[2];
		scheme = (int)results[3];
		cflU = results[4];

		if (results.size()>5)
		{
			quickChoice = (int)results[5];
		}

		switch (scheme) // Depending the choice of the user, creation of the appropriate object and realization of simulations
		{
		case 1:
		{
				  ExplicitUpWind scheme1(choice, step, grid,cflU);
				  scheme1.showScheme();
		}
			break;

		case 2:
		{
				  ImplicitUpWind scheme1(choice, quickChoice, step, grid, cflU);
				  scheme1.showScheme();
		}
			break;
		case 3:
		{
				  Lax_Wendroff scheme1(choice, step, grid, cflU);
				  scheme1.showScheme();
		}
			break;
		case 4:
		{
				  RichtmyerMS scheme1(choice, step, grid, cflU);
				  scheme1.showScheme();
		}
			break;
		}
	}
	catch (std::exception &e)
	{
		std::cerr << e.what() << std::endl;
	}
	
		return 0;
	
}