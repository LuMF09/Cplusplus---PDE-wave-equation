//Schemes.cpp
/*
Author: Lucien MAMAN
Date: 05/12/2016

Description: Implementation of all classes of Schemes.h file:
abstract Schemes,
Upwind Explicit, Upwind Explicit, Lax-Wendroff, Richtmyer multi-step all inherited from Schemes
*/
#include "Schemes.h"
#include "Inputs_Outputs.h"
#include "Maths_Functions.h"
#include "Errors.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

double const LEFT_BOUNDARY(-50.00), RIGHT_BOUNDARY(50.00), U(1.75); //x belong to [-50;50] and u=1.75

//SCHEMES ABSTRACT CLASS

//CONSTRUCTOR

//Initialization of all protected variables
Schemes::Schemes(int choice, double step, int gridSize, double cfl1) : x(gridSize), vANAL(gridSize), vNUM(), choice(choice), step(step), deltaT(0), deltaX(0), cfl(0), nbIt(0)
{
	if (gridSize == 1) //Exception division by 0
	{
		throw DivideZero("in Scheme constructor (see value of gridSize)");
	}

		deltaX = abs(RIGHT_BOUNDARY - LEFT_BOUNDARY) / (gridSize - 1);

		if (U == 0) //Exception division by 0
		{
			throw DivideZero(" in Scheme constructor (see value of U)");
		}

		deltaT = (cfl1 * deltaX)/U; 

		if (deltaT==0) //Exception division by 0
		{
			throw DivideZero(" in Scheme constructor (see value of delta T)");
		}

		nbIt = (int)std::round(step/deltaT);//round the value to have the closest optimised CFL

		if (nbIt==0) //Exception division by 0
		{
			throw DivideZero(" in Scheme constructor (see value of nbIt)");
		}
		deltaT = step / nbIt;

		if (deltaX == 0) //Exception division by 0
		{
			throw DivideZero(" in Scheme constructor (see value of deltaX)");
		}
		cfl=U*(deltaT / deltaX);

		for (int i = 0; i < nbIt; i++)
		{
			vNUM.push_back(std::vector<double>(gridSize));//create exactly the good size for vNUM
		}
			

	//Construction of the vector x with all positions

	x[0] = LEFT_BOUNDARY;
	x[(*this).getSizeGrid() - 1] = RIGHT_BOUNDARY;
	
	for (int i = 1; i < (*this).getSizeGrid() - 1; i++)
	{
		x[i] = x[i - 1] + deltaX;	
	}

	//Write in the Results CSV file the vector of positions
	Inputs_Outputs::writeVectorInCSV(x, (*this).getSizeGrid(), (*this).getStep(), "POSITION");
}

//ACCESSORS

double Schemes::getChoice() const
{
	return choice;
}
int Schemes::getSizeGrid() const
{
	return x.size();
}
double Schemes::getStep() const
{
	return step;
}

//FUNCTIONS 

void Schemes::calculAnalytical() // Compute and write in a CSV file the analytical solution
{
		std::cout << std::endl << "VALUES of f" << choice << " (Analytical solution) are calculated..." << std::endl;

		//Boundary conditions

		if (choice == 0) //f0 = 1/2 * (sign(x-Ut)+1)
		{
			vANAL[0] = 0.00;
			vANAL[(*this).getSizeGrid() - 1] = 1.00;
		}
		else if (choice == 1) //f1 = 1/2 * exp(-x^2)
		{
			vANAL[0] = 0.00;
			vANAL[(*this).getSizeGrid() - 1] = 0.00;
		}
		else throw BadInput(" in calculAnalytical(), see value of choice"); //Exception BadInput

		// Display the last line of the analytical solution in the CSV file Results

		for (int j = 1; j < (*this).getSizeGrid() - 1; j++){
			if (choice == 0)
			{
				vANAL[j] = 0.5*(Maths_Functions::sign(x[j] - U*((nbIt-1)*deltaT)) + 1.00);
			}
			else if(choice == 1)
			{
				vANAL[j] = 0.5*exp(-((x[j] - U*(*this).getStep())*(x[j] - U*(nbIt - 1)*deltaT)));
			}
			else throw BadInput(" in calculAnalytical(), see value of choice"); //Exception BadInput
		}
		Inputs_Outputs::writeVectorInCSV(vANAL, (*this).getSizeGrid(), (*this).getStep(), "ANALYTICAL");
	}

void Schemes::calculNorms() // Compute norm 0 (max error) and norm 2 (square root of the sum of the square of errors)
{
	double norm2(0.00), norm0(0.00), sumCarre(0.00), norm2Normalized(0.00);

		std::cout << "Norms are calculated ... " << std::endl;

		for (int i = 0; i < (*this).getSizeGrid(); i++)
		{
			sumCarre += (vANAL[i] - vNUM[nbIt-1][i])*(vANAL[i] - vNUM[nbIt-1][i]); // Sum of (error)^2

			if (norm0 < fabs(vANAL[i] - vNUM[nbIt-1][i])) // compare if value of previous max is higher
			{
				norm0 = fabs(vANAL[i] - vNUM[nbIt-1][i]);
			}
		}
		norm2 = sqrt(sumCarre); // Square root of sumCarre
		norm2Normalized = norm2 / (*this).getSizeGrid(); // we divide by the number of point to normalise the norm

		Inputs_Outputs::writeNormInCSV((*this).getStep(), norm0, norm2Normalized); // Write norms in "Norms" CSV file
}

void Schemes::initNumerical(std::vector<std::vector <double>>&v) // initialize a vector of vector with boundaries depending the choice made by the user
{
	
	for (int t = 0; t < nbIt; t++)
	{
		if (choice == 0)
		{
			v[t][0] = 0.00;
			v[t][(*this).getSizeGrid() - 1] = 1.00;
		}
		else if (choice == 1)
		{
			v[t][0] = 0.00;
			v[t][(*this).getSizeGrid() - 1] = 0.00;
		}
		else throw BadInput(" in initNumerical(), see value of choice"); //Exception BadInput
	}

	//Construct v from 1 to v[size of the grid - 2]
	for (int j = 1; j < (*this).getSizeGrid()-1; j++)
	{
		if (choice == 0)
		{
			v[0][j] = 0.5*(Maths_Functions::sign(x[j]) + 1.00);
		}
		else if(choice == 1)
		{
			v[0][j] = 0.5*exp(-(x[j] * x[j]));
		}
		else throw BadInput(" in initNumerical(), see value of choice");//Exception BadInput
	}
}

void Schemes::showScheme() {
	std::cout << std::endl << "PARAMETERS:" << std::endl << std::endl <<
		"Delta T= " << deltaT << std::endl <<
		"Delta X= " << deltaX << std::endl <<
		"Optimised CFL= " << cfl << std::endl << std::endl;

	//test if the scheme is stable
	(*this).isStable();
	std::cout << std::endl;

	//ANALYTICAL SOLYUTION
	(*this).calculAnalytical();

	//NUMERICAL SOLUTION
	(*this).calculNumerical();

	//NORMS
	(*this).calculNorms();

	

	std::cout << std::endl<<"Please find Analytical and Numerical solutions in Results CSV file" << std::endl <<
						    "Please find norms in Norms CSV file"<<std::endl;

	std::cout << std::endl;
}

	//EXPLICIT UpWind Scheme

	//CONSTRUCTOR

ExplicitUpWind::ExplicitUpWind(int choice, double step, int size, double cfl1) : Schemes(choice, step, size, cfl1) //The user choose f0 or f1, the time to run and the size of the grid
	{
	std::cout << std::endl << "The UPWIND scheme : forward in time and backward in space" << std::endl; //Present the scheme 
	}

	//FUNCTIONS

void ExplicitUpWind::calculNumerical() //Compute the NUMERICAL solution of the Explicit UPWIND scheme
{
		(*this).initNumerical(vNUM); //Initialize the vNUM vector of vector

			std::cout << "VALUES of f" << choice << " (Numerical solution) are calculated..." << std::endl;

			for (int t = 0; t <nbIt-1; t++)
			{
				for (int space = 1; space < (*this).getSizeGrid()-1; space++)
				{
					vNUM[t+1][space] = vNUM[t][space] - cfl*(vNUM[t][space] - vNUM[t][space - 1]); //Construction of vNUM following the scheme
				}
			}
			
			//Write in the CSV file the last line of the numerical solution
			Inputs_Outputs::writeNumericalInCSV(vNUM, "EXPLICIT UPWIND", (*this).getStep(),nbIt-1, (*this).getSizeGrid());
	}

void ExplicitUpWind::isStable() //Test if the scheme is stable 
{
	if (cfl<=1) // Stability condition
	{
		std::cout << std::endl << "This system is stable" << std::endl;
	}
	else std::cout << std::endl << "This system is not stable" << std::endl;
}

//IMPLICIT UPWIND SCHEME

//CONSTRUCTOR

ImplicitUpWind::ImplicitUpWind(int choice, int quickChoice, double step, int size, double cfl1) : Schemes(choice, step, size,cfl1), quickChoice(quickChoice) //The user choose f0 or f1, the time to run and the size of the grid and the way to compute
{
	std::cout << std::endl << "The IMPLICIT UPWIND scheme : forward in time and backward in space" << std::endl; // Present the scheme
}

//FUNCTIONS

void ImplicitUpWind::calculNumerical() //Compute the NUMERICAL solution of the Implicit UPWIND scheme
{
	(*this).initNumerical(vNUM); //Initialize the vNUM vector of vector

		std::cout << "VALUES of f" << choice << " (Numerical solution) are calculated..."<< std::endl;

		if (quickChoice==1) // Most efficient computations without matrix decomposition
		{
			for (int t = 0; t < nbIt-1; t++)
			{
				for (int space = 1; space < (*this).getSizeGrid() - 1; space++)
				{
					vNUM[t + 1][space] = (vNUM[t][space]) / (1 + cfl) + (cfl*(vNUM[t + 1][space - 1])) / (1 + cfl); //Construction of vNUM following the scheme
				}
			}
		}
		else if (quickChoice == 2)//THOMAS DECOMPOSIITON
		{
			std::vector<std::vector<double>> matrix((*this).getSizeGrid(), std::vector<double>((*this).getSizeGrid()));// matrix grid*grid
			std::vector<std::vector<double>> vIn;

			for (int i = 0; i < nbIt; i++)
			{
				vIn.push_back(std::vector<double>((*this).getSizeGrid()));//create exactly the good size for vNUM
			}

			(*this).initNumerical(vIn);

			double m(0);
	
			//construction of the matrix
			for (int raw = 0; raw < (*this).getSizeGrid(); raw++)
			{
				for (int col = 0; col < (*this).getSizeGrid(); col++)
				{
					if (raw == col)
					{
						matrix[raw][col] = 1 + cfl;
					}
					else if (raw == col + 1)
					{
						matrix[raw][col] = -cfl;
					}
					else
					{
						matrix[raw][col] = 0;
					}
					std::cout << matrix[raw][col] << " ";
				}			
				std::cout << std::endl;
			}

			for (int time = 1; time < nbIt; time++)
			{
				for (int k = 1; k < (*this).getSizeGrid()-1; k++)
				{
					vIn[time][k] = vIn[time - 1][k];
					m = matrix[k][k - 1] / matrix[k - 1][k - 1];// ak/bk-1 
					vIn[time][k] = vIn[time][k] - m*vIn[time][k - 1];//dk 
				}

				for (int s = (*this).getSizeGrid() - 1; s > 0; s--)
				{
					vNUM[time][s] = (vIn[0][s] / matrix[s][s]);
				}	
			}
		}
		else throw BadInput(" in calculNumerical() of Implicit upwind, see value of quickChoice"); //Exception BadInput

		//Write in the file the last line of the numerical solution
		Inputs_Outputs::writeNumericalInCSV(vNUM, "IMPLICIT UPWIND", (*this).getStep(), nbIt-1, (*this).getSizeGrid());
}

void ImplicitUpWind::isStable() //Test if the scheme is stable 
{
	std::cout << std::endl << "This system is always stable" << std::endl; //Inconditionnaly stable
}

//LAX-WENDROFF SCHEME

//CONSTRUCTOR

Lax_Wendroff::Lax_Wendroff(int choice, double step, int size, double cfl1) : Schemes(choice, step, size,cfl1) //The user choose f0 or f1, the time to run and the size of the grid
{
	std::cout << "The LAX WENDROFF scheme : forward in time and central in space" << std::endl; //Present the scheme 
}

//FUNCTIONS

void Lax_Wendroff::calculNumerical() //Compute the NUMERICAL solution of the LAX-WENDROFF scheme
{

	(*this).initNumerical(vNUM); //Initialize the vNUM vector of vector
	
		std::cout << "VALUES of f" << choice << " (Numerical solution) are calculated..." << std::endl;
		
		for (int t=0; t<nbIt-1; t++)
		{
			for (int space = 1; space < (*this).getSizeGrid() - 1; space++)
			{
				//Construction of vNUM following the scheme
				vNUM[t+1][space] = vNUM[t][space] - 0.5*cfl*(vNUM[t][space + 1] - vNUM[t][space - 1]) + 0.5*cfl*cfl*(vNUM[t][space + 1] - 2 * vNUM[t][space] + vNUM[t][space - 1]);
			}
		}

		//Write in the file the last line of the numerical solution
		Inputs_Outputs::writeNumericalInCSV(vNUM, "LAX-WENDROFF", (*this).getStep(), nbIt-1, (*this).getSizeGrid()); 
	}

void Lax_Wendroff::isStable() //Test if the scheme is stable 
{
	if (cfl<=1) // Stability condition
	{
		std::cout << std::endl << "This system is stable" << std::endl;
	}
	else std::cout << std::endl << "This system is not stable" << std::endl;
}

//RICHTMYER MULTI-STEP SCHEME

//CONSTRUCTOR

RichtmyerMS::RichtmyerMS(int choice, double step, int size, double cfl1) : Schemes(choice, step, size,cfl1)//The user choose f0 or f1, the time to run and the size of the grid
{
	std::cout << "The RICHTMYER MULTI-STEP scheme : 2 steps to construct" << std::endl; //Present the scheme 
}

//FUNCTIONS

void RichtmyerMS::calculNumerical() //Compute the NUMERICAL solution of the RICHTMYER MULTI-STEP scheme
{
	std::vector<std::vector<double>> vNumDemi;

	for (int i = 0; i < nbIt; i++)
	{
		vNumDemi.push_back(std::vector<double>((*this).getSizeGrid()));//create exactly the good size for vNumDemi
	}

	(*this).initNumerical(vNUM); //Initialize the vNUM vector of vector
	(*this).initNumerical(vNumDemi);
		std::cout << "VALUES of f" << choice << " (Numerical solution) are calculated..."<<std::endl;
		
		for (int t = 0; t<nbIt-1; t++)
		{
			for (int space = 1; space < (*this).getSizeGrid() - 1; space++) //Two steps to construct vNUM
			{
				//1st step
				vNumDemi[t + 1][space] = 0.5*(vNumDemi[t][space + 1] + vNumDemi[t][space - 1]) - (cfl / 4)*(vNumDemi[t][space + 1] - vNumDemi[t][space - 1]);
			}

			for (int space = 1; space < (*this).getSizeGrid() - 1; space++) 
			{
				//2nd step
				vNUM[t + 1][space] = vNumDemi[t][space] - (cfl / 2)*(vNumDemi[t + 1][space + 1] - vNumDemi[t + 1][space - 1]);
			}
		}

		//Write in the file the last line of the numerical solution
		Inputs_Outputs::writeNumericalInCSV(vNUM, "RICHTMYER", (*this).getStep(), nbIt-1, (*this).getSizeGrid());
}

void RichtmyerMS::isStable() //Test if the scheme is stable 
{
	if (cfl <= 2) // Stability condition
	{
		std::cout << "The scheme is stable" << std::endl;
	}
	else std::cout << "The scheme is not stable" << std::endl;
}