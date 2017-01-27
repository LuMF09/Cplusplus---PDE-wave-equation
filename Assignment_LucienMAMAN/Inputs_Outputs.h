//Inputs_Outputs.h
#ifndef HEADER_INPUTS_OUTPUTS
#define HEADER_INPUTS_OUTPUTS
/*!
* \file Inputs_Outputs.h
* \brief Set of functions which manage inputs and outputs of the software
* \author Lucien MAMAN
* \date 05/12/2016
*/

#include <vector>

	class Inputs_Outputs
	{
	public: 

		/*!
		*\brief Ask and verify inputs choose by the user in order to avoid invalid arguent errors
		* This static function does not belong to Schemes class because we need to call it without create a scheme object
		*
		*\return A vector with all inputs verified in order to use all functions without invalid arguments.
		*/
		static std::vector<double> askImputs();

		/*!
		*\brief Write in a CSV file the vector with its parameters
		* This static function does not belong to Schemes class because we need to call it without create a scheme object
		*		
		*  \param v: vector to write in the CSV file
		*  \param nbPoint: Size of the grid
		*  \param step: the time of the simulation
		*  \param a: Kind of vector written in the CSV file
		*/
		static void writeVectorInCSV(std::vector<double> &v, int nbPoint, double step, std::string a);

		/*!
		*\brief Write the numerical solution, a vector of vector, in a CSV file
		* This static function does not belong to Schemes class because we need to call it without create a scheme object
		*
		*  \param v: vector of vector to write in the CSV file
		*  \param a: Let the user know wich scheme is used in the simulation
		*  \param step: the time of the simulation
		*  \param t: the time corresponding to the line of the matrix to write in the CSV
		*  \param size: the size of the grid
		*/
		static void writeNumericalInCSV(std::vector<std::vector<double>> &v, std::string a, double step, int t, int size);

		/*!
		*\brief Write norm infiny and norm 2 in a CSV file
		* This static function does not belong to Schemes class because we need to call it without create a scheme object
		*
		*  \param step: time of simulation
		*  \param norm0: Write the norm infiny in the CSV file
		*  \param norm2Normalized: Write the norm 2 normalized in the CSV file
		*/
		static void writeNormInCSV(double step, double norm0, double norm2Normalized);
	};

#endif