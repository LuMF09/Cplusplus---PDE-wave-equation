//Schemes.h
#ifndef HEADER_SCHEMES
#define HEADER_SCHEMES
/*!
* \file Schemes.h
* \brief Structure of the different schemes to study a linear advection equation
* \author Lucien MAMAN
* \date 05/12/2016
*/

#include<vector>

/*! \class Schemes
* \brief Abstract class Schemes
*
*  Create and define different functions and attribute which will be use by the other inherit classes
*/
class Schemes
{
	protected:

		std::vector<double>x; /*!< Vector used to implement the different space positions  */
		std::vector<double> vANAL; /*!< Vector used to implement the analytical solution */
		std::vector<std::vector <double>> vNUM; /*!< Vector of vector used to implement the numerical solution */

		int choice; /*!< Set of equations choose by the user */
		int nbIt; /*!< Number of iterations to go to the correct time step */

		double step;/*!< Time of simulation choose by the user */
		double deltaT;/*!< Delta T (time), equals to (time of simulation)  / (grid Size - 1) */
		double deltaX;/*!< Delta X (space), equals to abs(Left boudary - right boundary)/ (grid Size - 1) */
		double cfl; /*!< Courant Friedrich Number, usefull to study the stability, equals to U*delta T/delta X */

	public:

		/*!
		*  \brief Constructor of the abstract class Schemes
		*
		*  \param choice: the choosen set of equations
		*  \param step: time to run
		*  \param size: number of discretisation points
		*  \param cfl: the choosen CFL number
		*/
		Schemes(int choice,double step, int size, double cfl);

		//ACCESSORS

		/*!
		*  \brief Const Accessor method wich provide the choice made by the user
		*
		*  \return The choice made by the user
		*/
		double getChoice() const;

		/*!
		*  \brief Const Accessor method wich provide the size of the grid chose by the user
		*
		*  \return The size of the grid
		*/
		int getSizeGrid() const;

		/*!
		*  \brief Const Accessor method wich provide the time step chose by the user
		*
		*  \return The time step made by the user
		*/
		double getStep() const;

		//FUNCTIONS

		/*!
		*  \brief Const method used to calculate the analytical solution
		*/
		void calculAnalytical();

		/*!
		*  \brief Method used to calculate norms 1 and 2
		*/
		void calculNorms();

		/*!
		*  \brief Method used to intialize a vector of vector with the initial conditions
		*
		*  \param &v : a copy of a vector of vector
		*/
		void initNumerical(std::vector<std::vector <double>>&v);

		/*!
		*  \brief Method called in the main function to show the analytical solution, the numerical solution, norms and to know is the scheme is stable or not
		*/
		void showScheme();

		//PURE VIRTUAL FUNCTIONS

		/*!
		*  \brief Pure virtual function used to construct the numerical solutions
		*/
		virtual void calculNumerical() = 0;	

		/*!
		*  \brief Pure virtual function Test if a scheme is stable or not
		*
		*/
		virtual void isStable() = 0;
};

/*! \class ExplicitUpWind
* \brief Class inherit from Schemes class
*
*  Used to make numerical simulation of the explicit upwind scheme: FTBS
*/
class ExplicitUpWind : public Schemes
{

public:

	/*!
	*  \brief Constructor of the ExplicitUpWind scheme
	*
	*  \param choice: the choosen set of equations
	*  \param step: time to run
	*  \param size: number of discretisation points
	*  \param cfl: the choosen CFL number
	*/
	ExplicitUpWind(int choice, double step, int size, double cfl);

	/*!
	*  \brief Definition of the pure virtual function: construct the numerical solution of this scheme
	*/
	void calculNumerical();

	/*!
	*  \brief Definition of the pure virtual function: Test if a scheme is stable or not
	*
	*/
	void isStable();
};

/*! \class ImplicitUpWind
* \brief Class inherit from Schemes class
*
*  Used to make numerical simulation of the implicite upwind scheme: FTBS
*/
class ImplicitUpWind : public Schemes
{
	private:
	int quickChoice;

	public:

		/*!
		*  \brief Constructor of the ImplicitUpWind scheme
		*
		*  \param choice: the choosen set of equations
		*  \quickChoice: Choice of the way to compute the explicit upwind scheme
		*  \param step: time to run
		*  \param size: number of discretisation points
		*  \param cfl: the choosen CFL number
		*/
		ImplicitUpWind(int choice, int quickChoice, double step, int size, double cfl);

		/*!
		*  \brief Definition of the pure virtual function: construct the numerical solution of this scheme
		*/
		void calculNumerical();

		/*!
		*  \brief Definition of the pure virtual function: Test if a scheme is stable or not
		*
		*/
		void isStable();
};

/*! \class Lax_Wendroff
* \brief Class inherit from Schemes class
*
*  Used to make numerical simulation of the Lax-Wendroff scheme: FTCS
*/
class Lax_Wendroff : public Schemes
{

	public:
		/*!
		*  \brief Constructor of the Lax Wendroff scheme
		*
		*  \param choice: the choosen set of equations
		*  \param step: time to run
		*  \param size: number of discretisation points
		*  \param cfl: the choosen CFL number
		*/
		Lax_Wendroff(int choice, double step, int size, double cfl);

		/*!
		*  \brief Definition of the pure virtual function: construct the numerical solution of this scheme
		*/
		void calculNumerical();

		/*!
		*  \brief Definition of the pure virtual function: Test if a scheme is stable or not
		*
		*/
		void isStable();
};

/*! \class RichtmyerMS
* \brief Class inherit from Schemes class
*
*  Used to make numerical simulation of the Richtmyer multi-step scheme in two steps
*/
class RichtmyerMS : public Schemes
{
//private:
//	std::vector<std::vector<double>> vNumDemi;
	public:
		/*!
		*  \brief Constructor of the Richtmyer multi-step scheme
		*
		*  \param choice: the choosen set of equations
		*  \param step: time to run
		*  \param size: number of discretisation points
		*  \param cfl: the choosen CFL number
		*/
		RichtmyerMS(int choice, double step, int size, double cfl);

		/*!
		*  \brief Definition of the pure virtual function: construct the numerical solution of this scheme
		*/
		void calculNumerical();

		/*!
		*  \brief Definition of the pure virtual function: Test if a scheme is stable or not
		*
		*/
		void isStable();
};
#endif