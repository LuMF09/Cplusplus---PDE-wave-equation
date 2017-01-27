#include "Errors.h"
#include <string>
#include <iostream>

//DivideZero to signal a division by 0

//CONSTRUCTOR
DivideZero::DivideZero(std::string a) : s("ERROR! Division by 0" + a){}

//FUNCTION WHAT()
char const* DivideZero::what() const // display a personalized message
{
	return s.c_str();
}

//TooLittleTab to signal if the table used to store inputs is bad implemented

//CONSTRUCTOR
TooLittleTab::TooLittleTab(std::string a) : s("ERROR! The tab created is too little" + a){}

//FUNCTION WHAT()
char const* TooLittleTab::what() const// display a personalized message
{
	return s.c_str();
}

//DocOpen to signal if the document in which we write results or norms is already open

//CONSTRUCTOR
DocOpen::DocOpen(std::string a) : s("ERROR! The document " + a + " is open. Please close it"){}

//FUNCTION WHAT()
char const* DocOpen::what() const // display a personalized message
{
	return s.c_str();
}

//BadInput to signal if the input expected does not exist

//CONSTRUCTOR
BadInput::BadInput(std::string a) : s("ERROR! None of the if condition is checked "+ a){}

//FUNCTION WHAT()
char const* BadInput::what() const // display a personalized message
{
	return s.c_str();
}