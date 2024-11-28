//I'm guessing you'll read this first, I used a lot of help on this because I thought it was very difficult. I spent probably 30+ hours on it and that's with help (mostly chat gpt). So I'm sorry I know you said not trust it but I did the best I could, for some reason c++ just isn't intuitive to me, although going through it I do feel more comfortable. I just learn better based of exmaples first then making my own code. Thanks for you patience in advance. 


#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <string>
#include <sstream>
#include <utility>
#include <cmath>
#include "CustomFunctions.h"
#include <limits>


using namespace std;

int main() {
   //filename
   std::string inputfile = "input2D_float.txt";
   std::string errorfile = "error2D_float.txt";
   std::string outputfile = "fitted_line.txt";
   std::string powerfile = "x^y.txt";
   std::string magnitudefile = "magnitudes.txt";

   //pairs
   std::vector<std::pair<float, float>> data;
   std::vector<float> errors;

   //read data
   readData(inputfile, data);
   readErrors(errorfile, errors);

   if(data.size() != errors.size()){
      std::cerr << "Mismatch between data and error points" << std::endl;
      return 1;
   }

   //x^y
   std::vector<float> results = CalculateXpowYData(data);

   //printing n lines
   //int n = 25;
   //printData(data, n);

   //calculating magnitudes
   //std::vector<float> Magnitudes = CalculateMagnitude(data);

   //std::cout<< "Magnitudes of the data points are:" << std::endl;
   //for(size_t i =0; i< Magnitudes.size(); ++i){
      //std::cout << "Magnitude of point " << i + 1 << " is: " << Magnitudes[i] <<std::endl;
   //}

int choice;
bool ContinueRunning = true;

while(ContinueRunning){
   std::cout << "Please select an operation" << std::endl;
   std::cout << "1. Print data" << std::endl;
   std::cout << "2. Calculate Magnitudes" << std::endl;
   std::cout << "3. Calculate line of best fit and Chi squared" << std::endl;
   std::cout << "4. Calculate x^y" << std::endl;
   std::cout << "5. Exit" << std::endl;

   std::cin >> choice;


   switch(choice){
      case 1: {
         int n;
         std::cout << "How many lines would you like to print?";
         std::cin >> n;
         print(data, n);
         break;
      }
      case 2:{
         vector<float> Magnitudes = CalculateMagnitude(data);
         print(Magnitudes);
         break;

      }
      case 3:{
         string outputFilename = "fitted_line.txt";
         BestFit(data, errors, outputfile);
         break;
      
      }  
      case 4:{
         vector<float> CalculateXpowY = CalculateXpowYData(data);
         print(data, results);
         break;
      }  
      case 5:{
         std::cout << "Exiting" << std::endl;
         ContinueRunning = false;
         break;

      }
   }

   if(ContinueRunning){
      char ContinueChoice;
      std::cout << "Do you want to perform another operation? Press y or n";
      std::cin >> ContinueChoice;

      if(ContinueChoice != 'y' && ContinueChoice != 'Y'){
         ContinueRunning = false;
         std::cout << "Exiting" << std::endl;
      }
   }  
}

   return 0;
}
