#include "CustomFunctions.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

//function to read data
void readData(const std::string& filename, std::vector<std::pair<float, float>>& data){
   std::ifstream data_file(filename);


 //Check for opening file correctly
 if(!data_file.is_open()) {
    std::cout << "Error opening file: " << filename << std::endl;
    return;
 }
 else{
    std::cout << "File: " << filename<< "opened successfully!" << std::endl;

 }

 std:: string line;
//creating pairs
 while(std::getline(data_file, line)){
    std::stringstream ss(line);
    std::string x_str, y_str;

      if(std::getline(ss, x_str, ',') && std::getline(ss, y_str, ',')){
         try{
            float x = std::stof(x_str);
            float y = std::stof(y_str);
            data.push_back(std::make_pair(x,y));
         }catch(const std::invalid_argument& e){
            std::cerr << "Skipping invalid data:" << line << std::endl;
            continue;
         }
      }
 }

data_file.close();
}

//reading errors
void readErrors(const std::string& filename, std::vector<float>& errors){
   std::ifstream error_file(filename);

   if(!error_file.is_open()){
      std::cout << "Error opening file: " << filename << std::endl;
      return;
   }

   std::string line;
   while(std::getline(error_file, line)){
      std::stringstream ss(line);
      std::string error_str;
      if(std::getline(ss, error_str, ',')){
         try{
            float error = std::stof(error_str);
            errors.push_back(error);
         }catch(const std::invalid_argument& e){
            std::cerr << "Skipping invalid error data: " << line << std::endl;
            continue;
         }
      }
   }

   error_file.close();
}



//calculating magntiude function
//function name
std::vector<float> CalculateMagnitude(const std::vector<std::pair<float, float>>& data){
   //entire group of magnitudes
   std::vector<float> Magnitudes;

   for(const auto& pair: data){
      //calculation of single magnitude
      float Magnitude = std::sqrt((pair.first*pair.first)+(pair.second * pair.second));
      //appending to make group of magnitudes 
      Magnitudes.push_back(Magnitude);
   }
   return Magnitudes;
}

//function to print
void printData(const std::vector<std::pair<float, float>>&data,int n){
   if (n> data.size()){
      std::cout << "Warning: Requested" << n << "lines, but only " << data.size() << "data points available." << std::endl;
      n = 5;
   }

   for (int i = 0; i < n; ++i){
      std::cout << "x: " << data[i].first << ", y:" << data[i].second << std::endl;
   }
}

//best fit function
//function name
void BestFit(const std::vector<std::pair<float, float>>& data, const std::vector<float>& errors, const std::string& outputFilename){
if(data.size() != errors.size()){
   std::cerr << "Data points and error points are not of equal size" << std::endl;
   return;
}

//new variables
//reference assignment sheet with least squares method for following
float sumX = 0, sumY = 0, sumXY = 0, sumXX =0;
int N = data.size();

//Sums from least squares

for(const auto& pair: data){
   sumX += pair.first;
   sumY += pair.second;
   sumXY += pair.first * pair.second;
   sumXX += pair.first * pair.first;
}

//y=mx + b, m=slope, b=intercept

float m = (N*sumXY - sumX *sumY)/ (N*sumXX - sumX * sumX);
float b = ((sumXX *sumY)- (sumXY * sumX)) / (N*sumXX - sumX * sumX);

//chi squared calculation look at assignment sheet for equation
float chi2 = 0;

for(int i = 0; i < N; ++i){
   float y_expected = m * data[i].first + b;
   float difference = data[i].second - y_expected;
   float error = errors[i];
   chi2 += (difference * difference) / (error * error);
}

//degrees of freedom (NDF)
int NDF = N-2; 

//chi squared NDF
float chi2_NDF = chi2 / NDF;

std::string lineEquation = "y = " + std::to_string(m) + "x + " + std::to_string(b);

print(lineEquation, chi2, NDF, chi2_NDF);

//saving to same fit file
std::ofstream outFile(outputFilename);
if(outFile.is_open()){
   outFile << lineEquation << std::endl;
   outFile << "Chi squared = " << chi2 << ", Degrees of Freedom(NDF) = " << NDF << std::endl;
   outFile.close();
   std::cout << "The equation for chi squared has been saved to: " << outputFilename << std::endl;
}

//printing equation
std::cout << "The fitted line is y = " << m << "x +" << b << std::endl;

//saving equation to file

//std::ofstream outFile(outputFilename);
//if(outFile.is_open()){
   //outFile << "The fitted line is y = " << m << "x +" << b<< std::endl;
   //outFile.close();
   //std::cout << "The equation has been saved to " << outputFilename << std::endl;
   //}
}

//creating printing overload for data points
void print(const std::vector<std::pair<float, float>>& data, int n){
   if (n > data.size()){
      std::cout << "Only " << data.size() << "lines of data are available" << std::endl;
      n = data.size();
   }


for (int i =0; i < n; ++i){
   std::cout << "x: " << data[i].first << ", y: " << data[i].second << std::endl;
   }
}

//overload for magnitudes
void print(const std::vector<float>& magnitudes){
   std::cout << "Magnitudes are:" << std::endl;

   std::ofstream outFile("magnitudes.txt");

   for(size_t i = 0; i < magnitudes.size(); ++i){
      outFile << "Magnitude of point " << i + 1 << " is " << magnitudes[i] << std::endl;
   }
   outFile.close();

   for(size_t i =0; i< magnitudes.size(); ++i){
      std::cout << "Magnitude of point " << i+1 << " is: " << magnitudes[i] << std::endl;
   }
}

//overload for chi
void print(const std::string& lineEquation, float chi2, int NDF, float chi2_NDF){
   std::cout << "The fitted line is:" << lineEquation << std::endl;
   std::cout << "Chi squared = " << chi2 << ", Degrees of Freedom(NDF) = " << NDF << std::endl;
   std::cout << "Chi squared per NDF = " << chi2_NDF << std::endl;
}

//need to round first
float CalculateXpowY(float x, float y){
   int roundedY = static_cast<int>(std::round(y));


//calculating x^y
   if(roundedY ==0){
    return 1.0f;
   }if(roundedY > 0){
      return x * CalculateXpowY(x, roundedY - 1);
   }else{
      return 1.0f / CalculateXpowY(x, -roundedY);
   }
}

std::vector<float>CalculateXpowYData(const std::vector<std::pair<float, float>>& data){
   std::vector<float> results;
   for(const auto& pair : data){
      float x = pair.first;
      float y = pair.second;
      results.push_back(CalculateXpowY(x, y));
   }
   return results;
}

void print(const std::vector<std::pair<float, float>>& data, const std::vector<float>&results){
   if(data.size() != results.size()){
      std::cerr << "Data and power results don't match" << std::endl;
      return;
   }

   std::ofstream outFile("x^y.txt");
   outFile << "x^y values";
   std::cout << "x^y calculations: " << std::endl;
   for(size_t i = 0; i < data.size(); ++i){
      float x = data[i].first;
      float y = data[i].second;
      float result = results[i];

      outFile << "For x = " << x << " and y = " << y << ", x^y = " << result << std::endl;
      std::cout << "For x = " << x << " and y = " << y << ", x^y = " << result << std::endl;
   }
   outFile.close();
}
