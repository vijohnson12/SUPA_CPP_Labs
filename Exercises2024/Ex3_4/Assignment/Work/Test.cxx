#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <limits>
#include <sstream>
#include <cstdlib>
#include "FiniteFunctions.h"
#include <filesystem> //To check extensions in a nice way
#include "gnuplot-iostream.h" //Needed to produce plots (not part of the course) 





int main(){



   //reading data
   //std::ifstream inputfile("MysteryData21112.txt");

   //if(!inputfile){
      //std::cerr << "Error opening file" << std::endl;
      //return 1;
   //}

   //reading data
   //std::string line;
   //while(std::getline(inputfile, line)){
      //std::cout << line << std::endl;
   //}



   // if(!filename.is_open()){
   //    std::cerr << "Error opening file" << std::endl;
   //    return 1;
   // }

   // double point;
   // while(data_file >> point){
   //    data.push_back(point);
   // }

   //data_file.close();



   //inputfile.close();

   FiniteFunction finiteFunc(-5.0, 5.0, "outputs/png/mystery_output");

   std::vector<double> dataPoints;
   std::ifstream dataFile("MysteryData21112.txt");

   if(!dataFile.is_open()){
      std::cerr << "Could not open file" << std::endl;
      return 1;
   }

   double value;
   while(dataFile >> value){
      dataPoints.push_back(value);
   }

   //dataFile.close();

   finiteFunc.plotFunction();
   int Nbins = 18;
   finiteFunc.plotData(dataPoints, Nbins, true);
   finiteFunc.printInfo();

   //Final part with everything on one graph
   //dist
   double mu = -1, sigma =1.5;
   double range_min = -10.0, range_max = 10.0;
   NormalDistribution normal(mu, sigma, range_min, range_max, "NormalDist");

   //samples
   int n_samples = 1000;
   double step_size = 0.5;
   std::vector<double> sampled_data = normal.metropolisSample(n_samples, step_size);

   //std::vector<double> mystery_data;
   for(int i =0; i < 500; ++i){
      dataPoints.push_back(range_min + static_cast<double>(rand()) / RAND_MAX * (range_max - range_min));
   }
   
   normal.plotFunction();
   normal.plotData(sampled_data, 18, false);
   normal.plotData(dataPoints, 18, true);

   // //Part 2.1
   // double mu = -1;
   // double sigma = 1.5;
   // double range_min = -10.0;
   // double range_max = 10.0;
   // std::string outfile = "NormalDist";

   // NormalDistribution normalDist(mu, sigma, range_min, range_max, outfile);
   
   // dataPoints = normalDist.metropolisSample(10000, 0.5);

   // for(int i =0; i < 10; ++i){
   //    std::cout << dataPoints[i] << std::endl;
   // }

   // normalDist.plotFunction();
   // normalDist.plotData(dataPoints, 18, true);

   // //Original graph
   // //Pretty sure it's this one
   // NormalDistribution normalDist(-1,1.5, -10.0, 10.0, "NormalDist");
   // normalDist.integral(1000);
   // normalDist.printInfo();
   // normalDist.plotFunction();
   // normalDist.plotData(dataPoints, 18, true);


   //Zero chance it's this one unless my formula was written that badly
   //Yeah, it was my code
   //Doesn't fit as well as normal
   CauchyLorentzDistribution cauchyDist(-1, 1.35, -10.0, 10.0, "CauchyLorentzDist");
   cauchyDist.integral(1000);
   cauchyDist.printInfo();
   cauchyDist.plotFunction();
   cauchyDist.plotData(dataPoints, 18, true);

   //isn't this one
   NegativeCrystalBallDistribution crystalBallDist(-1, 2.25, 1.0, 0.5, -10.0, 10.0, "NegativeCrystalBallDist");
   crystalBallDist.integral(1000);
   crystalBallDist.printInfo();
   crystalBallDist.plotFunction();
   crystalBallDist.plotData(dataPoints, 18, true);

   dataFile.close();

   
   return 0;
   
}