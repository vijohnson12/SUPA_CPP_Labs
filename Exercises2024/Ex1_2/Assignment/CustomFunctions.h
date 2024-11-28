#ifndef CUSTOMFUNCTIONS_H
#define CUSTOMFUNCTIONS_H

#include <vector>
#include <string>
#include <utility>

//functions
void readData(const std::string& filename, std::vector<std::pair<float, float>>& data);
void readErrors(const std::string& filename, std::vector<float>& errors);
std::vector<float> CalculateMagnitude(const std::vector<std::pair<float, float>>& data);
//void BestFit(const std::vector<std::pair<float, float>>& data, const std::string& outputFilename);
void BestFit(const std::vector<std::pair<float, float>>& data, const std::vector<float>&errors, const std::string& outputFile);
float CalculateXpowY(float x, float y);
std::vector<float> CalculateXpowYData(const std::vector<std::pair<float, float>>&data);

//printing
void print(const std::vector<std::pair<float, float>>& data, int n);
void print(const std::vector<float>& Magnitudes);
void print(const std::string& lineEquation, float chi2, int NDF, float chi2_NDF);
void print(const std::vector<std::pair<float, float>>& data, const std::vector<float>& results);

#endif 