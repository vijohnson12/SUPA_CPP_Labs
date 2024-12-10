//Again, I am so sorry, but good news, this took less time than the first lab
//I left everything in, so uncomment and recomment to pull up just the normaldist graph

#include <iostream>
#include <random>
#include <string>
#include <vector>
#include "FiniteFunctions.h"
#include <filesystem> //To check extensions in a nice way

#include "gnuplot-iostream.h" //Needed to produce plots (not part of the course) 

using std::filesystem::path;

//Empty constructor
FiniteFunction::FiniteFunction(){
  m_RMin = -5.0;
  m_RMax = 5.0;
  this->checkPath("DefaultFunction");
  m_Integral = 0;
}

//initialised constructor
FiniteFunction::FiniteFunction(double range_min, double range_max, std::string outfile){
  m_RMin = -10;
  m_RMax = 10;
  m_Integral = 0;
  this->checkPath(outfile); //Use provided string to name output files
}

//Plots are called in the destructor
//SUPACPP note: They syntax of the plotting code is not part of the course
FiniteFunction::~FiniteFunction(){
  Gnuplot gp; //Set up gnuplot object
  this->generatePlot(gp); //Generate the plot and save it to a png using "outfile" for naming 
}

/*
###################
//Setters
###################
*/ 
void FiniteFunction::setRangeMin(double RMin) {m_RMin = RMin;};
void FiniteFunction::setRangeMax(double RMax) {m_RMax = RMax;};
void FiniteFunction::setOutfile(std::string Outfile) {this->checkPath(Outfile);};

/*
###################
//Getters
###################
*/ 
double FiniteFunction::rangeMin() {return m_RMin;};
double FiniteFunction::rangeMax() {return m_RMax;};

/*
###################
//Function eval
###################
*/ 
double FiniteFunction::invxsquared(double x) {return 1/(1+x*x);};
double FiniteFunction::callFunction(double x) {return this->invxsquared(x);}; //(overridable)

/*
###################
Integration by hand (output needed to normalise function when plotting)
###################
*/ 
double FiniteFunction::integrate(int Ndiv){ //private
  //ToDo write an integrator
  double step = (m_RMax - m_RMin) / Ndiv;
  double sum = 0.0;
  for(int i = 0; i < Ndiv; ++i){
    double x1 = m_RMin + i * step;
    double x2 = m_RMin + (i + 1) * step;
    sum += 0.5 * (this->callFunction(x1) + this->callFunction(x2)) * step;
  }  
  return sum;
}
double FiniteFunction::integral(int Ndiv) { //public
  if (Ndiv <= 0){
    std::cout << "Invalid number of divisions for integral, setting Ndiv to 1000" <<std::endl;
    Ndiv = 1000;
  }
  if (m_Integral == NULL || Ndiv != m_IntDiv){
    m_IntDiv = Ndiv;
    m_Integral = this->integrate(Ndiv);
    return m_Integral;
  }
  else return m_Integral; //Don't bother re-calculating integral if Ndiv is the same as the last call
}

/*
###################
//Helper functions 
###################
*/
// Generate paths from user defined stem
void FiniteFunction::checkPath(std::string outfile){
 std::filesystem::path fp = outfile;
 m_FunctionName = fp.stem(); 
 m_OutData = m_FunctionName+".data";
 m_OutPng = m_FunctionName+".png";

 std::cout << "m_FunctionName: " << m_FunctionName << std::endl;
}

//Print (overridable)
void FiniteFunction::printInfo(){
  std::cout << "rangeMin: " << m_RMin << std::endl;
  std::cout << "rangeMax: " << m_RMax << std::endl;
  std::cout << "integral: " << m_Integral << ", calculated using " << m_IntDiv << " divisions" << std::endl;
  std::cout << "function: " << m_FunctionName << std::endl;
}

//Normal
NormalDistribution::NormalDistribution(double mu, double sigma, double range_min, double range_max, std::string outfile)
  : FiniteFunction(range_min, range_max, outfile), m_mu(mu), m_sigma(sigma){}

double NormalDistribution::callFunction(double x){
  return (1.0 / (m_sigma * std::sqrt(2 * M_PI))) * std::exp(-0.5 * std::pow((x - m_mu) / m_sigma, 2));
}

//Metropolis for Normal
std::vector<double> NormalDistribution::metropolisSample(int n_samples, double std_dev){
  std::vector<double> samples;
  std::random_device rd;
  std::mt19937 gen(rd());

  std::uniform_real_distribution<> uniform_dist(m_RMin, m_RMax);

  std::normal_distribution<> normal_dist(0.0, std_dev);

  //starting point
  double x_current = uniform_dist(gen);
  samples.push_back(x_current);

  for (int i = 1; i < n_samples; ++i){
    double x_candidate = x_current + normal_dist(gen);

    double f_current = this->callFunction(x_current);
    double f_candidate = this->callFunction(x_candidate);
    double A = std::min(f_candidate / f_current, 1.0);

    std::uniform_real_distribution<> acceptance_dis(0.0, 1.0);
    double T = acceptance_dis(gen);

    if(T < A){
      x_current = x_candidate;
    }
    samples.push_back(x_current);

  }
  return samples;
}



//Cauchy
CauchyLorentzDistribution::CauchyLorentzDistribution(double x0, double gamma, double range_min, double range_max, std::string outfile)
  : FiniteFunction(range_min, range_max, outfile), m_x0(x0), m_gamma(gamma){}

double CauchyLorentzDistribution::callFunction(double x){
  return (1.0 / ((M_PI * m_gamma) * (1 + std::pow((x - m_x0) / m_gamma, 2))));
}

//Negative
NegativeCrystalBallDistribution::NegativeCrystalBallDistribution(double x_bar, double sigma, double alpha, double n, double range_min, double range_max, std::string outfile)
  : FiniteFunction(range_min, range_max, outfile), m_x_bar(x_bar), m_sigma(sigma), m_alpha(alpha), m_n(n){
  m_N = 1.0 / m_sigma * (this->C(m_alpha, m_n) + this->D(m_alpha));
  }

double NegativeCrystalBallDistribution::callFunction(double x){
  double z = std::pow(x - m_x_bar,2) / m_sigma;
  if (z > -m_alpha) {
    return m_N * (std::exp(-0.5/m_sigma * std::pow(z, 2)));
  }else {
    double A = this->A(m_alpha, m_n);
    double B = this ->B(m_alpha, m_n);
    return m_N * (A * std::pow((B - z),-m_n));
  }
}

double NegativeCrystalBallDistribution::A(double alpha, double n){
  return (std::pow(n / std::abs(alpha), n)) * std::exp((std::abs(alpha), 2) / -2);
}

double NegativeCrystalBallDistribution::B(double alpha, double n){
  return n / std::abs(alpha) - std::abs(alpha);
}

double NegativeCrystalBallDistribution::C(double alpha, double n){
  return (n / std::abs(alpha)) * (1 / n-1) * std::exp(-std::pow(std::abs(alpha), 2) /2);
}

double NegativeCrystalBallDistribution::D(double alpha){
  return std::sqrt(M_PI / 2) * ( 1 + std::erf(std::abs(alpha)/ std::sqrt(2)));
}


/*
###################
//Plotting
###################
*/

//Hack because gnuplot-io can't read in custom functions, just scan over function and connect points with a line... 
void FiniteFunction::plotFunction(){
  m_function_scan = this->scanFunction(10000);
  m_plotfunction = true;
}

//Transform data points into a format gnuplot can use (histogram) and set flag to enable drawing of data to output plot
//set isdata to true (default) to plot data points in black, set to false to plot sample points in blue
void FiniteFunction::plotData(std::vector<double> &points, int Nbins, bool isdata){
  if (isdata){
    m_data = this->makeHist(points,Nbins);
    m_plotdatapoints = true;
  }
  else{
    m_samples = this->makeHist(points,Nbins);
    m_plotsamplepoints = true;
  }
}


/*
  #######################################################################################################
  ## SUPACPP Note:
  ## The three helper functions below are needed to get the correct format for plotting with gnuplot
  ## In theory you shouldn't have to touch them
  ## However it might be helpful to read through them and understand what they are doing
  #######################################################################################################
 */

//Scan over range of function using range/Nscan steps (just a hack so we can plot the function)
std::vector< std::pair<double,double> > FiniteFunction::scanFunction(int Nscan){
  std::vector< std::pair<double,double> > function_scan;
  double step = (m_RMax - m_RMin)/(double)Nscan;
  double x = m_RMin;
  //We use the integral to normalise the function points
  if (m_Integral == NULL) {
    std::cout << "Integral not set, doing it now" << std::endl;
    this->integral(Nscan);
    std::cout << "integral: " << m_Integral << ", calculated using " << Nscan << " divisions" << std::endl;
  }
  //For each scan point push back the x and y values 
  for (int i = 0; i < Nscan; i++){
    function_scan.push_back( std::make_pair(x,this->callFunction(x)/m_Integral));
    x += step;
  }
  return function_scan;
}

//Function to make histogram out of sampled x-values - use for input data and sampling
std::vector< std::pair<double,double> > FiniteFunction::makeHist(std::vector<double> &points, int Nbins){

  std::vector< std::pair<double,double> > histdata; //Plottable output shape: (midpoint,frequency)
  std::vector<int> bins(Nbins,0); //vector of Nbins ints with default value 0 
  int norm = 0;
  for (double point : points){
    //Get bin index (starting from 0) the point falls into using point value, range, and Nbins
    int bindex = static_cast<int>(floor((point-m_RMin)/((m_RMax-m_RMin)/(double)Nbins)));
    if (bindex<0 || bindex>Nbins){
      continue;
    }
    bins[bindex]++; //weight of 1 for each data point
    norm++; //Total number of data points
  }
  double binwidth = (m_RMax-m_RMin)/(double)Nbins;
  for (int i=0; i<Nbins; i++){
    double midpoint = m_RMin + i*binwidth + binwidth/2; //Just put markers at the midpoint rather than drawing bars
    double normdata = bins[i]/((double)norm*binwidth); //Normalise with N = 1/(Ndata*binwidth)
    histdata.push_back(std::make_pair(midpoint,normdata));
  }
  return histdata;
}

//Function which handles generating the gnuplot output, called in destructor
//If an m_plot... flag is set, the we must have filled the related data vector
//SUPACPP note: They syntax of the plotting code is not part of the course
void FiniteFunction::generatePlot(Gnuplot &gp){

  if (m_plotfunction==true && m_plotdatapoints==true && m_plotsamplepoints==true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "set style line 1 lt 1 lw 2 pi 1 ps 0\n";
    gp << "plot '-' with linespoints ls 1 title '"<<m_FunctionName<<"', '-' with points ps 2 lc rgb 'blue' title 'sampled data', '-' with points ps 1 lc rgb 'black' pt 7 title 'data'\n";
    gp.send1d(m_function_scan);
    gp.send1d(m_samples);
    gp.send1d(m_data);
  }
  else if (m_plotfunction==true && m_plotdatapoints==true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "set style line 1 lt 1 lw 2 pi 1 ps 0\n";
    gp << "plot '-' with linespoints ls 1 title '"<<m_FunctionName<<"', '-' with points ps 1 lc rgb 'black' pt 7 title 'data'\n";
    gp.send1d(m_function_scan);
    gp.send1d(m_data);
  }
  else if (m_plotfunction==true && m_plotsamplepoints==true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "set style line 1 lt 1 lw 2 pi 1 ps 0\n";
    gp << "plot '-' with linespoints ls 1 title '"<<m_FunctionName<<"', '-' with points ps 2 lc rgb 'blue' title 'sampled data'\n";
    gp.send1d(m_function_scan);
    gp.send1d(m_samples);
  }
  else if (m_plotfunction==true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "set style line 1 lt 1 lw 2 pi 1 ps 0\n";
    gp << "plot '-' with linespoints ls 1 title 'function'\n";
    gp.send1d(m_function_scan);
  }

  else if (m_plotdatapoints == true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "plot '-' with points ps 1 lc rgb 'black' pt 7 title 'data'\n";
    gp.send1d(m_data);
  }

  else if (m_plotsamplepoints == true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "plot '-' with points ps 2 lc rgb 'blue' title 'sampled data'\n";
    gp.send1d(m_samples);
  }
}
