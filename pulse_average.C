// This program takes the Y projection of
// h_daq and saves it in a separate TH1D histogram.
// SUch histogram has multiple peaks. To begin with the fitting
// process, here we try to average these waveforms.
// Once WF are averaged, we fit them with the newest
// electronics response function.

// IN PROGRESS: averaging waveforms!

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TH2.h>
#include <TH1.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TPad.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <TCanvas.h>
#include <TPad.h>
#include "TGraph.h"
#include <TGraphErrors.h>
#include <TSystem.h>
#include <math.h>
#include <TROOT.h>
#include <TAttLine.h>
#include <TLegend.h>
#include <TString.h>
#include <cmath> // For std::abs
#include "TFitResult.h"
#include <fstream>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <limits>

#include <fftw3.h>

#include "TPaveStats.h"
#include "TPaveLabel.h"

using std::cout;
using std::cin;
using std::endl;

std::vector<int> positivePeakBins;
std::vector<double> positivePeakValues;
std::vector<int> negativePeakStartBins;
std::vector<int> negativePeakEndBins;
std::vector<double> negativePeakStartValues;
std::vector<double> negativePeakEndValues;


int run = 21040;

///////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////DEFINING FIT FUNCTIONS/////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

// New electronics response function:

double response(double *x, double *par){
  Double_t t = x[0]-par[0];
  //Double_t t = par[0];
  Double_t A0 = par[1];
  Double_t tp = par[2];
  Double_t CT = 1./1.996;
  Double_t A = A0 * 2.7433/pow(tp*CT,4);
  Double_t p0 = 1.477/tp/CT;
  Double_t pr1 = 1.417/tp/CT;
  Double_t pr2 = 1.204/tp/CT;
  Double_t pi1 = 0.598/tp/CT;
  Double_t pi2 = 1.299/tp/CT;


  double k3 = par[3];
  double k4 = par[4];
  double k5 = par[5];
  double k6 = par[6];

  double value = A*((-(k3*k4) + pow(k4,2) + k3*k5 - k4*k5)/(exp(k4*t)*(k4 - k6)*(k4 - p0)*(pow(k4,2) + pow(pi1,2) - 2*k4*pr1 + pow(pr1,2))*(pow(k4,2) + pow(pi2,2) - 2*k4*pr2 + pow(pr2,2))) +
	     (-(k3*k5) + k3*k6 + k5*k6 - pow(k6,2))/(exp(k6*t)*(k4 - k6)*(k6 - p0)*(pow(k6,2) + pow(pi1,2) - 2*k6*pr1 + pow(pr1,2))*(pow(k6,2) + pow(pi2,2) - 2*k6*pr2 + pow(pr2,2))) +
	     (-(k3*k5) + k3*p0 + k5*p0 - pow(p0,2))/(exp(p0*t)*(k4 - p0)*(-k6 + p0)*(pow(p0,2) + pow(pi1,2) - 2*p0*pr1 + pow(pr1,2))*(pow(p0,2) + pow(pi2,2) - 2*p0*pr2 + pow(pr2,2))) +
     (pi1*((pow(pi1,2) + pow(pr1,2))*(2*k6*(pow(pi1,2) + pow(pr1,2))*(pr1 - pr2) + k6*p0*(-pow(pi1,2) + pow(pi2,2) - pow(pr1,2) + pow(pr2,2)) + (pow(pi1,2) + pow(pr1,2))*(pow(pi1,2) - pow(pi2,2) + (pr1 - pr2)*(2*p0 - 3*pr1 + pr2)) +
              k5*(2*pow(pi1,2)*(-2*pr1 + pr2) + p0*(pow(pi1,2) - pow(pi2,2) - 3*pow(pr1,2) + 4*pr1*pr2 - pow(pr2,2)) + 2*pr1*(pow(pi2,2) + 2*pow(pr1,2) - 3*pr1*pr2 + pow(pr2,2)) +
                 k6*(pow(pi1,2) - pow(pi2,2) + (pr1 - pr2)*(2*p0 - 3*pr1 + pr2)))) + k4*((pow(pi1,2) + pow(pr1,2))*(2*(pow(pi1,2) + pow(pr1,2))*(pr1 - pr2) + p0*(-pow(pi1,2) + pow(pi2,2) - pow(pr1,2) + pow(pr2,2))) +
              k5*(2*k6*(pow(pi1,2) + pow(pr1,2))*(pr1 - pr2) - k6*p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1,2) - pow(pr2,2)) + (pow(pi1,2) + pow(pr1,2))*(pow(pi1,2) - pow(pi2,2) + (pr1 - pr2)*(2*p0 - 3*pr1 + pr2))) +
              k6*(-pow(pi1,4) + pow(pi1,2)*(pow(pi2,2) - 2*pow(pr1,2) + 2*p0*pr2 + pow(pr2,2)) + pr1*(pr1*(pow(pi2,2) - pow(pr1,2) + pow(pr2,2)) - 2*p0*(pow(pi2,2) - pr1*pr2 + pow(pr2,2))))) +
           k3*(-((pow(pi1,2) + pow(pr1,2))*(4*pow(pi1,2)*pr1 - 2*pow(pi2,2)*pr1 - 4*pow(pr1,3) - 2*pow(pi1,2)*pr2 + 6*pow(pr1,2)*pr2 - 2*pr1*pow(pr2,2) + p0*(-pow(pi1,2) + pow(pi2,2) + 3*pow(pr1,2) - 4*pr1*pr2 + pow(pr2,2)) +
                   k6*(-pow(pi1,2) + pow(pi2,2) - (pr1 - pr2)*(2*p0 - 3*pr1 + pr2)))) + k5*(-pow(pi1,4) + pow(pi1,2)*(pow(pi2,2) - 4*p0*pr1 + 10*pow(pr1,2) + 2*p0*pr2 - 8*pr1*pr2 + pow(pr2,2)) +
                 k6*(2*pow(pi1,2)*(-2*pr1 + pr2) + p0*(pow(pi1,2) - pow(pi2,2) - 3*pow(pr1,2) + 4*pr1*pr2 - pow(pr2,2)) + 2*pr1*(pow(pi2,2) + 2*pow(pr1,2) - 3*pr1*pr2 + pow(pr2,2))) +
                 pr1*(2*p0*(pow(pi2,2) + 2*pow(pr1,2) - 3*pr1*pr2 + pow(pr2,2)) - pr1*(3*pow(pi2,2) + 5*pow(pr1,2) - 8*pr1*pr2 + 3*pow(pr2,2)))) +
              k4*(2*k6*(pow(pi1,2) + pow(pr1,2))*(pr1 - pr2) + k6*p0*(-pow(pi1,2) + pow(pi2,2) - pow(pr1,2) + pow(pr2,2)) + (pow(pi1,2) + pow(pr1,2))*(pow(pi1,2) - pow(pi2,2) + (pr1 - pr2)*(2*p0 - 3*pr1 + pr2)) +
                 k5*(2*pow(pi1,2)*(-2*pr1 + pr2) + p0*(pow(pi1,2) - pow(pi2,2) - 3*pow(pr1,2) + 4*pr1*pr2 - pow(pr2,2)) + 2*pr1*(pow(pi2,2) + 2*pow(pr1,2) - 3*pr1*pr2 + pow(pr2,2)) +
                    k6*(pow(pi1,2) - pow(pi2,2) + (pr1 - pr2)*(2*p0 - 3*pr1 + pr2))))))*cos(pi1*t) -
        ((pow(pi1,2) + pow(pr1,2))*((pow(pi1,2) + pow(pr1,2))*(p0*(pow(pi1,2) - pow(pi2,2) - pow(pr1 - pr2,2)) + pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(-3*pr1 + 2*pr2)) +
              k6*(pow(pi1,4) + (p0 - pr1)*pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) - pow(pi1,2)*(pow(pi2,2) - p0*pr1 + 2*p0*pr2 - 2*pr1*pr2 + pow(pr2,2))) +
              k5*(-pow(pi1,4) + (p0 - pr1)*pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(pow(pi2,2) - 3*p0*pr1 + 6*pow(pr1,2) + 2*p0*pr2 - 6*pr1*pr2 + pow(pr2,2)) +
                 k6*(p0*(pow(pi1,2) - pow(pi2,2) - pow(pr1 - pr2,2)) + pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(-3*pr1 + 2*pr2)))) +
           k4*((pow(pi1,2) + pow(pr1,2))*(pow(pi1,4) + (p0 - pr1)*pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) - pow(pi1,2)*(pow(pi2,2) - p0*pr1 + 2*p0*pr2 - 2*pr1*pr2 + pow(pr2,2))) +
              k5*((pow(pi1,2) + pow(pr1,2))*(p0*(pow(pi1,2) - pow(pi2,2) - pow(pr1 - pr2,2)) + pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(-3*pr1 + 2*pr2)) +
                 k6*(pow(pi1,4) + (p0 - pr1)*pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) - pow(pi1,2)*(pow(pi2,2) - p0*pr1 + 2*p0*pr2 - 2*pr1*pr2 + pow(pr2,2)))) +
              k6*((pow(pi1,2) + pow(pr1,2))*(pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(pr1 - 2*pr2)) -
                 p0*(pow(pi1,4) + pow(pr1,2)*(pow(pi2,2) + pow(pr1 - pr2,2)) - pow(pi1,2)*(pow(pi2,2) - 2*pow(pr1,2) + 2*pr1*pr2 + pow(pr2,2))))) +
           k3*((pow(pi1,2) + pow(pr1,2))*(-pow(pi1,4) + (p0 - pr1)*pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(pow(pi2,2) - 3*p0*pr1 + 6*pow(pr1,2) + 2*p0*pr2 - 6*pr1*pr2 + pow(pr2,2)) +
                 k6*(p0*(pow(pi1,2) - pow(pi2,2) - pow(pr1 - pr2,2)) + pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(-3*pr1 + 2*pr2))) +
              k5*(5*pow(pi1,4)*pr1 - 3*pow(pi1,2)*pow(pi2,2)*pr1 - 10*pow(pi1,2)*pow(pr1,3) + pow(pi2,2)*pow(pr1,3) + pow(pr1,5) - 2*pow(pi1,4)*pr2 + 12*pow(pi1,2)*pow(pr1,2)*pr2 - 2*pow(pr1,4)*pr2 - 3*pow(pi1,2)*pr1*pow(pr2,2) +
                 pow(pr1,3)*pow(pr2,2) - p0*(pow(pi1,4) + pow(pr1,2)*(pow(pi2,2) + pow(pr1 - pr2,2)) - pow(pi1,2)*(pow(pi2,2) + 6*pow(pr1,2) - 6*pr1*pr2 + pow(pr2,2))) +
                 k6*(-pow(pi1,4) + (p0 - pr1)*pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(pow(pi2,2) - 3*p0*pr1 + 6*pow(pr1,2) + 2*p0*pr2 - 6*pr1*pr2 + pow(pr2,2)))) +
              k4*((pow(pi1,2) + pow(pr1,2))*(p0*(pow(pi1,2) - pow(pi2,2) - pow(pr1 - pr2,2)) + pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(-3*pr1 + 2*pr2)) +
                 k6*(pow(pi1,4) + (p0 - pr1)*pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) - pow(pi1,2)*(pow(pi2,2) - p0*pr1 + 2*p0*pr2 - 2*pr1*pr2 + pow(pr2,2))) +
                 k5*(-pow(pi1,4) + (p0 - pr1)*pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(pow(pi2,2) - 3*p0*pr1 + 6*pow(pr1,2) + 2*p0*pr2 - 6*pr1*pr2 + pow(pr2,2)) +
                    k6*(p0*(pow(pi1,2) - pow(pi2,2) - pow(pr1 - pr2,2)) + pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(-3*pr1 + 2*pr2))))))*sin(pi1*t))/
	     (exp(pr1*t)*pi1*(pow(k4,2) + pow(pi1,2) - 2*k4*pr1 + pow(pr1,2))*(pow(k6,2) + pow(pi1,2) - 2*k6*pr1 + pow(pr1,2))*(pow(p0,2) + pow(pi1,2) - 2*p0*pr1 + pow(pr1,2))*
        (pow(pi1,4) - 2*pow(pi1,2)*(pow(pi2,2) - pow(pr1 - pr2,2)) + pow(pow(pi2,2) + pow(pr1 - pr2,2),2))) +
     (-(pi2*(k4*(-((pow(pi2,2) + pow(pr2,2))*(p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1,2) - pow(pr2,2)) - 2*(pr1 - pr2)*(pow(pi2,2) + pow(pr2,2)))) +
                k5*((pow(pi1,2) - pow(pi2,2) + (2*p0 + pr1 - 3*pr2)*(pr1 - pr2))*(pow(pi2,2) + pow(pr2,2)) + 2*k6*(pr1 - pr2)*(pow(pi2,2) + pow(pr2,2)) + k6*p0*(-pow(pi1,2) + pow(pi2,2) - pow(pr1,2) + pow(pr2,2))) +
                k6*(pow(pi2,4) - pow(pi2,2)*(2*p0*pr1 + pow(pr1,2) - 2*pow(pr2,2)) - pow(pi1,2)*(pow(pi2,2) + pr2*(-2*p0 + pr2)) - (pr1 - pr2)*pr2*(-2*p0*pr1 + pr2*(pr1 + pr2)))) +
             (pow(pi2,2) + pow(pr2,2))*((pow(pi1,2) - pow(pi2,2) + (2*p0 + pr1 - 3*pr2)*(pr1 - pr2))*(pow(pi2,2) + pow(pr2,2)) + 2*k6*(pr1 - pr2)*(pow(pi2,2) + pow(pr2,2)) + k6*p0*(-pow(pi1,2) + pow(pi2,2) - pow(pr1,2) + pow(pr2,2)) +
                k5*(k6*(pow(pi1,2) - pow(pi2,2) + (2*p0 + pr1 - 3*pr2)*(pr1 - pr2)) + p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1,2) - 4*pr1*pr2 + 3*pow(pr2,2)) -
                   2*(pow(pi2,2)*(pr1 - 2*pr2) + pr2*(pow(pi1,2) + pow(pr1,2) - 3*pr1*pr2 + 2*pow(pr2,2))))) +
             k3*((pow(pi2,2) + pow(pr2,2))*(k6*(pow(pi1,2) - pow(pi2,2) + (2*p0 + pr1 - 3*pr2)*(pr1 - pr2)) + p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1,2) - 4*pr1*pr2 + 3*pow(pr2,2)) -
                   2*(pow(pi2,2)*(pr1 - 2*pr2) + pr2*(pow(pi1,2) + pow(pr1,2) - 3*pr1*pr2 + 2*pow(pr2,2)))) +
                k5*(pow(pi2,4) - 2*p0*pow(pi2,2)*pr1 - pow(pi2,2)*pow(pr1,2) + 4*p0*pow(pi2,2)*pr2 + 8*pow(pi2,2)*pr1*pr2 - 2*p0*pow(pr1,2)*pr2 - 10*pow(pi2,2)*pow(pr2,2) + 6*p0*pr1*pow(pr2,2) + 3*pow(pr1,2)*pow(pr2,2) -
                   4*p0*pow(pr2,3) - 8*pr1*pow(pr2,3) + 5*pow(pr2,4) - pow(pi1,2)*(pow(pi2,2) + (2*p0 - 3*pr2)*pr2) + k6*p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1,2) - 4*pr1*pr2 + 3*pow(pr2,2)) -
                   2*k6*(pow(pi2,2)*(pr1 - 2*pr2) + pr2*(pow(pi1,2) + pow(pr1,2) - 3*pr1*pr2 + 2*pow(pr2,2)))) +
                k4*((pow(pi1,2) - pow(pi2,2) + (2*p0 + pr1 - 3*pr2)*(pr1 - pr2))*(pow(pi2,2) + pow(pr2,2)) + 2*k6*(pr1 - pr2)*(pow(pi2,2) + pow(pr2,2)) + k6*p0*(-pow(pi1,2) + pow(pi2,2) - pow(pr1,2) + pow(pr2,2)) +
                   k5*(k6*(pow(pi1,2) - pow(pi2,2) + (2*p0 + pr1 - 3*pr2)*(pr1 - pr2)) + p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1,2) - 4*pr1*pr2 + 3*pow(pr2,2)) -
                      2*(pow(pi2,2)*(pr1 - 2*pr2) + pr2*(pow(pi1,2) + pow(pr1,2) - 3*pr1*pr2 + 2*pow(pr2,2)))))))*cos(pi2*t)) +
        ((pow(pi2,2) + pow(pr2,2))*((pow(pi2,2) + pow(pr2,2))*(p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1 - pr2,2)) - (pow(pi1,2) + pow(pr1 - pr2,2))*pr2 + pow(pi2,2)*(-2*pr1 + 3*pr2)) +
              k6*(-pow(pi2,4) - (p0 - pr2)*pow(pr1 - pr2,2)*pr2 + pow(pi2,2)*(2*p0*pr1 + pow(pr1,2) - p0*pr2 - 2*pr1*pr2) + pow(pi1,2)*(pow(pi2,2) + pr2*(-p0 + pr2))) +
              k5*(-(pow(pi1,2)*pow(pi2,2)) + pow(pi2,4) - 2*p0*pow(pi2,2)*pr1 - pow(pi2,2)*pow(pr1,2) - p0*pow(pi1,2)*pr2 + 3*p0*pow(pi2,2)*pr2 + 6*pow(pi2,2)*pr1*pr2 - p0*pow(pr1,2)*pr2 + pow(pi1,2)*pow(pr2,2) -
                 6*pow(pi2,2)*pow(pr2,2) + 2*p0*pr1*pow(pr2,2) + pow(pr1,2)*pow(pr2,2) - p0*pow(pr2,3) - 2*pr1*pow(pr2,3) + pow(pr2,4) +
                 k6*(p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1 - pr2,2)) - (pow(pi1,2) + pow(pr1 - pr2,2))*pr2 + pow(pi2,2)*(-2*pr1 + 3*pr2)))) +
           k4*(-(k6*(pow(pi2,2) + pow(pr2,2))*((pow(pi1,2) + pow(pr1 - pr2,2))*pr2 + pow(pi2,2)*(-2*pr1 + pr2))) +
              k6*p0*(pow(pi2,4) + pow(pr1 - pr2,2)*pow(pr2,2) - pow(pi2,2)*(pow(pr1,2) + 2*pr1*pr2 - 2*pow(pr2,2)) + pow(pi1,2)*(-pow(pi2,2) + pow(pr2,2))) +
              (pow(pi2,2) + pow(pr2,2))*(-pow(pi2,4) - (p0 - pr2)*pow(pr1 - pr2,2)*pr2 + pow(pi2,2)*(2*p0*pr1 + pow(pr1,2) - p0*pr2 - 2*pr1*pr2) + pow(pi1,2)*(pow(pi2,2) + pr2*(-p0 + pr2))) +
              k5*((pow(pi2,2) + pow(pr2,2))*(p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1 - pr2,2)) - (pow(pi1,2) + pow(pr1 - pr2,2))*pr2 + pow(pi2,2)*(-2*pr1 + 3*pr2)) +
                 k6*(-pow(pi2,4) - (p0 - pr2)*pow(pr1 - pr2,2)*pr2 + pow(pi2,2)*(2*p0*pr1 + pow(pr1,2) - p0*pr2 - 2*pr1*pr2) + pow(pi1,2)*(pow(pi2,2) + pr2*(-p0 + pr2))))) +
           k3*((pow(pi2,2) + pow(pr2,2))*(-(pow(pi1,2)*pow(pi2,2)) + pow(pi2,4) - 2*p0*pow(pi2,2)*pr1 - pow(pi2,2)*pow(pr1,2) - p0*pow(pi1,2)*pr2 + 3*p0*pow(pi2,2)*pr2 + 6*pow(pi2,2)*pr1*pr2 - p0*pow(pr1,2)*pr2 +
                 pow(pi1,2)*pow(pr2,2) - 6*pow(pi2,2)*pow(pr2,2) + 2*p0*pr1*pow(pr2,2) + pow(pr1,2)*pow(pr2,2) - p0*pow(pr2,3) - 2*pr1*pow(pr2,3) + pow(pr2,4) +
                 k6*(p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1 - pr2,2)) - (pow(pi1,2) + pow(pr1 - pr2,2))*pr2 + pow(pi2,2)*(-2*pr1 + 3*pr2))) -
              k5*(p0*pow(pi1,2)*pow(pi2,2) - p0*pow(pi2,4) - 2*pow(pi2,4)*pr1 + p0*pow(pi2,2)*pow(pr1,2) - 3*pow(pi1,2)*pow(pi2,2)*pr2 + 5*pow(pi2,4)*pr2 - 6*p0*pow(pi2,2)*pr1*pr2 - 3*pow(pi2,2)*pow(pr1,2)*pr2 -
                 p0*pow(pi1,2)*pow(pr2,2) + 6*p0*pow(pi2,2)*pow(pr2,2) + 12*pow(pi2,2)*pr1*pow(pr2,2) - p0*pow(pr1,2)*pow(pr2,2) + pow(pi1,2)*pow(pr2,3) - 10*pow(pi2,2)*pow(pr2,3) + 2*p0*pr1*pow(pr2,3) + pow(pr1,2)*pow(pr2,3) -
                 p0*pow(pr2,4) - 2*pr1*pow(pr2,4) + pow(pr2,5) + k6*(-pow(pi2,4) + (p0 - pr2)*pow(pr1 - pr2,2)*pr2 + pow(pi1,2)*(pow(pi2,2) + (p0 - pr2)*pr2) + pow(pi2,2)*(2*p0*pr1 + pow(pr1,2) - 3*p0*pr2 - 6*pr1*pr2 + 6*pow(pr2,2)))) +
              k4*((pow(pi2,2) + pow(pr2,2))*(p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1 - pr2,2)) - (pow(pi1,2) + pow(pr1 - pr2,2))*pr2 + pow(pi2,2)*(-2*pr1 + 3*pr2)) +
                 k6*(-pow(pi2,4) - (p0 - pr2)*pow(pr1 - pr2,2)*pr2 + pow(pi2,2)*(2*p0*pr1 + pow(pr1,2) - p0*pr2 - 2*pr1*pr2) + pow(pi1,2)*(pow(pi2,2) + pr2*(-p0 + pr2))) +
                 k5*(-(pow(pi1,2)*pow(pi2,2)) + pow(pi2,4) - 2*p0*pow(pi2,2)*pr1 - pow(pi2,2)*pow(pr1,2) - p0*pow(pi1,2)*pr2 + 3*p0*pow(pi2,2)*pr2 + 6*pow(pi2,2)*pr1*pr2 - p0*pow(pr1,2)*pr2 + pow(pi1,2)*pow(pr2,2) -
                    6*pow(pi2,2)*pow(pr2,2) + 2*p0*pr1*pow(pr2,2) + pow(pr1,2)*pow(pr2,2) - p0*pow(pr2,3) - 2*pr1*pow(pr2,3) + pow(pr2,4) +
                    k6*(p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1 - pr2,2)) - (pow(pi1,2) + pow(pr1 - pr2,2))*pr2 + pow(pi2,2)*(-2*pr1 + 3*pr2))))))*sin(pi2*t))/
	     (exp(pr2*t)*pi2*(pow(pi1,4) - 2*pow(pi1,2)*(pow(pi2,2) - pow(pr1 - pr2,2)) + pow(pow(pi2,2) + pow(pr1 - pr2,2),2))*(pow(k4,2) + pow(pi2,2) - 2*k4*pr2 + pow(pr2,2))*(pow(k6,2) + pow(pi2,2) - 2*k6*pr2 + pow(pr2,2))*
       (pow(p0,2) + pow(pi2,2) - 2*p0*pr2 + pow(pr2,2))));


  if (t > 0)
    return value;
  else
    return 0;
}


///////////////////////////////////////////////////////////////////////////////////////////

// Ideal electronics response function:

double response_legacy(double *x, double *par){
  Double_t t = x[0]-par[0];
  Double_t A0 = par[1];
  Double_t tp = par[2];

  Double_t reltime = t / tp;
  Double_t gain = A0* 1.012;

  Double_t value = 4.31054 * exp(-2.94809 * reltime) * gain - 2.6202 * exp(-2.82833 * reltime) * cos(1.19361 * reltime) * gain -
           2.6202 * exp(-2.82833 * reltime) * cos(1.19361 * reltime) * cos(2.38722 * reltime) * gain +
           0.464924 * exp(-2.40318 * reltime) * cos(2.5928 * reltime) * gain +
           0.464924 * exp(-2.40318 * reltime) * cos(2.5928 * reltime) * cos(5.18561 * reltime) * gain +
           0.762456 * exp(-2.82833 * reltime) * sin(1.19361 * reltime) * gain -
           0.762456 * exp(-2.82833 * reltime) * cos(2.38722 * reltime) * sin(1.19361 * reltime) * gain +
           0.762456 * exp(-2.82833 * reltime) * cos(1.19361 * reltime) * sin(2.38722 * reltime) * gain -
           2.620200 * exp(-2.82833 * reltime) * sin(1.19361 * reltime) * sin(2.38722 * reltime) * gain -
           0.327684 * exp(-2.40318 * reltime) * sin(2.5928 * reltime) * gain +
           +0.327684 * exp(-2.40318 * reltime) * cos(5.18561 * reltime) * sin(2.5928 * reltime) * gain -
           0.327684 * exp(-2.40318 * reltime) * cos(2.5928 * reltime) * sin(5.18561 * reltime) * gain +
           0.464924 * exp(-2.40318 * reltime) * sin(2.5928 * reltime) * sin(5.18561 * reltime) * gain;


  if (t > 0 and t<10)
    return value;
  else
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////


// To scale ticks to microseconds

Double_t ScaleX(Double_t x){
  Double_t v;
  v = (0.512 * x);
  return v;
}

void ScaleAxis(TAxis *a, Double_t (*Scale)(Double_t)){

    if (!a) return; // just a precaution
    if (a->GetXbins()->GetSize())
        {
            // an axis with variable bins
            // note: bins must remain in increasing order, hence the "Scale"
            // function must be strictly (monotonically) increasing
            TArrayD X(*(a->GetXbins()));
            for(Int_t i = 0; i < X.GetSize(); i++) X[i] = Scale(X[i]);
            a->Set((X.GetSize() - 1), X.GetArray()); // new Xbins
        }
    else
        {
            // an axis with fix bins
            // note: we modify Xmin and Xmax only, hence the "Scale" function
            // must be linear (and Xmax must be greater than Xmin)
            a->Set( a->GetNbins(),
                    Scale(a->GetXmin()), // new Xmin
                    Scale(a->GetXmax()) ); // new Xmax
        }
    return;

    }

void ScaleXaxis(TH1 *h, Double_t (*Scale)(Double_t)){

    if (!h) return; // just a precaution
    ScaleAxis(h->GetXaxis(), Scale);
    return;

    }


/////////////////////////////////////////////////////////////////////////////////////////

// To find the peaks within each raw waveform, for each channel:

void FindPeaks(TH1D* hist, double positiveThreshold, double negativeThreshold,
               std::vector<int>& positivePeakBins, std::vector<double>& positivePeakValues,
               std::vector<int>& negativePeakBins, std::vector<double>& negativePeakValues) {

    positivePeakBins.clear();
    positivePeakValues.clear();
    negativePeakBins.clear();
    negativePeakValues.clear();

    bool inPositivePeakRegion = false;
    bool foundPositivePeak = false;
    double maxNegativePeakValue = std::numeric_limits<double>::max();
    int maxNegativePeakBin = -1;
    bool inNegativePeak = false;


    for (int i = 1; i <= hist->GetNbinsX(); ++i) {

        double current = hist->GetBinContent(i);
        double next = i < hist->GetNbinsX() ? hist->GetBinContent(i + 1) : std::numeric_limits<double>::max();
        double prev = i > 1 ? hist->GetBinContent(i - 1) : std::numeric_limits<double>::max();


        // Positive peak detection:
        if (!inPositivePeakRegion && current > positiveThreshold && current > prev && current > next) {

            inPositivePeakRegion = true;
            foundPositivePeak = true;
            // Reset for next negative peak
            maxNegativePeakValue = std::numeric_limits<double>::max();
            maxNegativePeakBin = -1;
        } else if (inPositivePeakRegion) {

            if (current <= positiveThreshold) {
                positivePeakBins.push_back(i); // push_back some stored 'max peak' bin
                positivePeakValues.push_back(current); // push_back some stored 'max peak' value
                inPositivePeakRegion = false;
            }
        }


        // Negative peak tracking AFTER a positive peak
        if (foundPositivePeak) {

            if (!inNegativePeak && current <= negativeThreshold) {
                // Entering a negative peak region
                inNegativePeak = true;
                maxNegativePeakValue = current;
                maxNegativePeakBin = i; // Provisionally, this might be updated
            } else if (inNegativePeak) {
                if (current < maxNegativePeakValue) {
                    // Found a new max negative value within the peak
                    maxNegativePeakValue = current;
                    maxNegativePeakBin = i;
                } else if (current == maxNegativePeakValue) {
                    // Extend the peak to the last bin with the same value
                    maxNegativePeakBin = i;
                } else if (current > maxNegativePeakValue || i == hist->GetNbinsX()) {
                    // Exiting the negative peak region or end of histogram
                    inNegativePeak = false;
                    negativePeakBins.push_back(maxNegativePeakBin);
                    negativePeakValues.push_back(maxNegativePeakValue);
                    foundPositivePeak = false; // Reset for the next sequence
                    // Prepare for the next negative peak
                    maxNegativePeakValue = std::numeric_limits<double>::max();
                    maxNegativePeakBin = -1;
                }
            }
        }

    } // END OF FOR(histogram bins)

    if (positivePeakBins.size()==0) { // If somehow positive peaks above the threshold weren't found, just get te maximum one.

        double max = hist->GetMaximum();
        double max_bin = hist->GetMaximumBin();

        positivePeakBins.push_back(max_bin);
        positivePeakValues.push_back(max);
    }


    // Print results
//     std::cout << "Positive Peaks: \n";
//     for (size_t i = 0; i < positivePeakBins.size(); ++i) {
//         std::cout << "Bin: " << positivePeakBins[i] << ", Value: " << positivePeakValues[i] << "\n";
//     }
//
//     std::cout << "\nNegative Peaks:\n";
//     for (size_t i = 0; i < negativePeakBins.size(); ++i) {
//         std::cout << "Bin: " << negativePeakBins[i] << ", Value: " << negativePeakValues[i] << "\n";
//     }

//     std::cout << "****************************************\n";


} //END FUNCTION

/////////////////////////////////////////////////////////////////////////////////////////



// BASED ON THE NUMBER OF ENTRIES OF HISTOGRAM: Noisy = millions of entries. Probably not a correct assumption, but good enough for now.
bool isHistogramTooNoisy(TH1D* hist, double maxAllowedEntries){

    if(hist->GetEntries() >= maxAllowedEntries){
        return true;
    }else return false;

}


/////////////////////////////////////////////////////////////////////////////////////////

//  Dealing with noisy histograms:

bool FindPeaksWithNoiseCheck(TH1D* hist, double positiveThreshold, double negativeThreshold,
                             std::vector<int>& positivePeakBins, std::vector<double>& positivePeakValues,
                             std::vector<int>& negativePeakBins, std::vector<double>& negativePeakValues,
                             double maxAllowedEntries){
    // Check if the histogram is too noisy
    if (isHistogramTooNoisy(hist, maxAllowedEntries)) {
        std::cout << "Histogram is too noisy and will be skipped." << std::endl;
        return false; // Histogram was too noisy and was skipped
    }

    // Call the original FindPeaks function if the histogram is not too noisy
    FindPeaks(hist, positiveThreshold, negativeThreshold, positivePeakBins, positivePeakValues, negativePeakBins, negativePeakValues);

    return true; // Histogram was processed
}


/////////////////////////////////////////////////////////////////////////////////////////

// Function to convert TH1F to std::vector<double>
std::vector<double> convert_hist_to_vector(TH1F* hist) {

    std::vector<double> vec(hist->GetNbinsX());
    for (int i = 0; i < hist->GetNbinsX(); ++i) {
        vec[i] = hist->GetBinContent(i + 1);
    }
    return vec;
}
/////////////////////////////////////////////////////////////////////////////////////////


// Function to align waveform. Using the external library FFTW (Fast Fourier Transform in the West) because apparently my ROOT didn't come with FFT tools...
std::vector<double> align_waveform(const std::vector<double>& waveform, double bin_size) {

    int N = waveform.size();
    double sample_rate = 1.0 / bin_size;

    // Setting up some varibles to use with the FFT-er:
    fftw_complex *in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_plan p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    // Copy data to fftw input
    for (int i = 0; i < N; ++i) {

        in[i][0] = waveform[i];
        in[i][1] = 0.0;
    }

    // Perform FFT
    fftw_execute(p);

    // Retrieve real and imaginary parts in the frequency domain:
    double *re = new double[N];
    double *im = new double[N];

    for (int i = 0; i < N; ++i) {

        re[i] = out[i][0];
        im[i] = out[i][1];
    }

    // Calculating magnitudes to identify dominant frequency:

    double *magnitude = new double[N];

    for (int i = 0; i < N; ++i) {
        magnitude[i] = sqrt(re[i] * re[i] + im[i] * im[i]);
    }
    int idx = std::distance(magnitude, std::max_element(magnitude, magnitude + N)); // index of max magnitude.

    double freq = idx * sample_rate / N; // DOminant frequency
    double phase_shift = std::atan2(im[idx], re[idx]); // phase shift

    // Calculate time offset
    double time_offset = phase_shift / (2 * M_PI * freq);

    // Apply phase shift to align the waveform
    for (int i = 0; i < N; ++i) {

        double phase = -2 * M_PI * i * time_offset / N;
        double re_shifted = re[i] * std::cos(phase) - im[i] * std::sin(phase);
        double im_shifted = re[i] * std::sin(phase) + im[i] * std::cos(phase);

        out[i][0] = re_shifted;
        out[i][1] = im_shifted; // Retrieve shifted values for inverse FFT
    }

    // Inverse FFT
    fftw_plan q = fftw_plan_dft_1d(N, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(q);


    // Saving our aligned waveforms:

    std::vector<double> aligned_waveform(N);

    for (int i = 0; i < N; ++i) {

        aligned_waveform[i] = in[i][0] / N;  // Normalize by N
    }

    // Cleanup
    fftw_destroy_plan(p);
    fftw_destroy_plan(q);
    fftw_free(in);
    fftw_free(out);
    delete[] re;
    delete[] im;
    delete[] magnitude;

    return aligned_waveform;
}



/////////////////////////////////////////////////////////////////////////////////////////



// Here is where everything happens:

void processWaveforms() {
    // Open the ROOT file
    TFile *file = new TFile("/home/karla/Documents/WireCell/wf_root/magnify-21040-1.root", "READ");

    // Create files to store the fit results:
    std::ofstream fit_results_file("fit_results_run_" + std::to_string(run) + ".txt");
    fit_results_file << "Ch. # \t t \t A_0 \t t_p \t k3 \t k4 \t k5 \t k6 \t Chi2/DOF" << endl;

    std::ofstream pre_fit_results_file("pre_fit_results_run_" + std::to_string(run) + ".txt");
    pre_fit_results_file << "Ch. # \t t \t A_0 \t t_p \t Chi2/DOF" << endl;

    std::ofstream noisy_hists_file("problems_" + std::to_string(run) + ".txt"); //second fit


    // Retrieve the TH2I histogram
    TH2I *hist2D = (TH2I*)file->Get("h_daq");
    if (!hist2D) {
        std::cerr << "Histogram not found!" << std::endl;
        return;
    }

    // For each channel:

    //for(int i : selectedChannels){
    for(int i = 1; i <= hist2D->GetNbinsX() ; ++i) {

        cout << "\n\n********* Ch. " << (i-1) << " *********\n" << endl;

        // Generate projY histogram name
        TString projYName = TString::Format("projY_ch_%d", (i-1));
        TH1D *projY = hist2D->ProjectionY(projYName, i, i); // add "e" option for errors
        //std::cout << "Number of Bins in ProjY(h_daq):  " << projY->GetNbinsX() << endl;

        // Getting the max and min values of the histogram to use as thresholds:
        double max = projY->GetMaximum();
        double min = projY->GetMinimum();


        // Giving a nice format to plots in case we want to plot raws:
//         projY->SetStats(kFALSE);
//         projY->GetXaxis()->SetTitle("Ticks");
//         projY->GetYaxis()->SetTitle("ADC counts");
//         projY->SetLineWidth(1);
//         projY->SetLineColor(kBlue);
//         projY->GetYaxis()->SetLabelSize(0.025);
//         projY->GetYaxis()->SetTitleSize(0.03);
//         projY->Draw();
//
//         gPad -> Print(Form("plots/raws/run_21040/WF_Raw_Ch_%i_Run_%i.png",(i-1),run));


        // Finding ALL peaks for each histogram:
        double positiveThreshold = max/2.0; // is this the best option for threshold?
        double negativeThreshold = min/2.0;

        std::vector<int> positivePeakBins; //An array with Positive Peak positions
        std::vector<double> positivePeakValues;
        std::vector<int> negativePeakBins; // An array with Negative Peak positions
        std::vector<double> negativePeakValues;

        // Find peaks:
        bool wasProcessed = FindPeaksWithNoiseCheck(projY, positiveThreshold, negativeThreshold, positivePeakBins, positivePeakValues, negativePeakBins, negativePeakValues, 1000000); //1000000 entries is taken as a noisy histogram here.

        if(wasProcessed==0){

            noisy_hists_file << i-1 << "\t (noisy)" << endl;
            continue; // If noisy, skip.

        }

        // Printing info about number and positions of peaks:
        cout << "Number of peaks: " << positivePeakBins.size() << endl;
        cout << "X-positions: " << endl;
        for (size_t j = 0; j < positivePeakBins.size(); ++j) std::cout << positivePeakBins[j] << " ";
        std::cout << std::endl;

        // Creating an array of histograms to store our peaks for each channel:
        std::vector<TH1F*> h_peaks;

        int nbefore = 15; // Number of samples before positive peak
        int nafter = 100; // Number of samples after positive peak

        double x_min; // This would be tick0
        double x_max;
        int nBins;


        for(int j=0; j<positivePeakBins.size(); ++j){ // For each peak, save a histogram

            if(positivePeakBins.size() == 1){ // If there is just one peak

                x_min = positivePeakBins[0]- nbefore;
                x_max = positivePeakBins[0] + nafter;

            }else{

                x_min = positivePeakBins[j]- nbefore;
                x_max = positivePeakBins[j] + nafter;

            }

            // NUmber of bins  in our individual waveforms/peaks:
            nBins = abs(x_max - x_min + 1);

            // Creating an array of histograms to store our peaks:
            h_peaks.push_back(new TH1F(Form("hist_%d_channel_%i", j, i), Form("Pulse %i, Channel %i", j, i), nBins, x_min, x_max));

            // Filling up individual pulse histograms:
            for (int k = x_min; k <= x_max; ++k) h_peaks[j]->Fill(k, projY->GetBinContent(k));

                /*h_peaks[j]->SetLineColor(kBlue);
                h_peaks[j]->GetXaxis()->SetTitle("Ticks");
                h_peaks[j]->GetYaxis()->SetTitle("ADC counts");
                h_peaks[j]->SetLineWidth(1);
                h_peaks[j]->SetLineColor(4);
                h_peaks[j]->GetYaxis()->SetLabelSize(0.025);
                h_peaks[j]->GetYaxis()->SetTitleSize(0.03);
                h_peaks[j]->Draw("HIST");
                gPad -> Print(Form("plots/singles_per_channel/run_21040/Ch_%i_WF_%i_Run_%i.png",(i-1),(j+1),run));*/

        } // After this FOR, individual pulse histograms have been created and stored in h_peaks[j].


        double bin_size = 0.512; // To calculate the sample rate.
        std::vector<std::vector<double>> aligned_waveforms; // To store the aligned pulses

        for (auto& hist : h_peaks) {// Aligning peaks:

            std::vector<double> waveform = convert_hist_to_vector(hist);
            std::vector<double> aligned_waveform = align_waveform(waveform, bin_size);

            aligned_waveforms.push_back(aligned_waveform);
        }

        // Averaging the aligned waveforms:

        int N = aligned_waveforms[0].size(); // NUmber of peaks in each raw.
        std::vector<double> averaged_waveform(N, 0.0); // To store averaged waveforms.

        for (const auto& aligned_waveform : aligned_waveforms) { // Adding entries
            for (int k = 0; k < N; ++k) {

                averaged_waveform[k] += aligned_waveform[k];
            }
        }

        for (int k = 0; k < N; ++k) { // Take the average

            averaged_waveform[k] /= aligned_waveforms.size();
        }

        // Let's plot them to see if this worked!

        TH1D* h_averaged = new TH1D(Form("h_averaged_ch_%d", i-1), Form("Averaged Waveform Channel %d", i-1), N, 0, N);

        for (int k = 0; k < N; ++k) {

            h_averaged->SetBinContent(k + 1, averaged_waveform[k]);
        }

        // Scaling X axis to be microseconds instead of ticks (512ns per tick)
        ScaleXaxis(h_averaged, ScaleX);
        h_averaged->GetXaxis()->SetRangeUser(0,70);

        h_averaged->SetLineColor(kBlue);
        h_averaged->GetXaxis()->SetTitle("Time (#mus)");
        h_averaged->GetYaxis()->SetTitle("ADC counts");
        h_averaged->SetLineWidth(1);
        h_averaged->SetLineColor(4);
        h_averaged->GetYaxis()->SetLabelSize(0.025);
        h_averaged->GetYaxis()->SetTitleSize(0.03);
        //h_averaged->Draw();
        //gPad -> Print(Form("plots/averaged/run_21040/averaged_WF_ch_%d_Run_%d.png", i-1,run));



        /////////////////////////////////////////////////////////////////////////////////
        // FITTING THE WAVEFORMS (IN THE WORKS!!):
        /////////////////////////////////////////////////////////////////////////////////

        /*
        // Pre fitting: Using old electronics response provided by Wenqiang.
        TF1 *f_ideal = new TF1("f_ideal", response_legacy, 0,100,3);
        f_ideal->SetNpx(1000);

        // Fitting using response provided by Xin:
        TF1 *f_response = new TF1("f_response", response, 0,100,7);
        //(t, A0, tp, k3, k4, k5, k6)
        //t = x[0] - par[0], A0 = par[1], tp = par[2].
        f_response->SetNpx(1000);


        f_ideal -> SetParameters(4.0, max, 2.2);
        // Gaussian, using old electronics response (t,A0,t_p).
        TFitResultPtr r0 =  h_averaged -> Fit("f_ideal","QS","",0,10.0);
//         cout << "Covariance Matrix Status: " << r0->CovMatrixStatus() << endl;
//         cout << "Probability of Chi2: " << r0->Prob() << endl;



        double chi2 = r0->Chi2();
        double dof = r0->Ndf();
        double chi2PerDof = dof > 0 ? chi2 / dof : 0; // Avoid division by zero

        // Writing pre-fit results in a file to plot them later:
        pre_fit_results_file << i-1 << "\t" << setprecision(4) << r0->Parameter(0) << "\t" << setprecision(5) << r0->Parameter(1) << "\t" << setprecision(4) << r0->Parameter(2) << "\t" << setprecision(4) << chi2PerDof << endl;


        if (min == 0 || i==1111 || i==1747 || i==2007 || i==2111 || i==2974) { // For channels with zero tails, OR WHERE FIT JUST DOESN'T WORK, try different parameters

            //cout << "** Channel with no negative peaks **" << endl;
            noisy_hists_file << i-1 << "\t (zero baseline)" << endl;

            //f_response->FixParameter(3,0.005);
            f_response->SetParameters(5.5,5500,2.2,0.1,0.1,0.0,0.0);

            f_response -> SetParName(0,"t");
            f_response -> SetParName(1,"A_{0}");
            f_response -> SetParName(2,"t_{p}");
            f_response -> SetParName(3,"k_{3}");
            f_response -> SetParName(4,"k_{4}");
            f_response -> SetParName(5,"k_{5}");
            f_response -> SetParName(6,"k_{6}");

            TFitResultPtr r = h_averaged->Fit(f_response,"QRS","",0.0,50.0);

            double chi2 = r->Chi2();
            double dof = r->Ndf();
            double chi2PerDof = dof > 0 ? chi2 / dof : 0;

            h_averaged->SetLineColor(38);
            h_averaged->Draw();

            // Create a new stats box if not found
            TPaveStats *st = (TPaveStats*)h_averaged->GetListOfFunctions()->FindObject("stats");
            if (!st) {
                st = new TPaveStats(0.58, 0.52, 0.88, 0.87, "brNDC");
                st->SetName("stats");
                st->SetBorderSize(1);
                h_averaged->GetListOfFunctions()->Add(st);
            }

            st->SetOptStat(1110); // Set the options for what to display
            st->SetOptFit(111);    // Show fit parameters as well

            st->SetTextSize(0.03);
            st->SetX1NDC(0.58); // new x start position
            st->SetX2NDC(0.88); // new x end position
            st->SetY1NDC(0.87); // new y start position
            st->SetY2NDC(0.52); // new y end position

            gPad->Modified();
            gPad -> Print(Form("plots/fitted/run_21040/zero_baseline/Fit_2_Ch_%i_Run_%i_VD.png",(i-1),run));
            gPad->Update();

            fit_results_file << i-1 << "\t" << setprecision(4) << r->Parameter(0) << "\t" << setprecision(5) << r->Parameter(1) << "\t" << setprecision(4) << r->Parameter(2) << "\t" << setprecision(4) << r->Parameter(3) << "\t"<< setprecision(4) << r->Parameter(4) << "\t" << setprecision(4) << r->Parameter(5) << "\t" << setprecision(4) << r->Parameter(6) << "\t" << setprecision(4) << chi2PerDof << "\t" << f_response->Eval(50) <<  endl;

        }else{

            //Feed pre-fit parameters to response function:
            f_response->SetParameters(f_ideal->GetParameter(0),f_ideal->GetParameter(1),f_ideal->GetParameter(2),0.1,0.1,0.03,0.03);
            //(t, A0, tp, k3, k4, k5, k6)

            // Setting some parameters for fits that don't want to work well:
            f_response->SetParLimits(1,f_ideal->GetParameter(1)*0.9,f_ideal->GetParameter(1)*1.1); //chi2 ~ 2.5e4


            f_response -> SetParName(0,"t");
            f_response -> SetParName(1,"A_{0}");
            f_response -> SetParName(2,"t_{p}");
            f_response -> SetParName(3,"k_{3}");
            f_response -> SetParName(4,"k_{4}");
            f_response -> SetParName(5,"k_{5}");
            f_response -> SetParName(6,"k_{6}");


            TFitResultPtr r = h_averaged->Fit(f_response,"QRS","",0,50.0);

            double chi2 = r->Chi2();
            double dof = r->Ndf();
            double chi2PerDof = dof > 0 ? chi2 / dof : 0;


            // Saving fit results into a file:
            fit_results_file << i-1 << "\t" << setprecision(4) << r->Parameter(0) << "\t" << setprecision(5) << r->Parameter(1) << "\t" << setprecision(4) << r->Parameter(2) << "\t" << setprecision(4) << r->Parameter(3) << "\t"<< setprecision(4) << r->Parameter(4) << "\t" << setprecision(4) << r->Parameter(5) << "\t" << setprecision(4) << r->Parameter(6) << "\t" << setprecision(4) << chi2PerDof << "\t" << f_response->Eval(50) <<  endl;



            h_averaged->SetLineColor(38);
            h_averaged->Draw();

            // Create a new stats box if not found
            TPaveStats *st = (TPaveStats*)h_averaged->GetListOfFunctions()->FindObject("stats");
            if (!st) {
                st = new TPaveStats(0.58, 0.52, 0.88, 0.87, "brNDC");
                st->SetName("stats");
                st->SetBorderSize(1);
                h_averaged->GetListOfFunctions()->Add(st);
            }

            st->SetOptStat(1110); // Set the options for what to display
            st->SetOptFit(111);    // Show fit parameters as well

            st->SetTextSize(0.03);
            st->SetX1NDC(0.58); // new x start position
            st->SetX2NDC(0.88); // new x end position
            st->SetY1NDC(0.87); // new y start position
            st->SetY2NDC(0.52); // new y end position

            gPad->Modified();
            gPad -> Print(Form("plots/fitted/run_21040/Fit_2_Ch_%i_Run_%i_VD.png",(i-1),run));
            gPad->Update();



        }


//         // Draw random fitted hists:
//         if (std::find(randomIntegers.begin(), randomIntegers.end(), i) != randomIntegers.end()) {
//
//             h_averaged->SetLineColor(38);
//             h_averaged->Draw();
//
//             // Create a new stats box if not found
//             TPaveStats *st = (TPaveStats*)h_averaged->GetListOfFunctions()->FindObject("stats");
//             if (!st) {
//                 st = new TPaveStats(0.58, 0.52, 0.88, 0.87, "brNDC");
//                 st->SetName("stats");
//                 st->SetBorderSize(1);
//                 h_averaged->GetListOfFunctions()->Add(st);
//             }
//
//             st->SetOptStat(1110); // Set the options for what to display
//             st->SetOptFit(111);    // Show fit parameters as well
//
//             st->SetTextSize(0.03);
//             st->SetX1NDC(0.58); // new x start position
//             st->SetX2NDC(0.88); // new x end position
//             st->SetY1NDC(0.87); // new y start position
//             st->SetY2NDC(0.52); // new y end position
//
//             gPad->Modified();
//             gPad->Update();
//             gPad->Print(Form("fitted/Fitted_Ch_%i_Run_%i_HD.png",(i-1),run));
//
//
//          }


        // Clean up
         delete f_response;
         delete f_ideal;

         */
         delete h_averaged;

        //std::cout << "End of Analysis for Channel " << (i-1) << " \n" << endl;

        /////////////////////////////////////////////////////////////////////////////////


    } //END FOR Channels


    // Clean up and close files
    file->Close();

    pre_fit_results_file.close();
    fit_results_file.close();
    noisy_hists_file.close();


    delete file;


} // END OF Function









