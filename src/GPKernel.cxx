#include "WCPLEEANA/GPKernel.h"

#include <math.h>
#include <iostream>

//----------------------------------------------------------------------
TMatrixDSym GPKernel::operator()(std::vector<GPPoint> pts, int dpar_idx) const
{
  TMatrixDSym sym(pts.size());
  for(int i = 0; i < (int)pts.size(); i++){
    for(int j = 0; j <= i; j++){
      double element = Element(pts[i], pts[j], dpar_idx);
      sym(i, j) = element;
      sym(j, i) = element;
    }
  }
  return sym;
}

//----------------------------------------------------------------------
TMatrixD GPKernel::KernelBlock(std::vector<GPPoint> pts1, std::vector<GPPoint> pts2) const
{
  TMatrixD block(pts1.size(), pts2.size());
  for(int i = 0; i < (int)pts1.size(); i++){
    for(int j = 0; j < (int)pts2.size(); j++){
      block(i, j) = Element(pts1[i], pts2[j]);
    }
  }
  return block;
}

//----------------------------------------------------------------------
TVectorD GPKernel::KernelDiag(std::vector<GPPoint> pts) const
{
  TVectorD diag(pts.size());
  for(int i = 0; i < (int)pts.size(); i++)
    diag[i] = Element(pts[i], pts[i]);
  return diag;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
double RBFKernel::Mag(GPPoint p1, GPPoint p2) const {
  double* p1x = p1.X();
  double* p2x = p2.X();
  double mag = 0;

  // Guard uses original coordinates before any log transform (B-02 fix).
  // Zero-length-scale dimensions with different values → infinite separation.
  for (int i=0;i<5;i++) { if (GPKernel::fPars[i]==0 && p1x[i]!=p2x[i]) { return 1e6; } }

  std::vector<double> parameters;
  double lp1x[5], lp2x[5];
  //use log scales for fractional smoothing (ie: 20%)
  for (int i=0;i<5;i++) {
    if (GPKernel::fPars[i]==0) {
      // Zero length scale: identical values guaranteed by guard above; skip.
      lp1x[i] = lp2x[i] = 0.0;
      parameters.push_back(1.0);  // unused; division avoided by lp1x==lp2x
    } else if (doLogScales[i]) {
      lp1x[i] = log(p1x[i]);
      lp2x[i] = log(p2x[i]);
      parameters.push_back(log(fPars[i]));
    } else {
      lp1x[i] = p1x[i];
      lp2x[i] = p2x[i];
      parameters.push_back(fPars[i]);
    }
  }
  //compute the non-euclidean distance between two points
  for (int i=0;i<5;i++) { if (lp1x[i]!=lp2x[i]) { mag += TMath::Power((lp1x[i]-lp2x[i])/parameters[i], 2); } }

  return mag;
};
//----------------------------------------------------------------------
double RBFKernel::Element(GPPoint pt1, GPPoint pt2, int dpar_idx) const
{
   assert(dpar_idx < 5 && "RBF Kernel has 6 hyperparameters, index is either -1 or 0-4");
   double coeff = fPars[5];
  
  double mag = RBFKernel::Mag(pt1,pt2);
  if(mag == 0.){
    if(dpar_idx >= 0 && dpar_idx < 5) {
      return 0.;
    } else {
      return coeff * (1 + fNoise);
    }
  }

  if(dpar_idx >=0 && dpar_idx < 5){
    double lscale = doLogScales[dpar_idx] ? log(fPars[dpar_idx]) : fPars[dpar_idx];
    if(doLogScales[dpar_idx]) {
      return coeff*exp(-mag/2)*TMath::Power((pt1.X()[dpar_idx] - pt2.X()[dpar_idx]),2)/TMath::Power(lscale,3);
    } else {
       return coeff*exp(-mag/2)*TMath::Power((pt1.X()[dpar_idx] - pt2.X()[dpar_idx]),2)/TMath::Power(lscale,2);
    }
  } else {
    return coeff*exp(-mag/2);
  }
}

//----------------------------------------------------------------------
void RBFKernel::SetThetas(std::vector<double> thetas)
{
  assert(thetas.size() == fPars.size() && "Theta vector should be of size 6 for RBF Kernel");
  fPars[0] = exp(thetas[0]);
  fPars[1] = exp(thetas[1]);
  fPars[2] = exp(thetas[2]);
  fPars[3] = exp(thetas[3]);
  fPars[4] = exp(thetas[4]);
  fPars[5] = exp(thetas[5]);
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
double RationalQuadraticKernel::Element(GPPoint pt1, GPPoint pt2, int dpar_idx) const
{
  assert(dpar_idx < 3 && "Rational Quadratic Kernel has only 3 hyperparameters, index is either -1, 0 or 1");
  
  double scale = fPars[0];
  double alpha = fPars[1];
  double coeff = fPars[2];
  if(Mag(pt1,pt2) == 0.){
    if(dpar_idx >= 0 && dpar_idx < 2) 
      return 0.;
    else
      return coeff * (1. + fNoise);
  }

  double dists = Mag(pt1,pt2);
  double base = dists/(2 * alpha * TMath::Power(scale, 2));
  base += 1.;

  double val = TMath::Power(base, -alpha);

    if(dpar_idx == 0)
    return coeff * dists * val/(base * TMath::Power(scale, 2));
  else if(dpar_idx == 1)
    return coeff * val * (-alpha * log(base) + dists / (2 * TMath::Power(scale, 2) * base));
  else
    return coeff * val;
}

//----------------------------------------------------------------------
void RationalQuadraticKernel::SetThetas(std::vector<double> thetas)
{
  assert(thetas.size() == fPars.size() && "Theta vector should be of size 3 for Rational Quadratic Kernel");
  fPars[0] = exp(thetas[0]);
  fPars[1] = exp(thetas[1]);
  fPars[2] = exp(thetas[2]);
}
