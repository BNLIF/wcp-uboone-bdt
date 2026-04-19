#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TDecompSVD.h"

#include "TMath.h"
#include "TFile.h"
#include <iostream>
#include <fstream>
#include <string>

const int nbins1_2D    = 1;	//1 Enu bins
const int nbins2_2D    = 9;	//9 theta bins
const int nbins3_2D    = 6;	//6 Pmu bins
const int nbins_tot_2D = 36;

const int nbins1_3D    = 4;	//4 Enu bins
const int nbins2_3D    = 9;	//9 theta bins
const int nbins3_3D    = 9;	//6 Pmu bins
const int nbins_tot_3D = 138;

std::vector<int> nbins_2D = { nbins1_2D, nbins2_2D, nbins3_2D };
std::vector<int> nbins_3D = { nbins1_3D, nbins2_3D, nbins3_3D };

//maps from signal bin index to Enu slice (not used in 2D)
int map_signal_dim1_2D[nbins_tot_2D] = { 0,0,0,0,0, 0,0,0,0,0,  0,0,0,0,0, 0,0,0,0,0,  0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0 };

//maps from signal bin index to costheta slice
int map_signal_dim2_2D[nbins_tot_2D] = { 0,0,0,
                                         1,1,1,1,
                                         2,2,2,2,
                                         3,3,3,
                                         4,4,4,4,
                                         5,5,5,5,
                                         6,6,6,6,
                                         7,7,7,7,7,
                                         8,8,8,8,8 };

//maps from signal bin index to Pmu min slice
int map_signal_dim3_2D[nbins_tot_2D] = { 0,1,2,
                                         0,1,2,3,
                                         0,1,2,3,
                                         0,  2,3,
                                         0,  2,3,4,
                                         0,  2,3,4,
                                         0,  2,3,4,
                                         0,  2,3,4,5,
                                         0,  2,3,4,5 };

//maps from signal bin index to Enu slice
int map_signal_dim1_3D[nbins_tot_3D] = { 0,0,0,0,0, 0,0,0,0,0,  0,0,0,0,0, 0,0,0,0,0,  0,0,0,0,0, 0,0,0,
                                         1,1,1,1,1, 1,1,1,1,1,  1,1,1,1,1, 1,1,1,1,1,  1,1,1,1,1, 1,1,1,1,1, 1,1,1,1,1,
                                         2,2,2,2,2, 2,2,2,2,2,  2,2,2,2,2, 2,2,2,2,2,  2,2,2,2,2, 2,2,2,2,2, 2,2,2,2,2, 2,2,2,2,2, 2,2,
                                         3,3,3,3,3, 3,3,3,3,3,  3,3,3,3,3, 3,3,3,3,3,  3,3,3,3,3, 3,3,3,3,3, 3,3,3 };

//maps from signal bin index to costheta slice
int map_signal_dim2_3D[nbins_tot_3D] = { 0,0,0,
                                         1,1,1,
                                         2,2,2,2,
                                         3,3,3,
                                         4,4,4,
                                         5,5,5,
                                         6,6,6,
                                         7,7,7,
                                         8,8,8,

                                         0,0,0,
                                         1,1,1,1,
                                         2,2,2,2,
                                         3,3,3,
                                         4,4,4,
                                         5,5,5,5,5,
                                         6,6,6,6,6,
                                         7,7,7,7,7,
                                         8,8,8,

                                         0,0,0,
                                         1,1,1,1,
                                         2,2,2,2,
                                         3,3,3,
                                         4,4,4,4,4,
                                         5,5,5,5,5,5,
                                         6,6,6,6,6,6,
                                         7,7,7,7,7,7,
                                         8,8,8,8,8,

                                         0,0,0,
                                         1,1,1,1,
                                         2,2,2,
                                         3,3,3,
                                         4,4,4,4,
                                         5,5,5,5,
                                         6,6,6,
                                         7,7,7,7,
                                         8,8,8,8,8   };

//maps from signal bin index to Pmu min slice
int map_signal_dim3_3D[nbins_tot_3D] = { 0,1,2,	//gives the Pmu bin low for each 138 bins
                                         0,1,2,
                                         0,1,2,3,
                                         0,  2,3,
                                         0,  2,3,
                                         0,  2,3,
                                         0,  2,3,
                                         0,  2,3,
                                         0,  2,3,

                                         0,1,2,
                                         0,1,2,3,
                                         0,1,2,3,
                                         0,  2,3,
                                         0,  2,3,
                                         0,  2,3,4,5,
                                         0,  2,3,4,5,
                                         0,  2,3,4,5,
                                         0,    3,  5,

                                         0,1,2,
                                         0,1,2,3,
                                         0,1,2,3,
                                         0,  2,3,
                                         0,  2,3,4,5,
                                         0,  2,3,4,5,6,
                                         0,  2,3,4,5,6,
                                         0,  2,3,  5,6,7,
                                         0,    3,  5,6,7,

                                         0,1,2,
                                         0,1,2,3,
                                         0,  2,3,
                                         0,  2,3,
                                         0,  2,3,  5,
                                         0,  2,3,  5,
                                         0,    3,  5,
                                         0,    3,  5,  7,
                                         0,        5,6,7,8 };

//maps from {Enu, costheta, Pmu} indices to signal bin index
std::vector<std::vector<std::vector<int>>> map_grid_signal_bin_2D = { { {  0,  1,  2,  2,  2,  2, },
                                                                        {  3,  4,  5,  6,  6,  6, },
                                                                        {  7,  8,  9, 10, 10, 10, },
                                                                        { 11, 11, 12, 13, 13, 13, },
                                                                        { 14, 14, 15, 16, 17, 17, },
                                                                        { 18, 18, 19, 20, 21, 21, },
                                                                        { 22, 22, 23, 24, 25, 25, },
                                                                        { 26, 26, 27, 28, 29, 30, },
                                                                        { 31, 31, 32, 33, 34, 35 } } };

//maps from {Enu, costheta, Pmu} indices to signal bin index
std::vector<std::vector<std::vector<int>>> map_grid_signal_bin_3D = { { {  0,  1,  2,  2,  2,  2,  2,  2,  2 },
                                                                        {  3,  4,  5,  5,  5,  5,  5,  5,  5 },
                                                                        {  6,  7,  8,  9,  9,  9,  9,  9,  9 },
                                                                        { 10, 10, 11, 12, 12, 12, 12, 12, 12 },
                                                                        { 13, 13, 14, 15, 15, 15, 15, 15, 15 },
                                                                        { 16, 16, 17, 18, 18, 18, 18, 18, 18 },
                                                                        { 19, 19, 20, 21, 21, 21, 21, 21, 21 },
                                                                        { 22, 22, 23, 24, 24, 24, 24, 24, 24 },
                                                                        { 25, 25, 26, 27, 27, 27, 27, 27, 27 } },

                                                                        { { 28, 29, 30, 30, 30, 30, 30, 30, 30 },
                                                                        { 31, 32, 33, 34, 34, 34, 34, 34, 34 },
                                                                        { 35, 36, 37, 38, 38, 38, 38, 38, 38 },
                                                                        { 39, 39, 40, 41, 41, 41, 41, 41, 41 },
                                                                        { 42, 42, 43, 44, 44, 44, 44, 44, 44 },
                                                                        { 45, 45, 46, 47, 48, 49, 49, 49, 49 },
                                                                        { 50, 50, 51, 52, 53, 54, 54, 54, 54 },
                                                                        { 55, 55, 56, 57, 58, 59, 59, 59, 59 },
                                                                        { 60, 60, 60, 61, 61, 62, 62, 62, 62 } },

                                                                      { { 63, 64, 65, 65, 65, 65, 65, 65, 65 },
                                                                        { 66, 67, 68, 69, 69, 69, 69, 69, 69 },
                                                                        { 70, 71, 72, 73, 73, 73, 73, 73, 73 },
                                                                        { 74, 74, 75, 76, 76, 76, 76, 76, 76 },
                                                                        { 77, 77, 78, 79, 80, 81, 81, 81, 81 },
                                                                        { 82, 82, 83, 84, 85, 86, 87, 87, 87 },
                                                                        { 88, 88, 89, 90, 91, 92, 93, 93, 93 },
                                                                        { 94, 94, 95, 96, 96, 97, 98, 99, 99 },
                                                                        {100,100,100,101,101,102,103,104,104 } },

                                                                      { {105,106,107,107,107,107,107,107,107 },
                                                                        {108,109,110,111,111,111,111,111,111 },
                                                                        {112,112,113,114,114,114,114,114,114 },
                                                                        {115,115,116,117,117,117,117,117,117 },
                                                                        {118,118,119,120,120,121,121,121,121 },
                                                                        {122,122,123,124,124,125,125,125,125 },
                                                                        {126,126,126,127,127,128,128,128,128 },
                                                                        {129,129,129,130,130,131,131,132,132 },
                                                                        {133,133,133,133,133,134,135,136,137 } } };
 

//gets the index for the next (or previous) bin along a given dimension.  Returns -1 if a boundary is hit.
int get_next_bin (int ndim, int start_bin, int dim, bool direction) {
  if (start_bin==-1) { return -1; }

  int nbins_dim = -1;
  std::vector<int> bin_3D = {};
  if      (ndim==2) {
    nbins_dim = nbins_2D[dim];
    bin_3D = { map_signal_dim1_2D[start_bin], map_signal_dim2_2D[start_bin], map_signal_dim3_2D[start_bin] };
  }
  else if (ndim==3) {
    nbins_dim = nbins_3D[dim];
    bin_3D = { map_signal_dim1_3D[start_bin], map_signal_dim2_3D[start_bin], map_signal_dim3_3D[start_bin] };
  }
//std::cout << "bin_3D = {" << bin_3D[0] << ", " << bin_3D[1] << ", " << bin_3D[2] << "}" << std::endl;
  while (true) {
    //increase/decrease bin[dim] and check whether you are out of bounds
    bin_3D[dim] = bin_3D[dim] + 1*direction - 1*(!direction);
    if (bin_3D[dim]<0 || bin_3D[dim]>=nbins_dim) { return -1; }

    //check whether the bin is different, otherwise continue
    int newbin = -1;
    if      (ndim==2) { newbin = map_grid_signal_bin_2D[bin_3D[0]][bin_3D[1]][bin_3D[2]]; }
    else if (ndim==3) { newbin = map_grid_signal_bin_3D[bin_3D[0]][bin_3D[1]][bin_3D[2]]; }
    if (newbin != start_bin) { return newbin; }
  }
  return -1;
}

//returns how many bins exist along the given dimension and direction
int get_nbins_next (int ndim, int start_bin, int dim, bool direction) {
  if (start_bin==-1) { return -1; }
std::cout << "ndim, start_bin, dim, dir = " << ndim << ",  " << start_bin << ",  " << dim << ",  " << direction << std::endl;

  int current_bin = start_bin;
  int nbins = -1;  
  while(current_bin!=-1) {
    nbins++;
std::cout << "current bin = " << current_bin << std::endl;
    current_bin = get_next_bin(ndim,current_bin,dim,direction);
  }
  return nbins;
}

/*
  void compute_2nd_derivatives (TMatrixD* C, int start_bin, int j, int k) {

      int prev_bin_j = get_next_bin(start_bin, j, false);
      int next_bin_j = get_next_bin(start_bin, j, true);
      int next_bin_k = get_next_bin(start_bin, k, true);

      Double_t epsilon = 1e-6;

    //second derivative in the same direction
    if (j==k) {
      if (prev_bin_j!=-1) {
        C_temp(start_bin,prev_bin_j) = C_temp(start_bin,prev_bin_j) + 1;
        C_temp(start_bin,start_bin ) = C_temp(start_bin,start_bin ) - 1 + epsilon;
      }
      if (next_bin_j!=-1) {
        C_temp(start_bin,next_bin_j) = C_temp(start_bin,next_bin_j) + 1;
        C_temp(start_bin,start_bin ) = C_temp(start_bin,start_bin ) - 1 + epsilon;
      }
    //second derivative in two separate directions
    } else {

      C_temp(start_bin,start_bin) = C_temp(start_bin,start_bin) + 1;
      if (next_bin_j!=-1) {
        C_temp(start_bin,next_bin_j) = C_temp(start_bin,next_bin_j) - 1 + epsilon;
        //std::cout << "start_bin, next_bin_j = " << start_bin << ",     " << next_bin_j << std::endl;
        int next_bin_j_next_bin_k = get_next_bin(next_bin_j, k, true);

        //std::vector<int> bin_3D = { map_signal_dim1_3D[next_bin_j], map_signal_dim2_3D[next_bin_j], map_signal_dim3_3D[next_bin_j] };
        //std::cout << "looking for bin after " << next_bin_j << ",  (" << bin_3D[0] << ", " << bin_3D[1] << ", " << bin_3D[2] << ")" << std::endl;
        if (next_bin_j_next_bin_k!=-1) {
          //std::cout << "found!" << std::endl;
          C_temp(start_bin,next_bin_j_next_bin_k) = C_temp(start_bin,next_bin_j_next_bin_k) + 0.5;
        }
      }
      if (next_bin_k!=-1) {
        C_temp(start_bin,next_bin_k) = C_temp(start_bin,next_bin_k) - 1 + epsilon;
        int next_bin_k_next_bin_j = get_next_bin(next_bin_k, j, true);
        if (next_bin_k_next_bin_j!=-1) {
          C_temp(start_bin,next_bin_k_next_bin_j) = C_temp(start_bin,next_bin_k_next_bin_j) + 0.5;
        }
      }
    }

  }
*/

//main function to compute C3 regularization matrix using the full 3D bin structure
TMatrixD C3_3D (int ndim) {

  Double_t epsilon = 1e-2;
  int nbins_tot = -1;
  if      (ndim==2) { nbins_tot = nbins_tot_2D; }
  else if (ndim==3) { nbins_tot = nbins_tot_3D; }

  TMatrixDSym C_sym(nbins_tot);
  for (int i=0;i<nbins_tot;i++) {
    for (int j=0;j<nbins_tot;j++) {
      C_sym(i,j) = 0;
    }
  }

  //take derivative along Enu, theta, and Pmu
  int dim_max = 3;
  int dim_min = dim_max-ndim;
  for (int dim=dim_min;dim<dim_max;dim++) {
    TMatrixD C_temp(nbins_tot,nbins_tot);

    //iterate over 36 (2D) or 138 (3D) truth signal bins and compute 3rd derivative for each
    for (int start_bin=0;start_bin<nbins_tot;start_bin++) {
      //set matrix to 0 initially
      for (int j=0;j<nbins_tot;j++) { C_temp(start_bin,j) = 0; }

      int prev_bin  = get_next_bin(ndim,start_bin, dim, false);
      int prev_bin2 = get_next_bin(ndim,prev_bin,  dim, false);
      int next_bin  = get_next_bin(ndim,start_bin, dim, true);
      int next_bin2 = get_next_bin(ndim,next_bin,  dim, true);

      int nbins_prev = get_nbins_next(ndim,start_bin,dim,false);
      int nbins_next = get_nbins_next(ndim,start_bin,dim,true);

      //Old C3 matrix
      if (nbins_prev>=2 && nbins_next>=2) {
        C_temp(start_bin,prev_bin2) = C_temp(start_bin,prev_bin2) - 1;
        C_temp(start_bin,prev_bin)  = C_temp(start_bin,prev_bin)  + 2;
        C_temp(start_bin,start_bin) = C_temp(start_bin,start_bin) + epsilon;
        C_temp(start_bin,next_bin)  = C_temp(start_bin,next_bin)  - 2;
        C_temp(start_bin,next_bin2) = C_temp(start_bin,next_bin2) + 1;
      } else if (nbins_prev==1 && nbins_next==1) {
        C_temp(start_bin,start_bin) = C_temp(start_bin,start_bin) + epsilon;
      } else if (nbins_next==0 || (nbins_next==1 && nbins_prev>=2)) {
        C_temp(start_bin,prev_bin2) = C_temp(start_bin,prev_bin2) - 1;
        C_temp(start_bin,prev_bin)  = C_temp(start_bin,prev_bin)  + 2;
        C_temp(start_bin,start_bin) = C_temp(start_bin,start_bin) - 1 + epsilon;
      } else if (nbins_prev==0 || (nbins_prev==1 && nbins_next>=2)) {
        C_temp(start_bin,start_bin) = C_temp(start_bin,start_bin) + 1 + epsilon;
        C_temp(start_bin,next_bin)  = C_temp(start_bin,next_bin)  - 2;
        C_temp(start_bin,next_bin2) = C_temp(start_bin,next_bin2) + 1;
      } else {
        std::cout << "error here" << std::endl;
      }
    }
std::cout << "here 2" << std::endl;

    //square C_temp for each dimension and then add to C_sym
    TMatrixD C_temp_t(nbins_tot,nbins_tot);
    C_temp_t.Transpose(C_temp);
    TMatrixD C_temp_2 = C_temp_t * C_temp;
    for (int i=0;i<nbins_tot;i++) {
      for (int j=0;j<nbins_tot;j++) {
        C_sym(i,j) = C_sym(i,j) + C_temp_2(i,j);
      }
    }

  }
std::cout << "here 3" << std::endl;

  //Decompose into C = V*D*V^T
  TMatrixDSymEigen C_eigen(C_sym);
  TMatrixD V = C_eigen.GetEigenVectors();
  TMatrixD V_t(nbins_tot,nbins_tot);
  V_t.Transpose(V);

 //take square root of D
  TVectorD V_eigen = C_eigen.GetEigenValues();
  TMatrixD D(nbins_tot,nbins_tot);
  for (int i=0;i<nbins_tot;i++) {
    for (int j=0;j<nbins_tot;j++) { D(i,j) = 0; }
    D(i,i) = TMath::Sqrt(V_eigen[i]);
  }

  TMatrixD C = V * D * V_t;

/*
  TFile *file = new TFile("/uboone/data/users/lcoopert/LEE/LEEana_xs_3D_feb24/wiener_svd/debug_wienerSVD_3D.root","RECREATE");
  file->cd();
  C.Write("C");
  file->Write();
  file->Close();
*/
  return C;
}

