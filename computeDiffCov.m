%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    MIT License
%  
%    Copyright (c) 2020 Abhinav Prakash
%  
%    Permission is hereby granted, free of charge, to any person obtaining a copy
%    of this software and associated documentation files (the "Software"), to deal
%    in the Software without restriction, including without limitation the rights
%    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
%    copies of the Software, and to permit persons to whom the Software is
%    furnished to do so, subject to the following conditions:
%  
%    The above copyright notice and this permission notice shall be included in all
%    copies or substantial portions of the Software.
%  
%    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
%    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
%    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
%    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
%    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
%    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
%    SOFTWARE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   AUTHOR: ABHINAV PRAKASH
%   DESCRIPTION:
%   Function to compute the difference in mean vector for two GPs, 
%   and the covariance matrix for the difference GP under null hypothesis.
%   INPUTS:
%       x1: Input variable matrix or vector for first set of points
%       y1: response vector for first set of points
%       x2: Input variable matrix or vector for second set of points   
%       y2: response vector for second set of points
%       xtest: Input variable matrix or vector for the test points at which the mean vector 
%              and covariance matrix of the difference has to be computed
%       params: Hyperparamters for computing the difference mean and covariance 
%   OUTPUT: A data structure with the following elements
%       mu: mean difference vector mu2 - mu1
%       K:  difference covariance matrix      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[out] = computeDiffCov(x1,y1,x2,y2,xtest,params)
   %accessing the hyperparameters from the data structure params
   theta = params.theta;
   sigma_f = params.sigma_f;
   sigma_n = params.sigma_n;
   beta = params.beta;

   KX1X1 = (sigma_f^2)*computeCorrelMat(x1,x1,theta); %covariance matrix for first function
   KX1X1 = KX1X1 + (eye(length(y1))*(sigma_n^2)); % adding noise
   opts.POSDEF = true;
   opts.SYM = true;
   invKX1X1 = linsolve(KX1X1,eye(size(KX1X1)),opts); %computing inverse
   clear KX1X1; %freeing memory
   KXTX1 = (sigma_f^2)*computeCorrelMat(xtest,x1,theta); %computing covariance matrix for test points and the first dataset
   mu1 = beta + (KXTX1*(invKX1X1*(y1-beta))); %computing mean vector conditioned on the first dataset 
   K1 = invKX1X1 * (KXTX1'); 
   K = KXTX1 * K1; %posterior covariance matrix for the first function
   clear invKX1X1 KXTX1;
   
   KX2X2 = (sigma_f^2)*computeCorrelMat(x2,x2,theta); %covariance matrix for second function
   KX2X2 = KX2X2 + (eye(length(y2))*(sigma_n^2)); % adding noise
   invKX2X2 = linsolve(KX2X2,eye(size(KX2X2)),opts); %computing inverse
   clear KX2X2; %freeing memory
   KXTX2 = (sigma_f^2)*computeCorrelMat(xtest,x2,theta); %computing covariance matrix for test points and the second dataset
   mu2 = beta + (KXTX2*(invKX2X2*(y2-beta))); %computing mean vector conditioned on the second dataset 
   K2 = KXTX2*invKX2X2 ; 
   K = K + (K2*(KXTX2')); %adding posterior covariance of second function to the first one
   clear KXTX2; %freeing memory

   KX2X1 = (sigma_f^2)*computeCorrelMat(x2,x1,theta); %computing cross-covariance between x2 and x1
   K = K - (2*K2*KX2X1*K1) ; %subtracting twice the posterior cross-covariance from the posterior covariance
   clear KX2X1; %freeing memory
   out.mu = mu2 - mu1 ; %computing the difference of the mean vectors
   out.K = (K + K')/2; %ensuring that diff matrix stays symmetric, numerically
end