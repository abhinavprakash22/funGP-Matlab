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
%   Function to estimate the hyperparameters from two datasets under the null
%   hypothesis that the functions generating the datasets are the same.
%   INPUTS:
%       x1: Input variable matrix or vector for first set of points
%       y1: response vector for first set of points
%       x2: Input variable matrix or vector for second set of points   
%       y2: response vector for second set of points
%   OUTPUT: A data structure with following elements
%       llval: optimal loglikelihood value
%       params: a data structure with estimated hyperparameters
%       exitflag: an integer stating the optimization termination status     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[out] = estimateHyperparameters(x1,y1,x2,y2) 
    x = [x1; x2]; %combine the input values
    y = [y1; y2]; %combine the response
    ncov = size(x,2); %number of input variables
    beta = mean(y); %initial value of beta to feed into optimization
    sigma_f = std(y)/sqrt(2); %initial value of sigma_f
    sigma_n = sigma_f; %initial value of sigma_n
    theta = std(x); %initial value of lengthscales
    par0 = [theta, sigma_f, sigma_n, beta]; % combine the initial values of all the hyperparameters in par0
    %define a function handle for the objective
    obj = @(par)computeGPloglik(x,y,struct('theta',par(1:ncov),'sigma_f',par(ncov+1),'sigma_n',par(ncov+2),'beta',par(ncov+3))); 
    %define options for optimization
    options = optimoptions('fminunc','Display','none','SpecifyObjectiveGradient',true);
    %optimize
    [sol, llval, exitflag] =  fminunc(obj,par0,options);
    out.llval = -llval; %convert negative loglikelihood to loglikelihood
    out.params.theta = abs(sol(1:ncov)); %set absolute value of theta in params data structure
    out.params.sigma_f = abs(sol(ncov+1)); %set absolute value of sigma_f in params data structure
    out.params.sigma_n = abs(sol(ncov+2)); %set absolute value of sigma_n in params data structure
    out.params.beta = sol(ncov+3); %set value of beta in params data structure   
    out.exitflag = exitflag;
end