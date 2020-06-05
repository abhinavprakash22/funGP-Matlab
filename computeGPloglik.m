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
%   Function to compute loglikelihood of a GP and its gradient with constant mean 
%   and ARD squared exponential covariance function
%   INPUTS:
%       x: Input variable matrix or vector 
%       y: response vector 
%       params: a data structure with the hyperparameters
%   OUTPUT:
%       llval: loglikelihood value
%       grval: gradient vector 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[llval, grval] = computeGPloglik(x, y, params)
    correlMat = computeCorrelMat(x,x,params.theta); %compute correlation matrix 
    covMat = ((params.sigma_f^2)*correlMat) + (eye(length(y))*(params.sigma_n^2)); %convert to covariance mat with noise
    cholCovMat = chol(covMat); %cholesky decomp of cov mat
    clear covMat; %free the memory
    diagChol = diag(cholCovMat); %get diag elements of chol mat
    oneVec = ones(length(y),1); %define a vector of ones
    %compute loglik value
    llval = ((1/2)*(y-(params.beta*oneVec))'*(cholCovMat\(cholCovMat'\(y-(params.beta*oneVec))))) + (sum(log(abs(diagChol)))) + (length(y)*log(2*pi)/2);
    % compute gradient if nargout > 1
    if nargout > 1
        solOneVec = cholCovMat\(cholCovMat'\oneVec); %computing inv([K + sigma_n^2])*oneVec
        invMat = cholCovMat\(cholCovMat'\eye(length(y))); %computing inv([K + sigma_n^2])
        clear cholCovMat; %free the memory
        %intermediate steps for computing gradient
        alpha = invMat*(y-(params.beta*oneVec)); 
        diffMat = (alpha*alpha') - invMat ; 
        clear invMat;
        grad_size = length(params.theta)+3; %length of gradient vector
        grval = zeros(grad_size, 1); %intiliazing gradient vector
        %computing gradients
        for i = 1:length(params.theta)
            delThetaMat = (params.sigma_f^2)*((((x(:,i) - (x(:,i)')).^2)./(params.theta(i)^3)).*correlMat);
            grval(i) = (-1/2)*trace(diffMat*delThetaMat);            
        end
        clear delThetaMat;
        delSigma_fMat = (2*params.sigma_f)*correlMat;
        clear correlMat;
        grval(length(params.theta)+1) =  (-1/2)*trace(diffMat*delSigma_fMat);
        clear delSigma_fMat;
        grval(length(params.theta)+2) = (-1/2)*trace(2*params.sigma_n*diffMat);
        grval(length(params.theta)+3) = (1/2)*((2*params.beta*oneVec'*solOneVec)-(y'*solOneVec)-(oneVec'*(alpha+(params.beta*solOneVec))));
        clear cholCovMat alpha oneVec;
    end
end