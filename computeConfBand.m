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
%   Function to compute the confidence band on difference of two GPs over a finite grid.
%   INPUTS:
%       K: Covariance matrix of the difference process over a finite grid.
%       confLevel: Significance level
%   OUTPUT:
%       band: positive vector of upper bounds for the confidence band. 
%              lower band can be computed as negative of the upper band.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[band] = computeConfBand(K,confLevel)
   
    K = (K + K')/2; %to make sure that the matrix is numerically symmetric
    [V,D] = eig(K); %eigen decomposition
    [d, sorted_index] = sort(diag(D),'descend'); %sort eigenvalues in decreasing order
    Ds = D(sorted_index,sorted_index); %sorted eigenvalue diagonal matrix
    Vs = V(:,sorted_index); %sorted eigenvectors
    %print largest and smallest eigenvalues
    %fprintf('Largest eigenvalue: %8.10f \n',max(d));
    %fprintf('Smallest eigenvalue: %8.10f \n',min(d));

    threshold = 0.001; %threshold to compute truncation number

    %find the number of eigenvalues greater than threshold*max(d)
    if (threshold*max(d)) > min(d)
        m = find(d<threshold*max(d),1);
    else
        m = size(K,1);
    end

    %print the tructation number and the smallest eigenvalue selected
    %fprintf('Number of eigenvalues selected: %d \n',m);
    %fprintf('Smallest eigenvalue selected: %f \n',d(m));


    R = sqrt(chi2inv(confLevel,m)); %compute the radius of coverage region with probability confLevel
    t_samples = 1000; %number of samples to be used in simulating the band at each point
    Z = zeros(t_samples,m); %matrix to strore the simulated standard normal random variables. 
    
    %sample random vectors till the number of samples reach t_samples
    t = 1;
    while t <= t_samples
        Z(t,:) = normrnd(0,1,m,1);
        if norm(Z(t,:))<=R
         t = t+1;
        end
    end

    Gx = Vs(:,1:m)*sqrt(Ds(1:m,1:m))*(Z'); %compute the simulated sample paths 
    band =  max(abs(Gx),[],2); %take the absolute maximum as the upper band.
end