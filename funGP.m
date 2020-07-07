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
%   main function for funGP algorithm
%   INPUTS:
%       x1: Input variable matrix or vector for first set of points
%       y1: response vector for first set of points
%       x2: Input variable matrix or vector for second set of points   
%       y2: response vector for second set of points
%       xtest: Input variable matrix or vector for the test points at which the mean vector 
%              and covariance matrix of the difference has to be computed
%       confLevel: Significance level for computing the difference band
%   OUTPUT: A data structure with the following elements
%       differ: a boolean stating whether the functions statistically differ or not
%       nPoints: number of points on which the functions differ
%       muDiff: a vector of pointwise difference in the mean vector
%       band: a matrix of two columns, first column for lower bounds
%             and the second column for upper bounds for the confidence band. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[out] = funGP(x1, y1, x2, y2, xtest, confLevel, option)
    default_option = struct;
    default_option.truncation.rel_threshold = 1e-6;
    default_option.zsamples = 1000;
    if exist('option','var')
        option_names = fieldnames(option);
        if any(~ismember(option_names,{'truncation','zsamples'}))
            error('option can only contain two fields: truncation and zsamples');
        end
        if ismember('truncation',option_names)
            truncation_type = fieldnames(option.truncation);
            if length(truncation_type) > 1 || ~ismember(truncation_type,{'rel_threshold', 'abs_threshold', 'trunc_number'})
                error("truncation must be specified either in terms of 'rel_threshold', 'abs_threshold', or 'trunc_number'.");
            end
        else 
            option.truncation = default_option.truncation;
        end
        if ~ismember('zsamples',option_names)
            option.zsamples = default_option.zsamples;
        end
    else
        option = default_option;
    end

    estimatedParams = estimateHyperparameters(x1,y1,x2,y2); %estimate the hyperparameters
    diffCov = computeDiffCov(x1,y1,x2,y2,xtest,estimatedParams.params); %compute the difference in mean vector and covariance mat
    uband = computeConfBand(diffCov.K, confLevel, option); %compute the positive band  
    nPoints = length(find(abs(diffCov.mu) > uband)); %find the number of points that are beyond the band
    %set differ to true if number of points beyond band is non zero
    if nPoints > 0
        differ = true;
    else
        differ = false;
    end
    %prepare the return data structure
    out.differ = differ;
    out.nPoints = nPoints;
    out.muDiff = diffCov.mu;
    out.band = [-uband, uband];
end
