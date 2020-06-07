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
%   An example file for implementing funGP algorithm for the following two cases:
%   1) function f1(x) = sin(2*pi*x) and f2(x) = 0.8*sin(2*pi*x)
%   2) both f1(x) and f2(x) are same as sin(2*pi*x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Example 1: when f1(x) and f2(x) are different

rng(1); %setting random seed

x1 = rand(100,1); %random input points for dataset1
y1 = sin(2*pi*x1) + normrnd(0,0.1,100,1); %generating response for dataset1 by adding noise to f1(x)

x2 = rand(100,1); %random input points for dataset2
y2 = 0.8*sin(2*pi*x2) + normrnd(0,0.1,100,1); %generating response for dataset2 by adding noise to f2(x)

%plotting the data
scatter(x1,y1);
hold on;
scatter(x2,y2);
hold off;

confLevel = 0.95; %setting confidence level
xtest = linspace(0,1,100)'; %regular 1-dim grid of 100 points

%running the algorithm
result = funGP(x1, y1, x2, y2, xtest, confLevel);

%printing out the results
if result.differ == true
    fprintf("The functions are statistically different.\n");
    fprintf("Number of points for which the functions are different are %d out of %d test points.\n",result.nPoints,size(xtest,1));
else
    fprintf("The functions are statistically the same.\n");
end

%plotting the difference along with confidence band
plot(xtest,result.muDiff,'LineWidth',2,'Color','blue');
hold on;
plot(xtest, result.band(:,1),'--','LineWidth',2,'Color','red');
plot(xtest, result.band(:,2),'--','LineWidth',2,'Color','red');
hold off;


%%Example 2: when f1(x) and f2(x) are the same

rng(2); %setting random seed

x1 = rand(100,1); %random input points for dataset1
y1 = sin(2*pi*x1) + normrnd(0,0.1,100,1); %generating response for dataset1 by adding noise to f1(x)

x2 = rand(100,1); %random input points for dataset2
y2 = sin(2*pi*x2) + normrnd(0,0.1,100,1); %generating response for dataset2 by adding noise to f2(x)

%plotting the data
scatter(x1,y1);
hold on;
scatter(x2,y2);
hold off;

confLevel = 0.95; %setting confidence level
xtest = linspace(0,1,100)'; %regular 1-dim grid of 100 points

%running the algorithm
result = funGP(x1, y1, x2, y2, xtest, confLevel);

%printing out the results
if result.differ == true
    fprintf("The functions are statistically different.\n");
    fprintf("Number of points for which the functions are different are %d out of %d test points.\n",result.nPoints,size(xtest,1));
else
    fprintf("The functions are statistically the same.\n");
end

%plotting the difference along with confidence band
plot(xtest,result.muDiff,'LineWidth',2,'Color','blue');
hold on;
plot(xtest, result.band(:,1),'--','LineWidth',2,'Color','red');
plot(xtest, result.band(:,2),'--','LineWidth',2,'Color','red');
hold off;


