function [x, yp, yc] = areaWeightedHistogram(x_in,a_in,n)
%% Function to calculate area weighted histograms
% Inputs:
% x_in - Values of data being assed (such as inclination angle)
% a_in - Weighting (such as face surface area)
% n - number of bin edges
%Outputs can be used in a plot/bar chart ex: plot(cax,y,yp);
% x - Values of bin centres
% yp - Probability of a value falling within a bin
% yc - Cumulative probability distribution 
try
    binEdges = linspace(min(x_in),max(x_in),n);
catch
    binEdges = linspace(0,1,n);
end

x = (binEdges(2:end)+binEdges(1:end-1))/2;
totalA = sum(a_in,"all","omitnan");
yp = zeros(size(x)); yc = zeros(size(x));
for i = 1:length(x)
    inds = x_in>=binEdges(i)&x_in<binEdges(i+1);
    yp(i) = sum(a_in(inds==1),"all","omitnan")/totalA;
    yc(i) = sum(yp,"all","omitnan");
end
yp = yp.*n;
end