function xrms = rmssd(x)
% function x = rmssd(x)
%
% Calculates the root-mean sucessive squared difference, a measure of the
% variability in a time-series

xrms = sqrt(sum((diff(x)).^2) ./ (length(x)-1));
