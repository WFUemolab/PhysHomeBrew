function x = winsor(x, varargin)
% function x = winsor(x, numstd)
%
% This function winsorizes the vector x and culls the lower and upper
% values > numstd and sets them equal to the value at numstd above and
% below the mean

if ~isempty(varargin)
    numstd = varargin{1};
else
    numstd = 3; %default = 3 sds
end

meanx = nanmean(x);
stdx = nanstd(x);

x(x>(meanx+numstd*stdx)) = meanx+numstd*stdx;
x(x<(meanx-numstd*stdx)) = meanx-numstd*stdx;
return
