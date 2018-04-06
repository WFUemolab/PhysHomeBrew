function x = binavg(dat, binsize)
% function x = binavg(dat, 60)
%
% function to average data in a vector into smaller bins - binsize is in
% the unit of the original data. So if you have 100s of data and want 20s
% bins, binsize = 20. If you have 100s of data and want 20 bins then
% binsize = 5.

if size(dat,1) > size(dat,2)
    dat = dat';
    flip = 1;
else
    flip = 0;
end

numbins = size(dat,2)./binsize;

for i = 1:floor(numbins)
    x(1,i) = nanmean(dat(1,i*binsize+1-binsize:i*binsize));
    
end

if numbins - floor(numbins) > 0, 
    i = i +1;
    finalbin = binsize*(numbins-floor(numbins)); 

    x(1,i) = nanmean(dat(1,i*finalbin+1-finalbin:i*finalbin));
    
end

if flip
    dat = dat';
end
