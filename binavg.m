function x = binavg(dat, binsize)

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
