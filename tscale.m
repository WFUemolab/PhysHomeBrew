function x = tscale(x)
%x = scale(x)
%
%Centers and scales column vectors to mean 50 and standard dev. = 10

x = x - repmat(nanmean(x),size(x,1),1);

x = x ./ repmat(nanstd(x)./10, size(x,1),1) + 50;
return