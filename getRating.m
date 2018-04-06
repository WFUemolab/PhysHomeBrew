function getRating(filestr)

d = dir(['raw/' filestr]);

for i = 1:length(d)
    [a,names{i},c] = fileparts(d(i).name);
end

cd rating;

for i = 1:length(names)
    [tmp(:,1), tmp(:,2)] = textread([names{i} '.rating.txt'], '%d%f', 'headerlines', 1);
    rate(1,i*2-1:i*2) = repmat(str2num(names{i}(4:8)),1,2);
    rate(2:size(tmp,1)+1, i*2-1:i*2) = tmp;
    clear tmp;
end

save rate rate
csvwrite('rateDataAll.txt', rate);
cd ..
