% check negative reaction times & other weird stuff and discard trials
wtrl = cell(size(block));
for fl = 1:size(block,2)
    [row,col] = find(isnan(trlsmp{fl}));
    wtrl{fl} = [wtrl{fl};row,col];
    u = unique(wtrl{fl}(:,1));
    c = histc(wtrl{fl}(:,1),u);
    if find(c>1)
        f = find(wtrl{fl}(:,1) == u(c>1));
        wtrl{fl}(f(2),:) = [];
    end
    [~, I] = sort(wtrl{fl}(:,1));
    wtrl{fl} = wtrl{fl}(I,:);
end