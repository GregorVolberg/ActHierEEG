function [cbcolor] = get_cbPalette()
% color blind - friendly palette
cbPalette = {'#999999', '#E69F00', '#56B4E9', ...
               '009E73', '#F0E442', '#0072B2', ...
               '#D55E00', '#CC79A7'};

cbcolor = nan(numel(cbPalette), 3);

for j = 1:numel(cbPalette)
str = cbPalette{j};
cbcolor(j,:) = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
end

end
