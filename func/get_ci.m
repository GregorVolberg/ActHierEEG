function [ci] = get_ci(arrayIn, rndIn) 
r2z = @(x) atanh(x); % r to fisher's z
z2r = @(x) tanh(x);  % fisher's z to r

rndcells  = arrayfun(@(x) z2r(shiftdim(squeeze(mean(r2z(arrayIn(:,rndIn{x},:)), 2))', 1)), ...
                       1:numel(rndIn), 'UniformOutput', false);
rndmat = reshape(cell2mat(rndcells), [size(rndcells{1},1), size(rndcells{1},2), numel(rndIn)]);

low  = arrayfun(@(x) squeeze(quantile(rndmat(x,:,:), 0.025, 3)), ...
                      1:3, 'UniformOutput', false);
high = arrayfun(@(x) squeeze(quantile(rndmat(x,:,:), 0.975, 3)), ...
                      1:3, 'UniformOutput', false);
                  
 ci(:,:,1) = reshape(cell2mat(low), [size(arrayIn, 3), size(arrayIn,1)]);
 ci(:,:,2) = reshape(cell2mat(high), [size(arrayIn, 3), size(arrayIn,1)]);
 end

