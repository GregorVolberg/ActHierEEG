function [rc_indices] = get_rc_indices(n_conditions)
mat = tril(ones(n_conditions, n_conditions), -1);
ind = find(mat == 1);
[r, c] = ind2sub(size(mat), ind);
rc_indices = [r, c]; 
end