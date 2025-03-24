function [ci] = bonnetwright(rs, alpha, n)
if size(rs,2) < size(rs,1)
    rs = rs';
end

b      = 3;
crit_z = norminv(1-(alpha/2));
L1 = 0.5*(log(1+rs) - log(1-rs)) - (((1+rs.^2/2).^(1/2) * crit_z) / ((n-b)^(1/2)));
L2 = 0.5*(log(1+rs) - log(1-rs)) + (((1+rs.^2/2).^(1/2) * crit_z) / ((n-b)^(1/2)));
cilow  = (exp(2*L1)-1) ./ (exp(2*L1)+1);
cihigh = (exp(2*L2)-1) ./ (exp(2*L2)+1);
ci = [cilow; cihigh]';
end
