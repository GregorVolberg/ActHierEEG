function [vec] = get_permvector(length_in)
num1 = randi(length_in);
tmpvec = [ones(num1, 1)*100; zeros(length_in - num1, 1)]; % 0 oder 100
vec = tmpvec(randperm(length_in)); % shuffle
end
