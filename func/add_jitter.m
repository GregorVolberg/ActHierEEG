function [jittered] = add_jitter(invec, jitterwidth)
rng(44);
jitter = -jitterwidth + (jitterwidth-(-jitterwidth)).*rand(length(invec),1);
jittered = invec + jitter;
end
