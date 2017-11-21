function [f] = force(x)

global k

% spring force
f = -k*x;

end

