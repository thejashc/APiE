function dzdt = integrator(t, z)
% The rule for updating variables

global Amat

dzdt = Amat*z;

end

