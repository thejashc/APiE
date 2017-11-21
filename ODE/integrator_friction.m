function dzdt = integrator_friction(t, z)
% The rule for updating variables

global m gamma friction forcing omg0 omg fmag

if (friction == 0 && forcing == 0)
    % no friction or forcing
    Amat = [0, 1;-omg0^2, 0];
    dzdt = Amat*z;
    
elseif (friction == 1 && forcing == 0)
    % friction only
    Amat = [0, 1;-omg0^2, (-gamma/m)];
    
    dzdt = Amat*z;
    
elseif (friction == 0 && forcing == 1)
    % forcing only
    Amat = [0, 1;-omg0^2, 0];
    Bmat = [0; (fmag/m)*cos(omg*t)];
    
    dzdt = Amat*z + Bmat;
    
elseif (friction == 1 && forcing == 1 )
    % friction and forcing
    Amat = [0, 1;-omg0^2, (-gamma/m)];
    Bmat = [0; (fmag/m)*cos(omg*t)];
    
    dzdt = Amat*z + Bmat;
end

