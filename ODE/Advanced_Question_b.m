% two bodies connected by a spring
clearvars;
clc;

% parameters
prd = 5;          % number of periods to simulate
n = 10000;          % number of intervals
A = 1;              % Amplitude
k = 1;              % Spring constant
mA = 2;             % mass of ball A
mB = 2;             % mass of ball B
omg0 = sqrt(k/mA);   % natural frequency
T = 2*pi*prd/omg0;  % Total time
typT = 2*pi/omg0;
Leq = 1;

% initialize arrays & initial values
t = linspace(0,T,n);
DeltaT = t(2) - t(1);
xA = zeros(length(t),3);
vA = zeros(length(t),3);
xB = zeros(length(t),3);
vB = zeros(length(t),3);
pot = zeros(length(t),1);
Len = zeros(length(t),1);

% values at t=0
% initial 3D position and velocity of the balls connected to spring
x0A = [0 0 0];
x0B = [1 1 1];
%normal = (x0A - x0B)/L;
% xTan = [0 0 sqrt(6)];
% vecTan = xTan - x0B;
% tangent = vecTan/sqrt(sum(vecTan.^2));
v0A = [+0.1 +0.1 +0.1];
v0B = [-0.1 -0.1 -0.1];
% v0A = 0.1*tangent;
% v0B = -1*v0A;

% position and velocity for balls A and B
vA(1,:) = v0A;
xA(1,:) = x0A;

vB(1,:) = v0B;
xB(1,:) = x0B;

L = sqrt(sum(( xA(1, :) - xB(1,:) ).^2));
pot(1,:) = 0.5*k*(L - Leq)^2;
Len(1) = L-Leq;

% update position and force
for i = 1:(length(t) - 1)
    
    if (i == 1)
        % Euler step
        xA(i+1,:) = xA(i,:) + vA(i,:)*DeltaT;
        xB(i+1,:) = xB(i,:) + vB(i,:)*DeltaT;
        
        vA(i+1,:) = (xA(i+1,:) - xA(i,:))/(DeltaT);
        vB(i+1,:) = (xB(i+1,:) - xB(i,:))/(DeltaT);
    else
        % force calculation
        L = sqrt(sum(( xA(i, :) - xB(i,:) ).^2));
        Len(i) = L-Leq;
        pot(i,:) = 0.5*k*(L - Leq)^2;
        normal = (xA(i,:) - xB(i,:))/L;
        fAB = -k*(L - Leq)*normal;
        fBA = -1.0*fAB;
        
        % update position for A and B
        xA(i+1,:) = 2*xA(i,:) - xA(i-1,:) + fAB*(DeltaT^2)/mA;
        xB(i+1,:) = 2*xB(i,:) - xB(i-1,:) + fBA*(DeltaT^2)/mB;
        
        % calculating the velocity for A and B
        vA(i, :) = (xA(i+1,:) - xA(i-1,:))/(2*DeltaT);
        vB(i, :) = (xB(i+1,:) - xB(i-1,:))/(2*DeltaT);
    end
    
end

% for time length(t)
L = sqrt(sum(( xA(end, :) - xB(end,:) ).^2));
Len(end) = L-Leq;
pot(end,:) = 0.5*k*(L - Leq)^2;
normal = (xA(end,:) - xB(end,:))/L;
fAB = -k*(L - Leq)*normal;
fBA = -1.0*fAB;

% update position for A and B
xA_fin = 2*xA(end,:) - xA(end-1,:) + fAB*(DeltaT^2)/mA;
xB_fin = 2*xB(end,:) - xB(end-1,:) + fBA*(DeltaT^2)/mB;

% calculating the velocity for A and B
vA(end, :) = (xA_fin - xA(end-1,:))/(2*DeltaT);
vB(end, :) = (xB_fin - xB(end-1,:))/(2*DeltaT);

% plotting results
post_processing_Advanced_Question_b;

%movie_Advanced_Question_b