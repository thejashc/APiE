% Harmonic Oscillator
clc; clearvars;

global A x0 v0 omg0 k m gamma friction forcing fmag omg

factor = 0.01:0.001:5;
factor = 1;

for expt = 1:length(factor)
    % parameters
    prd = 5;           % number of periods to simulate
    n = 10000;          % number of intervals
    A = 1;              % Amplitude
    k = 5;              % Spring constant
    m = 2;              % mass
    omg0 = sqrt(k/m);   % natural frequency
    T = 2*pi*prd/omg0;  % Total time
    typT = 2*pi/omg0;
    friction = 0;
    gamma = 1e-4;               % friction coefficient
    forcing = 0;
    fmag = 1.0;                % magnitude of force
    omg = factor(expt)*omg0;   % frequency of forcing
    
    opt = 1;
    
    % initialize arrays & initial values
    t = linspace(0,T,n);
    DeltaT = t(2) - t(1);
    x = zeros(length(t),1);
    v = zeros(length(t),1);
    w = zeros(length(t),1);
    
    % values at t=0
    % initial position and velocity
    x0 = 0;
    v0 =  +A*omg0;
    v(1) = v0;
    x(1) = x0;
    
    % w(1) using analytical solution
    % w(1) = -x0*omg0*sin(omg0*(-DeltaT/2)) + v0*cos(omg0*(-DeltaT/2));
    % w(1) using Forward Euler
    w(1) = v0 - (DeltaT/2)*(-(k/m)*x(1));
    
    %% Forward-Euler
    if (opt == 1)
        for i = 1:(length(t) - 1)
            v(i+1) = v(i) + force(x(i))*(DeltaT/m);
            x(i+1) = x(i) + v(i)*DeltaT;
        end
        
        % plotting the results
        post_processing(x,v,t);
        
        %% Verlet-Leap frog
    elseif (opt == 2)
        for i = 1:(length(t) - 1)
            
            if (friction == 0 && forcing == 0)
                % calculate the spring force
                fc = force(x(i));
                
                % advance half-step velocity
                w(i+1) = w(i) + (fc/m)*DeltaT;
                
                % advance position
                x(i+1) = x(i) + w(i+1)*DeltaT;
                
                % velocity calculation
                %v(i) = w(i);
                v(i) = (w(i+1)+w(i))/2;
                
            elseif (friction == 1 && forcing == 0)
                %% implicit formula
                %             % calculate the spring force
                %             fc = force(x(i));
                %
                %             % advance half-step velocity
                             a1 = (1 - gamma*DeltaT/(2*m))/(1 + gamma*DeltaT/(2*m));
                             a2 = (DeltaT*fc/m)/(1 + gamma*DeltaT/(2*m));
                             w(i+1) = a1*w(i) + a2;
                %
                %             % advance position
                %             x(i+1) = x(i) + w(i+1)*DeltaT;
                %
                %             % velocity calculation
                %             v(i) = (w(i+1)+w(i))/2;
                %% iterative procedure
                % calculate the spring force
                fc = force(x(i));
                
                % advance half-velocity using previous half-velocity
                w(i+1) = w(i) + (DeltaT/m)*(fc - gamma*w(i));
                
                % advance position
                x(i+1) = x(i) + DeltaT*w(i+1);
                
                % calculate approximate velocity
                if (i == 1)
                    x_ref = x0 - (DeltaT/2)*v0;
                    v(i) = (x(i+1) - x_ref)/(2*DeltaT);
                else
                    v(i) = (x(i+1)-x(i-1))/(2*DeltaT);
                end
                
                % recalculate the half-velocity using full velocity
                w(i+1) = w(i) + (DeltaT/m)*(fc - gamma*v(i));
                
                % advance position
                x(i+1) = x(i) + DeltaT*w(i+1);
                
                % calculate approximate velocity
                v(i) = (w(i+1) + w(i))/2;
                
            elseif (friction == 0 && forcing == 1)
                fc = force(x(i));
            elseif (friction == 1 && forcing == 1)
                fc = force(x(i), v(i));
            end
            
        end
        % additional step
        % increment i
        i = i+1;
        
        % calculate force
        f = force(x(i));
        
        % advance half-step velocity
        w_fin = w(i) + (f/m)*DeltaT;
        
        % velocity calculation
        v(i) = (w_fin+w(i))/2;
        
        % plotting the results
        post_processing(x,v,t);
        
        %% verlet
    elseif (opt == 3)
        
        for i = 1:(length(t) - 1)
            if (i == 1)
                x(i+1) = x(i) + v(i)*DeltaT;
                
                v(i+1) = (x(i+1) - x(i))/(DeltaT);
            else
                x(i+1) = 2*x(i) -x(i-1) + force(x(i))*(DeltaT^2)/m;
                
                % calculating the velocity
                v(i) = (x(i+1) - x(i-1))/(2*DeltaT);
            end
        end
        
        x_fin = 2*x(i+1) - x(i) + force(x(i+1))*(DeltaT^2)/m;
        
        % calculating the velocity
        v(i+1) = (x_fin - x(i))/(2*DeltaT);
        
        % plotting the results
        post_processing(x,v,t);
        
        %% Runge-Kutta using ode45 function
    elseif (opt == 4)
        range = linspace(0,T,n);
        %range = [0 T];
        z0 = [x0; v0];
        
        [t,z] = ode45(@integrator_friction, range, z0);
        
        % position and velocity
        x = z(:,1);
        v = z(:,2);
        
        % plotting the results
        post_processing(x, v, t);
        
        
    end
end
