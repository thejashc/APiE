function [] = post_processing(x,v, t)

%Plotting energies, positions, velocities and phase space trajectories
global x0 v0 omg0 k m friction gamma forcing fmag omg opt fwrite 

if (friction == 0 && forcing == 0)
    
    % analytical solutions for position and velocity
    x_anal =  x0*cos(omg0*t)     + (v0/omg0)*sin(omg0*t);
    v_anal = -x0*omg0*sin(omg0*t) + v0*cos(omg0*t);
    
    % position
    figure(1)
    ind = 1:10:length(t);
    x_anal_handpick = x_anal(ind);
    hold on
    plot(t(ind),x(ind),'ko');
    plot(t,x_anal,'r.');
    plot(t(ind), x_anal_handpick, 'ko','MarkerFaceColor','k')
    %csvwrite('./data/q5_euler_dt_2e-2_position_vs_time.dat',[t', x]);
    
    % velocity
    figure(2)
    v_anal_handpick = v_anal(ind);
    hold on
    plot(t(ind),v(ind),'k^','MarkerSize',3.5);
    plot(t,v_anal,'r^','MarkerSize',3.5);
    plot(t(ind), v_anal_handpick, 'k^','MarkerSize',3.5,'MarkerFaceColor','k')
    
    % energy
    E_init = 0.5*(k*x0^2 + m*v0^2)*ones(length(t),1);
    pot  = 0.5 * k * x.^2;
    kin  = 0.5 * m * v.^2;
    tot = pot + kin;
    figure(3)
    hold on
    plot(t,pot);
    plot(t,kin);
    plot(t,tot);
    plot(t,E_init)
    legend('Potential','Kinetic','Total','Initial')
    %csvwrite('./data/q11_ode45_energy_vs_time.dat',[t, kin, pot, tot]);
    %csvwrite('./data/q11_verletLF_energy_vs_time.dat',[t', kin, pot, tot]);
    csvwrite('./data/q11_euler_energy_vs_time.dat',[t', kin, pot, tot]);

    % plotting the phase-space plot for the spring mass system
    figure(4)
    a = sqrt((v0/omg0)^2 + (x0)^2);
    b = sqrt(v0^2 + (x0*omg0)^2);
    syms theta
    x_ellipse = a*cos(theta);
    y_ellipse = b*sin(theta);
    scatter(x, v, 'b.')
    axis square
    hold on;
    fplot(x_ellipse,y_ellipse)
    %csvwrite('./data/q5_euler_dt_2e-2_phase_space.dat',[v, x]);
    
elseif (friction == 1 && forcing == 0)
    %% post processing with friction
    % analytical solutions for position and velocity
    alpha = gamma/(2*m);
    omg_new = sqrt(omg0^2 - alpha^2);
    x_anal = (x0*cos(omg_new*t)     + (v0/omg_new)*sin(omg_new*t)).*exp(-alpha*t);
    v_anal = (-x0*omg_new*sin(omg_new*t) + v0*cos(omg_new*t)).*exp(-alpha*t) -alpha*exp(-alpha*t).*(x0*cos(omg_new*t)+(v0/omg_new)*sin(omg_new*t));
    
    % position
    figure(1)
    plot(t, x_anal, 'r', 'LineWidth', 2)
    hold on;
    plot(t,x,'k.');
    legend('Analytical', 'Numerics')
    % writing
    if (opt == 2 && fwrite == 1)
        csvwrite('./data/q12_verletLF_position_vs_time.dat', [t', x]);
    elseif (opt == 4 && fwrite == 1)
        csvwrite('./data/q12_ode45_position_vs_time.dat', [t, x]);
    end
    
    % velocity
    figure(2)
    plot(t,v_anal, 'r', 'LineWidth', 2)
    hold on;
    plot(t,v,'k.');
    legend('Analytical', 'Numerics')
    
    % energy
    E_init = 0.5*(k*x0^2 + m*v0^2);
    pot  = 0.5 * k * x.^2;
    kin  = 0.5 * m * v.^2;
    tot = pot + kin;
    figure(3)
    plot(t, E_init*exp(-(gamma/m)*t))
    hold on;
    plot(t, tot)
    legend('Analytical', 'Numerics')
    % writing
    if (opt == 2 && fwrite == 1)
        csvwrite('./data/q12_verletLF_energy_vs_time.dat', [t', tot]);
    elseif (opt == 4 && fwrite == 1)
        csvwrite('./data/q12_ode45_energy_vs_time.dat', [t, tot]);
    end
    
elseif (friction == 0 && forcing == 1)
    %% post processing with forcing
    % analytical solutions for position and velocity
    tan_phi = (-v0/omg0)/(x0 - (fmag/m)*(1/(omg0^2 - omg^2)));
    phi = atan(tan_phi);
    B = -(v0/omg0)*(1/sin(phi));
    x_anal = (fmag/m)*(1/(omg0^2 - omg^2))*cos(omg*t) + B*cos(omg0*t + phi);
    
    % position
    figure(1)
    plot(t,x,'k.');
    hold on;
    plot(t,x_anal, 'r')
    
    % velocity
    figure(2)
    plot(t,v,'k.');
    hold on;
    
elseif (friction == 1 && forcing == 1)
    %% post processing with both friction and forcing
    % analytical solutions for position and velocity
    phi = atan((gamma/m)*omg/(omg0^2 - omg^2));
    Amp = (fmag/m)*(1/( (omg0^2 - omg^2)^2 + (omg*(gamma/m))^2 )^(1/2));
    x_anal = Amp*cos(omg*t - phi);
    
    % write data to a file
    %     fname = ['./data/A_f_omg_omg0_gamma_',num2str(gamma),'_m_',num2str(m),'.dat'];
    %     fileID = fopen(fname,'a');
    %     fmt = '%6f %6f %6f %6f\n';
    %     fprintf(fileID,fmt, [Amp fmag omg omg0]);
    %     fclose(fileID);
    
    %     % position
    figure(1)
    plot(t,x,'k.');
    hold on;
    plot(t,x_anal, 'r')
    % writing
    if (opt == 2 && fwrite == 1)
        csvwrite('./data/q15_verletLF_position_vs_time.dat', [t', x]);
    elseif (opt == 4 && fwrite == 1)
        csvwrite('./data/q15_ode45_position_vs_time.dat', [t, x]);
    end
    
    % velocity
    figure(2)
    plot(t,v,'k.');
    hold on;
end

end

