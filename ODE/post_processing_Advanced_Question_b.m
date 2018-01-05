% plotting the total energy of the system
kinA = 0.5*mA*sum(vA.^2, 2);
kinB = 0.5*mB*sum(vB.^2, 2);
tot = kinA + kinB + pot;

figure(1)
plot(t, pot)
hold on
plot(t, kinA)
plot(t, kinB)
plot(t, tot)
xlabel('t')
ylabel('E')
legend('pot','kinA','kinB','tot')
title('Energy vs time');
grid on;

% plotting position
figure(2)
plot(t, xA(:,1))
hold on
plot(t, xB(:,1))
xlabel('t')
ylabel('x(t)')
legend('A','B')
title('position for A and B in x-dirn');
grid on;

figure(3)
plot(t, xA(:,2))
hold on
plot(t, xB(:,2))
xlabel('t')
ylabel('y(t)')
legend('A','B')
title('position for A and B in y-dirn');
grid on;

figure(4)
plot(t, xA(:,3))
hold on
plot(t, xB(:,3))
xlabel('t')
ylabel('z(t)')
legend('A','B')
title('position for A and B in z-dirn');
grid on;

figure(5)
com_vel = (mA*vA + mB*vB)/(mA + mB);
com_mom = mA*(vA - com_vel) + mB*(vB - com_vel);
plot(t, com_mom(:,1));
hold on;
plot(t, com_mom(:,2));
plot(t, com_mom(:,3));
xlabel('t')
ylabel('mom')
legend('momx','momy','momz')
title('Linear momentum of com');
grid on;

figure(6)
com_vA = vA - com_vel;
com_vB = vB - com_vel;
com_kin = 0.5*(mA*sum(com_vA.^2, 2) + mB*sum(com_vB.^2, 2));
com_tot = com_kin + pot;
plot(t, com_kin)
hold on;
plot(t, pot)
plot(t, com_tot)
xlabel('t')
ylabel('E')
legend('kin','pot','tot')
title('energy of center of mass');
grid on;

figure(7)
mu = (mA*mB)/(mA + mB);
xAB = xB - xA;
vAB = vB - vA;
com_ang_mom = 2*mu^2*cross(xAB, vAB);
plot(t, com_ang_mom(:,1));
hold on;
plot(t, com_ang_mom(:,2));
plot(t, com_ang_mom(:,3));
xlabel('t')
ylabel('angmom')
legend('angmomx','angmomy','angmomz')
title('Angular momentum center of mass');
grid on;

% filewriting
csvwrite('./data/Adv_b_com_energy.dat',[t', com_kin, pot, com_tot]);
csvwrite('./data/Adv_b_com_lin_mom.dat',[t', com_mom]);
csvwrite('./data/Adv_b_com_ang_mom.dat',[t', com_ang_mom]);