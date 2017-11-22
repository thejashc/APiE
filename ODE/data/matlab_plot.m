clearvars;
clc;

X = importdata('A_f_omg_omg0_gamma_1_m_2.dat');
Y = importdata('A_f_omg_omg0_gamma_2_m_2.dat');
Z = importdata('A_f_omg_omg0_gamma_0.5_m_2.dat');

% analytical
syms omg
m = 2;
k = 5;
omg0 = sqrt(k/m);
A_omg = @(gamma, omg) (1/m)*((omg0^2 - omg.^2).^2 + ((gamma/m)*omg).^2).^(-0.5);

ind = 1:50:length(X(:,1));

plot(X(:,3),X(:,1)./X(:,2),'r.')
hold on;
plot(X(ind,3),A_omg(1, X(ind,3)), 'ro')
plot(Y(:,3),Y(:,1)./Y(:,2),'b.')
plot(Y(ind,3),A_omg(2, Y(ind,3)), 'bo')
plot(Z(:,3),Z(:,1)./Z(:,2),'g.')
plot(Z(ind,3), A_omg(0.5, Z(ind,3)), 'go')

