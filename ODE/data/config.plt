#!gnuplot

#set terminal postscript eps size 3.5,2.62 enhanced color \
#    font 'Helvetica,10' linewidth 1

#set terminal pngcairo size 1366,768 enhanced font 'Helvetica,12' linewidth 1

set terminal epslatex size 3.5,2.62 color colortext font 'Helvetica,6' linewidth 1
set key samplen 4.0 spacing 2.75 font ",4"

# analytical solution
m = 2.0
k = 5.0
A = 1.0
Gamma = 0.2
fmag = 1.0
omg0 = sqrt(k/m)
tp = 2*pi/omg0
x0 = 0
v0 = A*omg0
pos(x) = (v0/omg0)*sin(omg0*x) + x0*cos(omg0*x)
vel(x) = v0*cos(omg0*x) - x0*omg0*sin(omg0*x)
energy(x) = 0.5*(m*vel(x)**2 + k*pos(x)**2)

# friction
alpha = (Gamma/(2.0*m))
omg_fr = sqrt( omg0**2.0 - alpha**2)
tp_fr = 2*pi/omg_fr
pos_fr(x) = (x0*cos(omg_fr*x) + (v0/omg_fr)*sin(omg_fr*x))*exp(-alpha*x)
vel_fr(x) = (-x0*omg_fr*sin(omg_fr*x) + v0*cos(omg_fr*x))*exp(-alpha*x) - alpha*exp(-alpha*x)*(pos_fr(x))
energy_fr(x) = 0.5*(m*vel_fr(x)**2 + k*pos_fr(x)**2)
tau = (1/(2.0*alpha))

# friction + forcing
omg = 0.5*omg0 
phi = atan((Gamma/m)*omg/(omg0**2 - omg**2))
Amp = (fmag/m)*(1/( (omg0**2 - omg**2)**2 + (omg*(Gamma/m))**2 )**0.5)
pos_fr_f(x) = Amp*cos(omg*x - phi)

b = sqrt((v0**2 + (x0*omg0)**2)/(omg0**2))
a = sqrt(v0**2 + (x0*omg0)**2)
set sample 1e4
