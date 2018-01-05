#! gnuplot

# Question 5 
# (a) position vs time
load "config.plt"
set output "../plots/q5_euler_dt_2e-4_position_vs_time.tex"

set xlabel "$t/t_p$"
set ylabel "$x(t/t_p)$"
set grid
set key bottom right
unset title
set datafile separator ','

plot 'q5_euler_dt_2e-4_position_vs_time.dat' u 1:2 every 100 w p pt 7 ps 0.5 lc 'black' t'Euler, $\Delta t=2\times 10^{-4}$',\
      pos(x) w l lc 'red' t'Analytical solution'

reset

# (b) Energy  vs time
load "config.plt"
set output "../plots/q5_euler_dt_2e-4_energy_vs_time.tex"

set xlabel "$t/t_p$"
set ylabel "$E(t/t_p)$"
set grid
set key top right
unset title
set datafile separator ','

plot[][0:4] 'q5_euler_dt_2e-4_energy_vs_time.dat' u ($1/tp):($2) every 100 w lp pt 7 ps 0.5 lc 'red' t'$E_{kin}$',\
     'q5_euler_dt_2e-4_energy_vs_time.dat' u ($1/tp):($3) every 100 w lp pt 7 ps 0.5 lc 'black' t'$E_{pot}$',\
     'q5_euler_dt_2e-4_energy_vs_time.dat' u ($1/tp):($4) every 1000 w lp pt 7 ps 0.5 lc 'blue' t'$E_{tot}$',\
     energy(x) w l lc 'green' t'$E_{init}$'

reset
# (c) phase space
load "config.plt"
set output "../plots/q5_euler_dt_2e-4_phase_space.tex"

set xlabel "$v(t)$"
set ylabel "$x(t)$"
set grid
set key top right
set datafile separator ','

set parametric
set trange [0:2*pi]
# parametric equation for ellipse
fx(t) = a*sin(t)
fy(t) = b*cos(t) 

set sample 1e4

plot 'q5_euler_dt_2e-4_phase_space.dat' u 1:2 every 100 w p pt 7 ps 0.5 lc 'blue' t'phase space',\
     fx(t), fy(t) w l lc 'black' t'Analytical Solution'

reset

# Question 6
# Energy comparison
load "config.plt"

set output "../plots/q6_euler_compare_dt_energy_vs_time.tex"

set xlabel "$t/t_p$"
set ylabel "$E(t/t_p)$"
set grid
set key bottom right
unset title
set datafile separator ','

plot[:][2:3] 'q5_euler_dt_2e-3_energy_vs_time.dat' u ($1/tp):($4) every 100 w p pt 7 ps 0.5 lc 'red' t'$\Delta t=2\times 10^{-3}$',\
     'q5_euler_dt_2e-4_energy_vs_time.dat' u ($1/tp):($4) every 1000 w p pt 7 ps 0.5 lc 'blue' t'$\Delta t=2\times 10^{-4}$',\
     'q5_euler_dt_2e-5_energy_vs_time.dat' u ($1/tp):($4) every 10000 w p pt 7 ps 0.5 lc 'green' t'$\Delta t=2\times 10^{-5}$',\
     energy(x) w l lc 'black' t'Analytical solution'

reset

# Question 7
# Energy comparison
load "config.plt"

set output "../plots/q7_verletLF_compare_dt_energy_vs_time.tex"

set xlabel "$t/t_p$"
set ylabel "$E(t/t_p)$"
set grid
set key top right
unset title
set datafile separator ','

plot[:][2:3] 'q7_verletLF_dt_2e-3_energy_vs_time.dat' u ($1/tp):($4) every 100 w p pt 7 ps 0.5 lc 'red' t'$\Delta t=2\times 10^{-3}$',\
     'q7_verletLF_dt_2e-4_energy_vs_time.dat' u ($1/tp):($4) every 1000 w p pt 7 ps 0.5 lc 'blue' t'$\Delta t=2\times 10^{-4}$',\
     'q7_verletLF_dt_2e-5_energy_vs_time.dat' u ($1/tp):($4) every 10000 w p pt 7 ps 0.5 lc 'green' t'$\Delta t=2\times 10^{-5}$',\
     energy(x) w l lc 'black' t'Analytical solution'

reset

# Question 8
# Energy comparison
load "config.plt"

set output "../plots/q8_verletLF_compare_dt_energy_vs_time_intermediate_velocities.tex"

set xlabel "$t/t_p$"
set ylabel "$E(t/t_p)$"
set grid
set key top right
unset title
set datafile separator ','

plot[:][2.49:2.51] 'q8_verletLF_dt_2e-3_energy_vs_time_intermediate_velocity_1.dat' u ($1/tp):($4) every 20 w lp pt 7 ps 0.5 lc 'red' t'$v[i\Delta t] \approx w[i]$',\
     'q8_verletLF_dt_2e-3_energy_vs_time_intermediate_velocity_2.dat' u ($1/tp):($4) every 100 w lp pt 7 ps 0.5 lc 'blue' t'$v[i\Delta t] \approx 0.5*( w[i+1] + w[i] )$',\
     energy(x) w l lc 'black' t'Analytical solution'

reset

# Question 10
load "config.plt"

set output "../plots/q10_ode45_timestep.tex"

set xlabel '$t/t_p$'
set ylabel '$\Delta t$'
set grid
unset key
unset title
set logscale xy
set format "%g"
set datafile separator ','

plot 'q10_ode45_timestep.dat' u 1:2  w lp pt 7 ps 0.5 lc 'red' t'$v[i\Delta t] \approx w[i]$',\

reset

# Question 11
# Energy comparison different schemes
load "config.plt"

set output "../plots/q11_compare_schemes_energy_vs_time.tex"

set xlabel "$t/t_p$"
set ylabel "$E(t/t_p)$"
set grid
set key top right
unset title
set datafile separator ','

plot 'q11_ode45_energy_vs_time.dat' u ($1/tp):($4)  w l lc 'red' t'RK4',\
     'q11_verletLF_energy_vs_time.dat' u ($1/tp):($4)  w l lc 'blue' t'Verlet leap-frog',\
     'q11_euler_energy_vs_time.dat' u ($1/tp):($4)  w l lc 'green' t'Euler',\
     energy(x) w l lc 'black' t'Analytical solution'

reset

# Question 12
# scheme comparison for case of friction -- position
load "config.plt"

set output "../plots/q12_ode45_vs_verletLF_position_vs_time.tex"

set xlabel "$t$"
set ylabel "$x(t)$"
set grid
set key top right
unset title
set datafile separator ','

plot 'q12_ode45_position_vs_time.dat' u 1:2  every 100 w p lc 'red' t'RK4',\
     'q12_verletLF_position_vs_time.dat' u 1:2 every 100 w p lc 'blue' t'Verlet leap-frog',\
     pos_fr(x) w l lc 'black' t'Analytical solution'

reset

# scheme comparison for case of friction -- position
load "config.plt"

set output "../plots/q12_ode45_vs_verletLF_energy_vs_time.tex"

set xlabel "$t$"
set ylabel "$E_{tot}(t)$"
set grid
set key top right
unset title
set datafile separator ','

plot 'q12_ode45_energy_vs_time.dat' u 1:2  every 200 w p lc 'red' t'RK4',\
     'q12_verletLF_energy_vs_time.dat' u 1:2 every 200 w p lc 'blue' t'Verlet leap-frog',\
     energy_fr(x) w l lc 'black' t'Analytical solution'

reset

# Question 13
# energy decay interpretation
load "config.plt"

set output "../plots/q13_ode45_vs_verletLF_energy_decay.tex"
set xlabel "$t$"
set ylabel "$E_{tot}(t)$"
set grid
set key top right
unset title
set datafile separator ','

plot 'q12_ode45_energy_vs_time.dat' u 1:2  every 200 w p lc 'red' t'RK4',\
     energy_fr(x) w l lc 'black' lw 2 t'Analytical solution',\
     energy_fr(0)*exp(-x/tau) w l lc 'violet' lw 2 t'$\exp(-t/\tau)$'

reset

# Question ecay interpretation
# position vs time for the forcing + friction case
load "config.plt"

set output "../plots/q15_verletLF_position_vs_time.tex"
set xlabel "$t$"
set ylabel "$x(t)$"
set grid
set key top right
unset title
set datafile separator ','

set arrow from 40.0,-1.5 to 40.0,1.5 nohead lc 'blue' linetype -10 dashtype 1 lw 3.5

plot 'q15_verletLF_position_vs_time.dat' u 1:2  every 1 w l linetype 10 dashtype 1 lw 3.5 lc 'red' t'Verlet leap-frog',\
     pos_fr_f(x) w l lc 'black' lw 2 t'asymptotic solution (long time)',\
     fmag w l lc 'green' linetype 10 dashtype 3 lw 3.5 notitle,\
     -1.0*fmag w l lc 'green' linetype 10 dashtype 3 lw 3.5 notitle
   
reset
# position vs time for the forcing + friction case -- transient
load "config.plt"

set output "../plots/q15_verletLF_position_vs_time_transient.tex"
set xlabel "$t$"
set ylabel "$x(t)$"
set grid
set key top right
unset title
set datafile separator ','


plot[0:40] 'q15_verletLF_position_vs_time.dat' u 1:2  every 1 w l linetype 10 dashtype 1 lw 3.5 lc 'red' t'Verlet leap-frog',\
     pos_fr_f(x) w l lc 'black' lw 2 t'asymptotic solution (long time)',\
     fmag w l lc 'green' linetype 10 dashtype 3 lw 3.5 notitle,\
     -1.0*fmag w l lc 'green' linetype 10 dashtype 3 lw 3.5 notitle
   
reset

# Contour for varying gamma
load "config.plt"

set output "../plots/f_omg_contour_gamma_vary.tex"

set samples 2500,2500
set isosamples 2, 7

#set pm3d
#set hidden3d

set xlabel '$\omega/\omega_{0}$'
set ylabel '$F$'
set zlabel '$A$'
set key top left
set grid

set xrange [1e-3:2.*omg0]
set yrange [0:fmag]
splot[][][0:] (y/m)*(1./( (omg0**2.0 - x**2.0)**2.0 + (x*(0.5/m))**2.0 )**0.5) w lines lc rgb '#f89441' lw 2.5 t'$\gamma/m\omega_{0} = 0.158$',\
              (y/m)*(1./( (omg0**2.0 - x**2.0)**2.0 + (x*(0.8/m))**2.0 )**0.5) w lines lc rgb '#a82296'  lw 2.5 t'$\gamma/m\omega_{0} = 0.252$'

reset

# maximum omg for a given F
load "config.plt"

set output "../plots/q15_max_omg_ratio_vs_F_gamma_compare.tex"
set ylabel '$\omega_{max}/\omega_{0}$'
set xlabel '$F$'
set grid
set key top right
unset title
set datafile separator ','


plot[][0.98:1.0] 'q15_max_omg_ratio_vs_F_gamma_0_5.dat' u 1:2  every 1 w l linetype 10 dashtype 1 lw 3.5 lc rgb '#f89441' t'$\gamma/m\omega_{0} = 0.158$',\
	   'q15_max_omg_ratio_vs_F_gamma_0_8.dat' u 1:2  every 1 w l linetype 10 dashtype 1 lw 3.5 lc rgb '#a82296' t'$\gamma/m\omega_{0} = 0.252$',\
   
reset
# Question Appendix A
load "config.plt"
set output "../plots/log_A_f_vs_log_omg.tex"

set xlabel '$\omega$'
set ylabel '$\frac{A}{F}$'
set grid
set logscale x
set logscale y
set key top right

# analytical solution
A_omg(gamma, x) = (1/m)*((omg0**2 - x**2)**2 + (x*(gamma/m))**2 )**(-0.5)

set sample 1e4
set arrow from omg0,1e-3 to omg0,1e4 nohead lc 'blue' linetype -10 dashtype 1 lw 2.5
set label '$\omega = \omega_{0}$' at omg0+1.0,5e-3 center

p[0.1:10][] 'A_f_omg_omg0_gamma_0.0001_m_2.dat' u ($3):($1/$2) every 30 w p pt 7 ps 0.5 lc rgb '#f89441' t'$\gamma = 10^{-4}$',\
  'A_f_omg_omg0_gamma_0.5_m_2.dat' u ($3):($1/$2) every 30 w p pt 7 ps 0.5 lc rgb '#a82296' t'$\gamma = 0.5$',\
  'A_f_omg_omg0_gamma_1_m_2.dat' u ($3):($1/$2) every 30 w p pt 7 ps 0.5 lc rgb '#000004'  t'$\gamma = 1.0$',\
  'A_f_omg_omg0_gamma_2_m_2.dat' u ($3):($1/$2) every 30 w p pt 7 ps 0.5 lc 'red' t'$\gamma = 2.0$',\
  [0.001:8] A_omg(1e-4, x) w l ls 1 lc rgb '#f89441' lw 2.5 notitle,\
  [0.001:8] A_omg(0.5, x) w l ls 2 lc rgb '#a82296' lw 2.5 notitle,\
  [0.001:8] A_omg(1.0, x) w l ls 3 lc rgb '#000004' lw 2.5 notitle,\
  [0.001:8] A_omg(2.0, x) w l ls 4 lc 'red' lw 2.5 notitle

reset

# Question Appendix B
# energy in the center of mass frame
load "config.plt"

set output "../plots/Adv_b_com_energy.tex"
set ylabel '$E(t)$'
set xlabel '$t$'
set grid
set key top right
unset title
set datafile separator ','

plot[][0:0.5] 'Adv_b_com_energy.dat' u 1:2  w l lw 3.5 lc rgb '#f89441' t'$E_{kin, com}$',\
     'Adv_b_com_energy.dat' u 1:3  w l lw 3.5 lc rgb '#a82296' t'$E_{pot, com}$',\
     'Adv_b_com_energy.dat' u 1:4  w l lw 3.5 lc rgb '#000004' t'$E_{tot, com}$'

reset

# linear momentum in the center of mass frame
load "config.plt"

set output "../plots/Adv_b_com_lin_mom.tex"
set ylabel '$\bm{p}_{tot, com}(t)$'
set xlabel '$t$'
set grid
set key top right
unset title
set datafile separator ','

plot 'Adv_b_com_lin_mom.dat' u 1:2  w l lw 3.5 lc rgb '#f89441' t'$p_{x, com}$',\
     'Adv_b_com_lin_mom.dat' u 1:3  w l lw 3.5 lc rgb '#a82296' t'$p_{y, com}$',\
     'Adv_b_com_lin_mom.dat' u 1:4  w l lw 3.5 lc rgb '#000004' t'$p_{z, com}$'

reset
# angular momentum in the center of mass frame
load "config.plt"

set output "../plots/Adv_b_com_ang_mom.tex"
set ylabel '$\bm{L}_{tot, com}(t)$'
set xlabel '$t$'
set grid
set key top right
unset title
set datafile separator ','

plot 'Adv_b_com_ang_mom.dat' u 1:2  w l lw 3.5 lc rgb '#f89441' t'$L_{x, com}$',\
     'Adv_b_com_ang_mom.dat' u 1:3  w l lw 3.5 lc rgb '#a82296' t'$L_{y, com}$',\
     'Adv_b_com_ang_mom.dat' u 1:4  w l lw 3.5 lc rgb '#000004' t'$L_{z, com}$'
