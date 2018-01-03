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

# Question 15
load "config.plt"
set output "../plots/log_A_f_vs_log_omg.tex"

set xlabel "{/Symbol w}"
set ylabel "A/f"
set grid
set logscale x
set logscale y
set key bottom left 

# analytical solution
A_omg(gamma, x) = (1/m)*((omg0**2 - x**2)**2 + (x*(gamma/m))**2 )**(-0.5)

set sample 1e4

p 'A_f_omg_omg0_gamma_0.0001_m_2.dat' u ($3):($1/$2) every 20 w p pt 7 ps 0.5 t'{/Symbol g} = 10^{-4}',\
  'A_f_omg_omg0_gamma_0.5_m_2.dat' u ($3):($1/$2) every 20 w p pt 7 ps 0.5 t'{/Symbol g} = 0.5',\
  'A_f_omg_omg0_gamma_1_m_2.dat' u ($3):($1/$2) every 20 w p pt 7 ps 0.5 t'{/Symbol g} = 1.0',\
  'A_f_omg_omg0_gamma_2_m_2.dat' u ($3):($1/$2) every 20 w p pt 7 ps 0.5 t'{/Symbol g} = 2.0',\
  [0.001:8] A_omg(1e-4, x) w l ls 1 t'{/Symbol g} = 10^{-4}, Analytical',\
  [0.001:8] A_omg(0.5, x) w l ls 2 t'{/Symbol g} = 0.5, Analytical',\
  [0.001:8] A_omg(1.0, x) w l ls 3 t'{/Symbol g} = 1.0, Analytical',\
  [0.001:8] A_omg(2.0, x) w l ls 4 t'{/Symbol g} = 2.0, Analytical'

reset
