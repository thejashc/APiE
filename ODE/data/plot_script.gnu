#! gnuplot

set terminal postscript eps size 3.5,2.62 enhanced color \
    font 'Helvetica,10' linewidth 1

set output "log_A_f_vs_log_omg.eps"

set xlabel "{/Symbol w}"
set ylabel "A/f"
set grid
set logscale x
set logscale y
set key bottom left 

# analytical solution
m = 2.0
k = 5.0
omg0 = sqrt(k/m)
A_omg(gamma, x) = (1/m)*((omg0**2 - x**2)**2 + (x*(gamma/m))**2 )**(-0.5)

set sample 1e4

p 'A_f_omg_omg0_gamma_0.0001_m_2.dat' u ($3):($1/$2) every 20 w p pt 7 ps 0.25 t'{/Symbol g} = 10^{-4}',\
  'A_f_omg_omg0_gamma_0.5_m_2.dat' u ($3):($1/$2) every 20 w p pt 7 ps 0.25 t'{/Symbol g} = 0.5',\
  'A_f_omg_omg0_gamma_1_m_2.dat' u ($3):($1/$2) every 20 w p pt 7 ps 0.25 t'{/Symbol g} = 1.0',\
  'A_f_omg_omg0_gamma_2_m_2.dat' u ($3):($1/$2) every 20 w p pt 7 ps 0.25 t'{/Symbol g} = 2.0',\
  [0.001:8] A_omg(1e-4, x) w l ls 1 t'{/Symbol g} = 10^{-4}, Analytical',\
  [0.001:8] A_omg(0.5, x) w l ls 2 t'{/Symbol g} = 0.5, Analytical',\
  [0.001:8] A_omg(1.0, x) w l ls 3 t'{/Symbol g} = 1.0, Analytical',\
  [0.001:8] A_omg(2.0, x) w l ls 4 t'{/Symbol g} = 2.0, Analytical'

