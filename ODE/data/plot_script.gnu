#! gnuplot

set terminal postscript eps size 3.5,2.62 enhanced color \
    font 'Helvetica,20' linewidth 2

set output "log_A_f_vs_log_omg.eps"

set xlabel "{/Symbol w}"
set ylabel "A/f"
set grid
#set logscale x
#set logscale y

p 'A_f_omg_omg0_gamma_0.5_m_2.dat' u ($3):($1/$2) w p t'{/Symbol g} = 0.5',\
  'A_f_omg_omg0_gamma_1_m_2.dat' u ($3):($1/$2) w p t'{/Symbol g} = 1.0',\
  'A_f_omg_omg0_gamma_2_m_2.dat' u ($3):($1/$2) w p t'{/Symbol g} = 2.0',\
