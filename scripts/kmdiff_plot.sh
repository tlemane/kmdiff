#!/bin/bash

if ! command -v gnuplot &> /dev/null; then
  echo "'gnuplot' not found."
  exit
fi

function usage
{
  echo "usage: kmdiff_plots.sh <kmdiff_dir> <nb_controls> <nb_cases> <pc1> <pc2>"
  exit
}

if [[ "$1" == "help" || "$#" -ne 6 ]]; then
  usage
fi

filename=$1
nb_controls=$2
nb_cases=$3
pc1=$4
pc2=$5
output=$6

pc1_=$((pc1+1))
pc2_=$((pc2+1))
xlabel=$(head -n${pc1_} ${filename} | tail -n1)
ylabel=$(head -n${pc2_} ${filename} | tail -n1)

gnuplot_script="
set object rectangle from screen 0,0 to screen 1,1 behind fillcolor rgb 'black' fillstyle solid noborder;
set key textcolor rgb 'white';
set border lw 3 lc rgb 'white';
set grid lc rgb 'white';
set terminal svg size 1024,768;
set output '${output}';
set xlabel \"PCA ${pc1} - ${xlabel}%\" textcolor rgb 'white';
set ylabel \"PCA ${pc2} - ${ylabel}%\" textcolor rgb 'white';
plot '${filename}' every ::1::${nb_controls} using ${pc1}:${pc2} title \"Controls\" pt 6 ps 0.8 lc rgb 'green', \
     '${filename}' every ::${nb_controls} using ${pc1}:${pc2} title \"Cases\" pt 4 ps 0.8 lc rgb 'red';
"

gnuplot -e "${gnuplot_script}"

