#plot "Z.out" w l
#pause -1
#reset

set term pdf color solid
set output "A.pdf"
set xr [0:1500]
set xl "Wavelnumber (1/cm)"
set yl "Power"
set grid x
plot\
   "FT.c111_4x2_0809121414.vel" using 1:4 title "Whole Slab" w l, \
   "FT.c111_4x2_0809121414_bulktop.vel" using 1:4 title "Bulk" w l, \
   "FT.c111_4x2_0809121414_toplayerXYZ.vel" using 1:4 title "Top Layer" w l

set output "B.pdf"

plot \
   "FT.c111_4x2_0809121414_toplayerXYZ.vel" using 1:4 title "Top Layer" w l, \
   "FT.c111_4x2_0809121414_toplayerX.vel" using 1:4 title "Top layer, [-1 1 0]" w l, \
   "FT.c111_4x2_0809121414_toplayerY.vel" using 1:4 title "Top layer, [1 1 -2]" w l, \
   "FT.c111_4x2_0809121414_toplayerZ.vel" using 1:4 title "Top layer, [1 1 1]" w l

