set xlabel 'Temperature'
set ylabel 'Time'
set zlabel 'MSD'
a = 0.1
f(x,y) = a*x * y**2
fit f(x,y) 'tmp/avg-msd.dat' using 1:2:3:(1) via a
splot 'tmp/avg-msd.dat' using 1:2:3, f(x,y)
print 'f(x,y)=',a,' * x * y**2'
