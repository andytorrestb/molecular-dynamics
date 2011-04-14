set xlabel 'Temperature'
set ylabel 'Energy'
set logscale x

a = 0.02
b = 1
#f(x) = 1/b*log(x/a)
#fit f(x) 'tmp/temp-response.dat' using 1:2 via b, a
f(x) = log(x/a)
fit f(x) 'tmp/temp-response.dat' using 1:2 via a

plot f(x), 'tmp/temp-response.dat' using 1:2
print 'E(T)=1/',b,'*log(T/',a,')'
print 'T(E)=',a,'*exp(',b,'*E)'
