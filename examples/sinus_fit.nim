import gnuplot, os
import random
import math
import "../genopt"
import "../genopt/utils"
import "../genopt/genoptions"

options.starting_sigma = 2.0
options.maxiter = 50000

var
    x,y, noisy, fit: seq[float]
    n = 200
    p: Population


x = linspace(0,2*Pi,n)
for i in 0..<n:
    y.add(sin(2*Pi*x[i])) 
    noisy.add(y[i] + (rand(2.0)-1.0))


var 
    # sinus paramters [Amp, Freq, phase, c]
    initial = @[0.5,0.5,0.5,0.5]


        
p = genpopulation(initial,20)
echo(p)

proc calc_fitness(s: seq[float]): float =
    for i in 0..<len(noisy):
        result += (noisy[i] - (s[0]*sin(2*Pi*x[i]*s[1]+s[2]) + s[3]) )^2
    result /= float(n)

fitnessfunction = calc_fitness

var
    res, t, conv: seq[float]

(res, t, conv) = iterate(p)

echo(res)

for i in 0..<n:
    fit.add((res[0]*sin(2*Pi*x[i]*res[1]+res[2]) + res[3]))

cmd "set size square"
cmd "set terminal png size 1000,1000"
#
plot x, y, "original"

plot x, noisy, "noisy"

plot x, fit, "fit"
cmd "set output 'sine.png'"
plot x, fit # workaround to get plotting to work
