import gnuplot, os
import genopt
import genopt/utils
import genopt/genoptions

options.maxiter = 2000

var
    x,y,z: seq[float]
    p: Population

var 
    target = zeros(10)
        
p = genpopulation(target,20)
echo(p)

proc calc_fitness(s: seq[float]): float =
    var diff: seq[float]

    for i in 0..<len(s):
        diff.add(target[i] - s[i])
    for d in diff:
        result += d*d

fitnessfunction = calc_fitness

var
    res, t, conv: seq[float]

(res, t, conv) = iterate(p)

echo(res)

cmd "set size square"
cmd "set logscale xy"
cmd "set terminal png size 500,500"
cmd "set output 'convergence.png'"
plot t, conv