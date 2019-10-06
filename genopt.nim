## Module for genetic/evolutionary optimization
##
##
## Examples
## ========
##
## Simple Example:
##
##
## Optimize all numbers to zero
## ---------------------------
##
## .. code-block:: nim
##
##  import gnuplot, os
##
##  var
##      x,y,z: seq[float]
##      p: Population
##
##  x = linspace(1,2,10)
##  y = linspace(2,3,10)
##
##  mutate(x,1.0)
##  echo(x)
##
##  var 
##      target = zeros(10)   
##
##  p = genpopulation(target,20)
##  echo(p)
##
##  proc calc_fitness(s: seq[float]): float =
##      var diff: seq[float]
##
##      for i in 0..<len(s):
##          diff.add(target[i] - s[i])
##      for d in diff:
##          result += d*d
##
##  fitnessfunction = calc_fitness
##
##  var
##      res, t, conv: seq[float]
##
##  (res, t, conv) = iterate(p)
##
##    echo(res)
##
##  cmd "set size square"
##  cmd "set logscale xy"
##  cmd "set terminal png size 500,500"
##  cmd "set output 'convergence.png'"
##  plot t, conv
##


import math
import strformat
import random
import times
import stats
import genopt/utils


const crossover_size = 0.7
const mutation_rate = 0.5
const report_every = 50
const max_iter = 5000
const starting_sigma = 4.0

const max_val = 10.0
const min_val = -10.0


type 
    FitnessFunction* = proc(s: seq[float]): float

    Population* = seq[seq[float]]


var
    fitnessfunction*: proc(s: seq[float]): float
    handleunfit*: proc(unfit: seq[float]): seq[float]


proc handleunfit_proc(unfit: seq[float]): seq[float] =
    return unfit

handleunfit = handleunfit_proc



proc recombine(seq1,seq2 : seq[float]) : seq[float] =
    var 
        k: int
        n_crossover: int
        alpha: float

    result = zeros(len(seq1))
    n_crossover = int( float(len(seq1)) * crossover_size)
    #n_crossover = len(seq1)
    for i in 0..<n_crossover:
        k = rand(len(seq1)-1)
        alpha = 0.75#rand(1.0)#0.5
        result[k] = alpha * seq1[k] + (1 - alpha) * seq2[k]
        #result2[k] = alpha * seq2[k] + (1 - alpha) * seq1[k]


proc mutate(s: var seq[float], sigma: float, mutation_rate = mutation_rate) =
    proc rev_sigmoidal(x: float): float =
        result = 0.8 / (exp((-x+0.5)*15.0) + 1.0) + 0.1

    var 
        mutation = 0.0

    for i in 0..<len(s) :
        if rand(1.0) < mutation_rate:
            #mutation = np.random.normal() * sigma
            #if mutation > sigma * 1.0:
            #    mutation = sigma
            #if mutation < -sigma * 1.0:
            #    mutation = -sigma
            #s[i] += mutation

            mutation = (rand(1.0) - 0.5)*2.0 * sigma
            s[i] += mutation

            #if mutation > 0.0:
            #    s[i] += rev_sigmoidal(mutation)
            #else:
            #    s[i] += -rev_sigmoidal(-mutation)

            #mutation = 1.0 - (rand(0.3) - 0.15) * sigma
            #s[i] *= mutation


proc recombine_population(population: var Population) =
        # n_recombination = int(population.shape[1]/3)
        # n_recombination = int(population.shape[1])
    var
        n,k,l: int
        i_top,i_med: int
        child: seq[float]

    n = len(population)
    #for i in 0..4:
    #    population[^(i+1)] = child        

    for i in countdown(n-1,4):
        i_top = rand(4)
        #i_med = rand(int(  n / 2 ))
        i_med = rand(n-1)
        #echo(itop, " ", i_med)
        #echo((2 * i+1), " ",(2 * i+2))
        child = recombine(population[i_top], population[i_med])
        #child = recombine(population[k], population[l])
        population[i] = child


proc mutate_population(population: var Population, sigma: float)=
    # for i in prange(population.shape[1]):
    #    population[:, i] = mutate(population[:, i], sigma)

    for i in 4..<len(population):
        # if i < int(population.shape[1]/3):
            if i < 6:
                mutate(population[i], sigma / 10)  #
            elif i < 10:
                mutate(population[i], sigma / 2)  #
            else:
                mutate(population[i], sigma)  #


proc check_limits(population: var Population)=
    for j in 0..<len(population):
        for i in 0..<len(population[0]):
            if population[j][i] < min_val:
                population[j][i] = min_val
            if population[j][i] > max_val:
                population[j][i] = max_val    


proc calc_fitness(population: Population): seq[float] =
    result = zeros(len(population))
    for i in 0..<len(population):
        result[i] = fitnessfunction(population[i])


proc genpopulation*(initial: seq[float],n: int, sigma=(max_val-min_val)/1.5): Population =
    for i in 0..<n:
        result.add(initial)
        mutate(result[i],sigma,1.0)
    check_limits(result)
    #echo(result)


proc iterate*(p: var Population): (seq[float],seq[float],seq[float]) =
    var
        logpoints: seq[int]
        checkpoints: seq[int]
        convergence: seq[float]
        t: seq[float]
        fitness: seq[float]
        mean_fitness: float
        sorted_ind: seq[int]
        linreg : RunningRegress
        indices: seq[int]
        slope:float
        correlation: float

    # At which points to report on the progress?
    logpoints = arange(report_every, max_iter, report_every)
    checkpoints = arange(50, max_iter, 50)

    # List containing the convergence data for each iteration
    convergence = zeros(max_iter)

    # List containing the required time for each iteration
    t = zeros(max_iter)

    var
        sigma = starting_sigma
    
    let start_time = now()


    for i in 0..<max_iter:

        mutate_population(p, sigma)
        check_limits(p)

        fitness = calc_fitness(p)
        
        sorted_ind = argsort(fitness)

        fitness = fitness[sorted_ind]
        p = p[sorted_ind]

        for i in 0..<len(p):
            if classify(fitness[i]) == fcNan:
                #p[i] = handleunfit(p[i]) 
                #fitness[i] = 100
                p[i] = recombine(p[rand(4)],p[rand(4)])
                fitness[i] = 100

        recombine_population(p)

        mean_fitness = 0
        for i in 0..<len(fitness):
            mean_fitness += fitness[i]
        mean_fitness /= float(len(fitness))

        convergence[i] = fitness[0]

        t[i] = (float( inMicroseconds(now() - start_time) )/1e6)

        if fitness[0] < 0.0000000001:
            break

        if i < 500:
            sigma = starting_sigma   # 1
        #elif i < 1000:
        #    sigma = starting_sigma / 2 # 0.5
        else:
        #if i>100:
            if i in checkpoints:
                indices = arange(i - 100, i, 1)
                linreg.push(t[indices], convergence[indices])
                slope = linreg.slope()
                correlation = abs(linreg.correlation())
                linreg.clear()
                if slope > 0 and correlation > 0.1:
                    if sigma > starting_sigma * 0.000001:
                        sigma *= 0.95

                if slope > 0 and correlation < 0.01:
                    if sigma < starting_sigma:
                        sigma *= 1.05
 
       
        if i in logpoints:
            echo(fmt"i: {i:6d}, best: {fitness[0]:4.6f}, sigma: {sigma:4.6f}, corr: {correlation:1.6f}, slope: {slope:1.6f},")

    result = (p[0], t, convergence)


when isMainModule:

    import gnuplot, os

    var
        x,y,z: seq[float]
        p: Population

    x = linspace(1,2,10)
    y = linspace(2,3,10)

    mutate(x,1.0)
    echo(x)

    var 
        target = zeros(10)
        
    #for i in 0..<len(target):
    #    target[i] = 5.0

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
