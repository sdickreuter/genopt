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
import genopt/genoptions



type 
    ## Type for proc used to calculate the fitness of am individual of the population
    FitnessFunction* = proc(s: seq[float]): float

    ## Type which holds the the whole population with every individual being a seq[float]
    Population* = seq[seq[float]]


var
    ## global variable which holds the options for the algorithm
    options*: Genoptions

    ## global variable which points to the fitnessfunction to be used
    fitnessfunction*: proc(s: seq[float]): float
    
    #handleunfit*: proc(unfit: seq[float]): seq[float]


options = initGenoptions()


#proc handleunfit_proc(unfit: seq[float]): seq[float] =
#    return unfit

#handleunfit = handleunfit_proc



proc recombine(seq1,seq2 : seq[float]) : seq[float] =
    ## recombines two seqs and generates a "child" 
    ## seq1 and seq2 are combined such that properties from seq1 dominate in the child
    var 
        k: int
        n_crossover: int
        alpha: float

    result = zeros(len(seq1))
    n_crossover = int( float(len(seq1)) * options.crossover_size)
    for i in 0..<n_crossover:
        k = rand(len(seq1)-1)
        alpha = 0.75#rand(1.0)#0.5
        result[k] = alpha * seq1[k] + (1 - alpha) * seq2[k]
        #result2[k] = alpha * seq2[k] + (1 - alpha) * seq1[k]


proc mutate(s: var seq[float], sigma: float, mutation_rate = options.mutation_rate) =
    ## mutate a seq, which means to add a random amount to each element of the seq
    ## sigma defines the amplitude of the random amount while mutation_rate 
    ##defines how many elements are mutated 
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


proc recombine_population(population: var Population) =
    ## recombine individual of whole population
    ## fitter individuals don't get recombined and dominate the
    ## properties of the children
    var
        n,k,l: int
        i_top,i_med: int
        child: seq[float]

    n = len(population)   

    for i in countdown(n-1,4):
        i_top = rand(4)
        #i_med = rand(int(  n / 2 ))
        i_med = rand(n-1)
        child = recombine(population[i_top], population[i_med])
        #child = recombine(population[k], population[l])
        population[i] = child


proc mutate_population(population: var Population, sigma: float)=
    ## mutate population
    ## fitter individuals don't get mutated or mutated to a lesser extend (smaller sigma)
    for i in 4..<len(population):
        # if i < int(population.shape[1]/3):
            if i < 6:
                mutate(population[i], sigma / 10)  #
            elif i < 10:
                mutate(population[i], sigma / 2)  #
            else:
                mutate(population[i], sigma)  #


proc check_limits(population: var Population)=
    ## check limits of a population and set values accordingly
    for j in 0..<len(population):
        for i in 0..<len(population[0]):
            if population[j][i] < options.min_val:
                population[j][i] = options.min_val
            if population[j][i] > options.max_val:
                population[j][i] = options.max_val    


proc calc_fitness(population: Population): seq[float] =
    ## calculate fitness of every individual of a population using the globally
    ## defined "fitnessfunction"
    result = zeros(len(population))
    for i in 0..<len(population):
        result[i] = fitnessfunction(population[i])


proc genpopulation*(initial: seq[float],n: int, sigma=(options.max_val-options.min_val)/1.5): Population =
    ## generate new population
    ## takes a initial individual and mutates it to generate a broad spectrum of individuals
    for i in 0..<n:
        result.add(initial)
        mutate(result[i],sigma,1.0)
    check_limits(result)
    #echo(result)


proc iterate*(p: var Population): (seq[float],seq[float],seq[float]) =
    ## iterates the genetic optimization until maxiter is reached or the fitness reaches a threshold
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
    logpoints = arange(options.report_every, options.max_iter, options.report_every)
    
    # At which points to calculate slope and correlation of the convergence and adjust sigma
    checkpoints = arange(50, options.max_iter, 50)

    # List containing the convergence data for each iteration
    convergence = zeros(options.max_iter)

    # List containing the required time for each iteration
    t = zeros(options.max_iter)

    var
        sigma = options.starting_sigma
    
    let start_time = now()


    for i in 0..<options.max_iter:

        mutate_population(p, sigma)
        check_limits(p)

        fitness = calc_fitness(p)
        
        # calculate and save sorted indices 
        sorted_ind = argsort(fitness)

        # rearange fitness and population p by according to fitness
        fitness = fitness[sorted_ind]
        p = p[sorted_ind]

        # check for unfit (fitness == fcNan) individuals and replace them with new ones
        for i in 0..<len(p):
            if classify(fitness[i]) == fcNan:
                #p[i] = handleunfit(p[i]) 
                p[i] = recombine(p[rand(4)],p[rand(4)])
                fitness[i] = 100

        recombine_population(p)

        # calculate mean fitness
        #mean_fitness = 0
        #for i in 0..<len(fitness):
        #    mean_fitness += fitness[i]
        #mean_fitness /= float(len(fitness))

        convergence[i] = fitness[0]

        t[i] = (float( inMicroseconds(now() - start_time) )/1e6)
        
        # fitness threshold for stopping the iteration
        if fitness[0] < 0.0000000001:
            break

        if i < 500:
            sigma = options.starting_sigma   # 1
        #elif i < 1000:
        #    sigma = starting_sigma / 2 # 0.5
        else:
        #if i>100:
            if i in checkpoints:
                
                # calculate slope and correlation of the last 100 fitness-values
                indices = arange(i - 100, i, 1)
                linreg.push(t[indices], convergence[indices])
                slope = linreg.slope()
                correlation = abs(linreg.correlation())
                linreg.clear()
                
                # modify sigma to hopefully get better convergence
                if slope > 0 and correlation > 0.1:
                    if sigma > options.starting_sigma * 0.000001:
                        sigma *= 0.95

                if slope > 0 and correlation < 0.01:
                    if sigma < options.starting_sigma:
                        sigma *= 1.05
 
        # output info
        if i in logpoints:
            echo(fmt"i: {i:6d}, best: {fitness[0]:4.6f}, sigma: {sigma:4.6f}, corr: {correlation:1.6f}, slope: {slope:1.6f},")

    # return best individual as well as time and convergence data
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
