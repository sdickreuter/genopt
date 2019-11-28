
type
    Genoptions* = object
        crossover_size*: float
        mutation_rate*: float
        report_every*: int
        max_iter*: int
        starting_sigma*: float
        alpha*: float
        mutate_percent*: bool

        max_val*: float
        min_val*: float


proc initGenoptions*(): Genoptions =
    result.crossover_size = 0.7
    result.mutation_rate = 0.5
    result.report_every = 50
    result.max_iter = 5000
    result.starting_sigma = 4.0
    result.alpha = 0.66
    result.mutate_percent = true

    result.max_val = 10.0
    result.min_val = -10.0
