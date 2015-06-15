import tree

time_disturbance = 11.3
slope = 1.3

# These are the 'X' values:
lma = [0.018, 0.09485587]
# These are the 'Y' values that happen to be quite good:
seed_rain = [450.587114143259, 1.31521896285294]

m = tree.TreeModel(time_disturbance, slope)

# First, compute fitness for a given set of species traits and
# densities, at equilibrium; this takes a while.
m.equilibrium(lma, seed_rain, verbose=True)

# This computes the true fitness values at the vector lma:
w = m.fitness(lma)
w

# This will construct an approximate fitness landscape.  By default
# this will use our GP approach, which samples a bunch of points
# sequentially.  It takes a while and is fairly chatty.  We can turn
# it down if that gets annoying.
w_approx = m.fitness_approximate(lma)
w_approx

# Once that's been generated then subsequent approximate fitness
# calculations are fast.
#
# Generate a series of 50 points to estimate fitness for:
lma2 = m.grid(50)
# And compute the fitness
w2 = m.fitness_approximate(lma2)
w2

# At this point I get baffled because in R I'd want to plot these, but
# plotting in Python is a black art to me :)

# Run an evolutionary assembly
m.evolve(20, lma, seed_rain, verbose=True)
m.fitness_approximate(lma2)

m.evolve(20, verbose=True)
m.fitness_approximate(lma2)
