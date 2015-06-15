import rpy2.robjects as robjects
robjects.r("library(plant)")
robjects.r("library(plant.assembly)")

class TreeModel():
    def __init__(self, time_disturbance, slope):
        self.time_disturbance = time_disturbance
        self.slope = slope
        self.sys = None
    def evolve(self, nsteps=20, lma=[], seed_rain=[], verbose=False):
        """NOTE: This starts from an *empty* community"""
        f = robjects.r["for_python_evolve"]
        self.sys = f(self.time_disturbance, self.slope, nsteps,
                     robjects.FloatVector(lma),
                     robjects.FloatVector(seed_rain),
                     verbose)
    def equilibrium(self, lma, seed_rain, nsteps=20, verbose=False):
        assert len(lma) == len(seed_rain)
        f = robjects.r['for_python_equilibrium']
        self.sys = f(self.time_disturbance,
                     self.slope,
                     robjects.FloatVector(lma),
                     robjects.FloatVector(seed_rain),
                     nsteps,
                     verbose)
    def fitness(self, lma):
        f = robjects.r['for_python_fitness']
        res = f(self.sys, robjects.FloatVector(lma))
        return r2list(res)
    def fitness_approximate(self, lma):
        f = robjects.r['for_python_fitness_approximate']
        res = f(self.sys, robjects.FloatVector(lma))
        return r2list(res)
    def grid(self, n=50):
        """helper to generate a bunch of equally spaced points (in log
    space) within the current range"""
        res = robjects.r['seq_log_range'](self.bounds, n)
        return r2list(res)
    @property
    def bounds(self):
        return r2list(self.sys.rx2('bounds'))

# Not really clear what the right way to do this is.  With numpy:
#   import rpy2.robjects.numpy2ri as rpyn
#   vector=rpyn.ri2numpy(vector_R)
def r2list(obj):
    return [i for i in obj]
