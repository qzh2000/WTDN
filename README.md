# The meaning of symbols in the code
      
        # model parameter        
        self.I = 11 # collection centres index
        self.J = 37 # sorting centers index
        self.K = 30 # incinerators index
        self.W = 30 # landfills index
        self.L = 3 # facility sizes index
        self.N = 8 # Uncertainty index
        self.M = 100000 # A large number
        self.C1 = [50, 100, 300] # Capacities of sorting centers
        self.C2 = [50, 100, 150] # Capacities of incinerators
        self.C3 = [50, 100, 150] # Capacities of landfills
        self.Cbar1 = 24 # Transport capacity of each truck on route (c, s)
        self.Cbar2 = 48 # Transport capacity of each truck on route (s, i) and (s, l)
        self.d = 4.7755 * 1e-6 # DALYs per person caused by transport waste
        self.d1 = [0.07, 0.14, 0.28] # DALYs per person caused by constructing sorting centers with size t
        self.d2 = [5.95, 11.9, 17.85] # DALYs per person caused by constructing incinerators with size t
        self.d3 = [3.89, 7.78, 11.66] # DALYs per person caused by constructing landfills with size t
        self.tc1 = [] # Transport costs for each fully loaded truck on routes (c, s)
        self.tc2 = [] # Transport costs for each fully loaded truck on routes (s, i)
        self.tc3 = [] # Transport costs for each fully loaded truck on routes (s, l)
        self.p1 = [] # Affected population near routes (c, s)
        self.p2 = [] # Affected population near routes (s, i)
        self.p3 = [] # Affected population near routes (s, l)
        self.P1 = [] # Affected population near collection center with size t;
        self.P2 = [] # Affected population near incinerator with size t;
        self.P3 = [] # Affected population near landfill with size t;
        self.c1 = [] # Nominal land-use stress ratio of constructing collection center with size t
        self.c2 = [] # Nominal land-use stress ratio of constructing incinerator with size t
        self.c3 = [] # Nominal land-use stress ratio of constructing landfill with size t
        self.cn1 = [] # Fluctuating land-use stress ratio of constructing collection center with size t
        self.cn2 = [] # Fluctuating land-use stress ratio of constructing incinerator with size t
        self.cn3 = [] # Fluctuating land-use stress ratio of constructing landfill with size t
        self.D = [458, 266, 485, 395, 192, 834, 230, 147, 169, 332, 501]  # Amount of waste at collection center
        self.sigma = ([0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]) # parameter in ambuguity set
        self.C = 2.5 # epsilon
        self.Gamma = 5.3 # Fluctuation range parameter

        # initialize variables
        self.y1 = {} # 0-1 variables, if a sortin centers with size t is built
        self.y2 = {} # 0-1 variables, if a incinerator with size t is built
        self.y3 = {} # 0-1 variables, if a landfill with size t is built
        self.f1 = {} # Amount of waste transported on routes (c, s)
        self.f2 = {} # Amount of waste transported on routes (s, i)
        self.f3 = {} # Amount of waste transported on routes (s, l)
        self.x1 = {} # Waste transport trucks required on routes (c, s)
        self.x2 = {} # Waste transport trucks required on routes (s, i)
        self.x3 = {} # Waste transport trucks required on routes (s, l).
 
        # Auxiliary variable about SP
        self.maxshipcost = {}

        # Auxiliary variable about lower level model
        self.alpha = {}
        self.beta = {}
        self.gamma = {}
        self.delt = {}
        self.epsilon = {}
        self.eta = {}
        self.theta = {}
        self.lamda = {}
        self.Omega1 = {}
        self.Omega2 = {}
        self.Omega3 = {}

        # Auxiliary variable about ambiguous chance constraint
        self.c = {}
        self.u10 = {}
        self.v10 = {}
        self.u1 = {}
        self.v1 = {}
        self.u2 = {}
        self.v2 = {}
        self.u3 = {}
        self.v3 = {}
        self.beta1 = {}
        self.beta2 = {}
        self.beta3 = {}
        self.beta4 = {}
        self.beta5 = {}
        self.beta6 = {}
