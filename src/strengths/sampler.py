class SamplerBase:
    def __init__(self):
        self.t = []
        self.x = []
        
    def requires_sample(self, t):
        return NotImplementedError()
    
    def sample(self, t, state):
        self.t.append(t)
        self.x.append(state)
        
class OnTimeSampler(SamplerBase):
    def __init__(self, t_sample):
        SamplerBase.__init__(self)
        self.t_sample = t_sample
        self.pos = 0
    
    def requires_sample(self, t):
        if self.pos == len(self.t_sample):
            return False
        sampling_is_required = False
        while t>=self.t_sample[self.pos]:
            sampling_is_required = True
            self.pos+=1
            if self.pos == len(self.t_sample):
                break
        return sampling_is_required

class OnIntervalSampler(SamplerBase):
    def __init__(self, interval):
        SamplerBase.__init__(self)
        self.interval = interval
        self.pos = -1
    
    def requires_sample(self, t):
        current_pos = int(t/self.interval)
        if current_pos>self.pos:
            self.pos = current_pos
            return True
        return False
        
class OnIterationSampler(SamplerBase):
    def __init__(self):
        SamplerBase.__init__(self)
    
    def requires_sample(self, t):
        return True
        
class ManualSampler(SamplerBase):
    def __init__(self):
        SamplerBase.__init__(self)
    
    def requires_sample(self, t):
        return False

def create_sampler(
    policy,
    t_sample,
    interval,
    t_max
    ):
    if policy == "on_t_sample":  return OnTimeSampler(t_sample)
    if policy == "on_interval":  return OnIntervalSampler(interval)
    if policy == "on_iteration": return OnIterationSampler()
    if policy == "no_sampling":  return ManualSampler()
