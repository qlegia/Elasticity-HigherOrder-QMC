
class weight(object):
    def __init__(self, wp):
        assert(len(wp) <= len(self.paramlist))
        for i in range(len(wp)):
            setattr(self,self.paramlist[i],wp[i])

####################################################################################################
## PROD
####################################################################################################

class PROD_tz(weight):
    def __init__(self,wp):
        self.name = "PROD"
        self.archivename = "PROD_TZ"
        self.paramlist = ["theta","zeta"]
        super(PROD_tz,self).__init__(wp)
    def __call__(self,j):
        assert(self.zeta >= 1)
        return self.theta * double(j)**(-self.zeta)
    def __str__(self):
        return self.name+"_t%.3e_z%.3e"%(self.theta,self.zeta)
    def nice_string(self):
        return self.archivename+"_t%s_z%s"%(str(self.theta),str(self.zeta))

class PROD_INV_tz(PROD_tz):
    def __init__(self, wp):
        super(PROD_INV_tz,self).__init__(wp)
        self.name = "PROD_INV"
        self.archivename = "PROD_INV"

class PROD_Haar(weight):
    def __init__(self,wp):
        self.name = "PROD_Haar"
        self.archivename = "PROD_Haar"
        self.paramlist = ["c1","c2","beta"]
        super(PROD_Haar,self).__init__(wp)
    def __call__(self,j):
        lk = lambda j: (int(log2(j)), j-2**int(log2(j)))
        l,k = lk(j)
        return self.c1 / (1. + self.c2 * 2.**(self.beta*l))
    def __str__(self):
        return self.name+"_c1%.3e_c2%.3e_b%.3e"%(self.c1,self.c2,self.beta)
    def nice_string(self):
        return self.archivename+"_c1%s_c2%s_b%s"%(str(self.c1),str(self.c2),str(self.beta))


####################################################################################################
## SPOD
####################################################################################################

class SPOD_tz(weight):
    """Additional parameter q is from (|nu|+q)!"""
    def __init__(self,wp):
        self.name = "SPOD"
        self.archivename = "SPOD_TZ"
        self.paramlist = ["theta","zeta","q"]
        super(SPOD_tz,self).__init__(wp)
    def __call__(self,j):
        assert(self.zeta >= 1)
        return self.theta * double(j)**(-self.zeta)
    def __str__(self):
        ret = self.name+"_t%.3e_z%.3e"%(self.theta,self.zeta)
        if self.q != 0: ret += "_q%d"%self.q
        return ret
    def nice_string(self):
        ret = self.archivename+"_t%s_z%s"%(str(self.theta),str(self.zeta))
        if self.q != 0: ret += "_q%d"%self.q
        return ret

class SPOD_2dstaircase(weight):
    def __init__(self,wp):
        self.name = "SPOD_2dstaircase"
        self.paramlist = ["t","q","zeta"]
        super(SPOD_2dstaircase,self).__init__(wp)
    def __str__(self):
        ret = self.name+"_t%d"%self.t
        if self.q != 0: ret += "_q%d"%self.q
        return ret
    def nice_string(self):
        return str(self)
    ## would need indices for call
    #def __call__(self,j):
    #    return double(j)**(-self.zeta)

weight_map = {
    "PROD_TZ": PROD_tz,
    "PROD_INV_TZ": PROD_INV_tz,
    "PROD_HAAR": PROD_Haar,
    "SPOD_TZ": SPOD_tz,
    "SPOD_2DSTAIRCASE": SPOD_2dstaircase
}

def get_weight_params(data):
    c = weight_map[data["weight"].upper()]([])
    wp = []
    for p in c.paramlist:
        if "weight_%s"%p in data:
            wp.append(data["weight_%s"%p])
    return wp

def get_weight(w,wp):
    w = w.upper()
    if not w in weight_map.keys():
        raise Exception("No weight type found for: "+w)
    return weight_map[w](wp)

if __name__ == "__main__":
    w = PROD_tz([0.1,3])
    print(str(w))
    print( w.theta, w.zeta)

    w = SPOD_tz([0.1,3,0])
    print(str(w))
    w = SPOD_tz([0.1,3,1])
    print(str(w))

    w = SPOD_2dstaircase([1,2,3])
    print(str(w))

