import config as cf
import numpy as np
import math
import time

map_ped = []
evt_list = []
hits_list = []
tracks2D_list = []
tracks3D_list = []

data = np.zeros((cf.n_View, cf.n_ChanPerView, cf.n_Sample)) #view, vchan


"""
the mask will be used to differentiate background (True for noise processing) from signal (False for noise processing)
at first everything is considered background (all at True)
"""
mask = np.ones((cf.n_View, cf.n_ChanPerView, cf.n_Sample), dtype=bool)

"""
alive_chan mask intends to not take into account broken channels
True : not broken
False : broken
"""
alive_chan = np.ones((cf.n_View, cf.n_ChanPerView, cf.n_Sample), dtype=bool)


ped_rms = np.zeros((cf.n_View, cf.n_ChanPerView)) 
ped_mean = np.zeros((cf.n_View, cf.n_ChanPerView)) 

def reset_event():
    data[:,:,:] = 0.
    mask[:,:,:] = True
    ped_rms[:,:] = 0.
    ped_mean[:,:] = 0.

    hits_list.clear()
    tracks2D_list.clear()
    tracks3D_list.clear()

    [x.reset() for x in map_ped]


class pdmap: 
    def __init__(self, view, vchan):
        self.view  = view
        self.vchan = vchan
        self.raw_ped   = 0.
        self.raw_rms   = -1.
        self.evt_ped   = 0.
        self.evt_rms   = -1.
        
    def set_raw_pedestal(self, ped, rms):
        self.raw_ped = ped
        self.raw_rms = rms

    def set_evt_pedestal(self, ped, rms):
        self.evt_ped = ped
        self.evt_rms = rms
        
    def reset(self):
        self.raw_ped = 0.
        self.raw_rms = -1.
        self.evt_ped = 0.
        self.evt_rms = -1.

    def get_ana_chan(self):
        return self.view, self.vchan

class event:
    def __init__(self, date, evt_nb, evt_glob, t_s, t_ms):
        self.date         = date
        self.evt_nb       = evt_nb
        self.evt_nb_glob  = evt_glob

        self.time_s    = t_s
        self.time_ms   = t_ms

        self.nHits      = np.zeros((cf.n_View), dtype=int)
        self.nClusters  = np.zeros((cf.n_View), dtype=int)
        self.nTracks2D  = np.zeros((cf.n_View), dtype=int)
        self.nTracks3D  = 0
        
    def __eq__(self, other):
        return (self.date, self.evt_nb, self.time_s, self.time_ns) == (other.date, other.evt_nb, other.time_s, other.time_ns)

    def dump(self):
        print("DATE ",self.date, " EVENT ", self.evt_nb)
        print("Taken at ", time.ctime(self.time_s), " + ", self.time_ms, " ms ")
        print(" TS = ", self.time_s, " +", self.time_ms)

    def dump_reco(self):
        print("\n~.~.~.~.~.~.~ Reconstruction Summary ~.~.~.~.~.~.~")
        for iv in range(cf.n_View):
            print(" View ", iv, " : ")
            print("\tNb of Hits :", self.nHits[iv])
            print("\tNb of Clusters :", self.nClusters[iv])
            print("\tNb of 2D tracks :", self.nTracks2D[iv])
        print("\n")
        print(" --> Nb of 3D tracks ", self.nTracks3D)
        print("\n")


class hits:
    def __init__(self, view, channel, start, stop, charge_int, max_t, max_adc, min_t, min_adc):
        self.view    = view
        self.channel = channel
        self.start   = start
        self.stop    = stop

        self.max_t   = max_t
        self.max_adc = max_adc
        self.min_t   = min_t
        self.min_adc = min_adc

        self.charge_int  = charge_int
        self.charge_max  = 0.
        self.charge_min  = 0.
        self.charge_pv   = 0. #peak-valley

        self.cluster = -1 #cluster
        self.X       = -1
        self.Z       = -1
        self.matched = -9999
        

    def __lt__(self,other):
        #"""sort hits by increasing channel and increasing Z"""
        #return (self.X < other.X) or (self.X== other.X and self.Z < other.Z)

        """ sort hits by decreasing Z and increasing channel """
        return (self.Z > other.Z) or (self.Z == other.Z and self.X < other.X)

    def hit_positions(self, v):
        self.X = self.channel*cf.ChanPitch
        if(self.view == 0):
            self.Z = cf.Anode_Z - self.max_t*cf.n_Sampling*v*0.1
        else:
            self.Z = cf.Anode_Z - self.min_t*cf.n_Sampling*v*0.1

    def hit_charge(self):
        self.charge_int *= cf.n_Sampling / cf.ADCtofC

        self.charge_max = (self.max_adc) * cf.n_Sampling / cf.ADCtofC
        self.charge_min = (self.min_adc) * cf.n_Sampling / cf.ADCtofC
        self.charge_pv  = (self.max_adc - self.min_adc) * cf.n_Sampling / cf.ADCtofC

    def set_match(self, ID):
        self.matched = ID

    def set_cluster(self, ID):
        self.cluster = ID

    def get_charges(self):
        return (self.charge_int, self.charge_max, self.charge_min, self.charge_pv)

class trk2D:
    def __init__(self, ID, view, ini_slope, ini_slope_err, x0, y0, q0, chi2, cluster):
        self.trackID = ID
        self.view    = view
    
        self.ini_slope       = ini_slope
        self.ini_slope_err   = ini_slope_err
        self.end_slope       = ini_slope
        self.end_slope_err   = ini_slope_err

        self.nHits      = 1
        self.nHits_dray = 0

        self.path    = [(x0,y0)]
        self.dQ      = [q0]

        self.chi2_fwd    = chi2
        self.chi2_bkwd   = chi2

        self.drays   = []
        
        self.tot_charge_int = q0[0]
        self.tot_charge_max = q0[1]
        self.tot_charge_min = q0[2]
        self.tot_charge_pv  = q0[3]

        self.dray_charge_int = 0.
        self.dray_charge_max = 0.
        self.dray_charge_min = 0.
        self.dray_charge_pv = 0.

        self.len_straight = 0.
        self.len_path = 0.

        self.matched = -1
        self.cluster = cluster
        
    def __lt__(self,other):
        """ sort tracks by decreasing Z and increasing channel """
        return (self.path[0][1] > other.path[0][1]) or (self.path[0][1] == other.path[0][1] and self.path[0][0] < other.path[0][0])


    def add_drays(self, x, y, q):
        self.drays.append((x,y,q[0], q[1], q[2], q[3]))
        self.dray_charge_int += q[0]
        self.dray_charge_max += q[1]
        self.dray_charge_min += q[2]
        self.dray_charge_pv  += q[3]
        self.nHits_dray += 1
        self.remove_hit(x, y, q)


    def remove_hit(self, x, y, q):
        pos = -1
        for p,t in enumerate(self.path):
            if(t[0] == x and t[1] == y and self.dQ[p]==q):
                pos = p
                break
        #print("removing at position ", pos)
        if(pos >= 0):
            self.path.pop(pos)
            self.dQ.pop(pos)
            self.nHits -= 1
            self.tot_charge_int -= q[0]
            self.tot_charge_max -= q[1]
            self.tot_charge_min -= q[2]
            self.tot_charge_pv  -= q[3]
        else:
            print("?! cannot remove hit ", x, " ", y, " ", q)
            
    def add_hit(self, x, y, q):
        self.nHits += 1
        
        self.len_path += math.sqrt( pow(self.path[-1][0]-x, 2) + pow(self.path[-1][1]-y,2) )
        #beware to append (x,y) after !
        self.path.append((x,y))
        self.dQ.append(q)
        self.tot_charge_int += q[0]
        self.tot_charge_max += q[1]
        self.tot_charge_min += q[2]
        self.tot_charge_pv  += q[3]
        self.len_straight = math.sqrt( pow(self.path[0][0]-self.path[-1][0], 2) + pow(self.path[0][1]-self.path[-1][1], 2) )


    def add_hit_update(self, slope, slope_err, x, y, q, chi2):
        self.end_slope = slope
        self.end_slope_err = slope_err
        self.nHits += 1
        self.len_path += math.sqrt( pow(self.path[-1][0]-x, 2) + pow(self.path[-1][1]-y,2) )

        #beware to append (x,y) after !
        self.path.append((x,y))
        self.dQ.append(q)
        self.chi2 = chi2
        #self.tot_charge += q
        self.tot_charge_int += q[0]
        self.tot_charge_max += q[1]
        self.tot_charge_min += q[2]
        self.tot_charge_pv  += q[3]
        self.len_straight = math.sqrt( pow(self.path[0][0]-self.path[-1][0], 2) + pow(self.path[0][1]-self.path[-1][1], 2) )


    def update_forward(self, chi2, slope, slope_err):
        self.chi2_fwd = chi2
        self.end_slope = slope
        self.end_slope_err = slope_err

    def update_backward(self, chi2, slope, slope_err):
        self.chi2_bkwd = chi2
        self.ini_slope = slope
        self.ini_slope_err = slope_err

    def reset_path(self, path, dQ):
        self.path = path
        self.dQ = dQ
        self.finalize_track()
        

    def finalize_track(self):
        self.nHits = len(self.path)

        self.tot_charge_int = sum(i for i,j,k,l in self.dQ)
        self.tot_charge_max = sum(j for i,j,k,l in self.dQ)
        self.tot_charge_min = sum(k for i,j,k,l in self.dQ)
        self.tot_charge_pv  = sum(l for i,j,k,l in self.dQ)

        self.nHits_dray = len(self.drays)
        self.dray_charge_int = sum(k for i,j,k,l,m,n in self.drays)
        self.dray_charge_max = sum(l for i,j,k,l,m,n in self.drays)
        self.dray_charge_min = sum(m for i,j,k,l,m,n in self.drays)
        self.dray_charge_pv  = sum(n for i,j,k,l,m,n in self.drays)


        self.len_straight = math.sqrt( pow(self.path[0][0]-self.path[-1][0], 2) + pow(self.path[0][1]-self.path[-1][1], 2) )        
        self.len_path = 0.
        for i in range(self.nHits-1):
            self.len_path +=  math.sqrt( pow(self.path[i][0]-self.path[i+1][0], 2) + pow(self.path[i][1]-self.path[i+1][1],2) )
            
            

    def dist(self, other, i=-1, j=0):
        return math.sqrt(pow( self.path[i][0] - other.path[j][0], 2) + pow(self.path[i][1] - other.path[j][1], 2))



    def slope_comp(self, other):#, sigcut):
        """ check if both tracks have the same slope direction """
        if(self.end_slope * other.ini_slope < 0.):
            return 9999. #False

        """ if slope error is too low, re-assign it to 5 percent """
        if(self.end_slope_err == 0 or math.fabs(self.end_slope_err/self.end_slope) < 0.05):
            end_err = math.fabs(self.end_slope*0.05)
        else:
            end_err = self.end_slope_err

        if(other.ini_slope_err == 0 or math.fabs(other.ini_slope_err/other.ini_slope) < 0.05):
            ini_err = math.fabs(other.ini_slope*0.05)
        else:
            ini_err = other.ini_slope_err

        #return (math.fabs( self.end_slope - other.ini_slope) < (sigcut*end_err + sigcut*ini_err))

        return math.fabs( self.end_slope - other.ini_slope) / (end_err + ini_err)


    def x_extrapolate(self, other, rcut):
        xa, za = self.path[-1][0], self.path[-1][1]
        xb, zb = other.path[0][0], other.path[0][1]

        return ( math.fabs( xb - (xa+(zb-za)*self.end_slope)) < rcut) and (math.fabs( xa-(xb+(za-zb)*other.ini_slope)) < rcut)

    def z_extrapolate(self, other, rcut):
        xa, za = self.path[-1][0], self.path[-1][1]
        xb, zb = other.path[0][0], other.path[0][1]
        if(self.end_slope == 0 and other.ini_slope == 0) : 
            return True

        if(self.end_slope == 0):
            return (math.fabs(za - zb - (xa-xb)/other.ini_slope) < rcut)
        elif( other.ini_slope == 0):
            return ( math.fabs(zb - za - (xb-xa)/self.end_slope) < rcut)
        else:
            return ( math.fabs(zb - za - (xb-xa)/self.end_slope) < rcut) and (math.fabs(za - zb - (xa-xb)/other.ini_slope) < rcut)


    def joinable(self, other, dcut, sigcut, rcut):
        if(self.view != other.view): 
            return False
        if( self.dist(other) < dcut and self.slope_comp(other) <  sigcut and self.x_extrapolate(other, rcut) and self.z_extrapolate(other, rcut)):            
            return True


    def merge(self, other):
        self.nHits += other.nHits
        self.nHits_dray += other.nHits_dray
        self.chi2 += other.chi2 #should be refiltered though

        self.tot_charge_int += other.tot_charge_int
        self.tot_charge_max += other.tot_charge_max
        self.tot_charge_min += other.tot_charge_min
        self.tot_charge_pv  += other.tot_charge_pv 

        self.dray_charge_int += other.dray_charge_int
        self.dray_charge_max += other.dray_charge_max
        self.dray_charge_min += other.dray_charge_min
        self.dray_charge_pv  += other.dray_charge_pv 

        self.len_path += other.len_path 
        self.len_path += self.dist(other)
        self.matched = -1
        self.drays.extend(other.drays)

        if(self.path[0][1] > other.path[0][1]):
               self.ini_slope = self.ini_slope
               self.ini_slope_err = self.ini_slope_err
               self.end_slope = other.end_slope
               self.end_slope_err = other.end_slope_err
               
               self.path.extend(other.path)
               self.dQ.extend(other.dQ)
               self.len_straight = math.sqrt( pow(self.path[0][0]-self.path[-1][0], 2) + pow(self.path[0][1]-self.path[-1][1], 2) )

        else:
               self.ini_slope = other.ini_slope
               self.ini_slope_err = other.ini_slope_err
               self.end_slope = self.end_slope
               self.end_slope_err = self.end_slope_err

               self.path = other.path + self.path
               self.dQ = other.dQ + self.dQ
               self.len_straight = math.sqrt( pow(self.path[0][0]-self.path[-1][0], 2) + pow(self.path[0][1]-self.path[-1][1],2) )        
        

    def mini_dump(self):
        print("[",self.view,"] from (%.1f,%.1f)"%(self.path[0][0], self.path[0][1]), " to (%.1f, %.1f)"%(self.path[-1][0], self.path[-1][1]), " N = ", self.nHits, " L = %.1f/%.1f"%(self.len_straight, self.len_path), " Q = ", self.tot_charge_int )
               

class trk3D:
    def __init__(self, tv0, tv1):

        self.chi2    = 0.5*(tv0.chi2 + tv1.chi2)
        self.momentum = -1

        self.nHits_v0   = tv0.nHits
        self.nHits_v1   = tv1.nHits

        self.len_straight_v0 = -1#tv0.len_straight
        self.len_straight_v1 = -1#tv1.len_straight

        self.len_path_v0 = -1 #tv0.len_path
        self.len_path_v1 = -1 #tv1.len_path

        self.tot_charge_int_v0 = -1
        self.tot_charge_int_v1 = -1
        self.tot_charge_max_v0 = -1
        self.tot_charge_max_v1 = -1
        self.tot_charge_min_v0 = -1
        self.tot_charge_min_v1 = -1
        self.tot_charge_pv_v0  = -1
        self.tot_charge_pv_v1  = -1

        self.ini_theta = -1
        self.end_theta = -1
        self.ini_phi = -1
        self.end_phi = -1

        self.t0_corr = 0.
        self.z0_corr = 0.

        self.ini_x = tv0.path[0][0]
        self.ini_y = tv1.path[0][0]
        self.ini_z = 0.5*(tv0.path[0][1] + tv1.path[0][1])

        self.end_x = tv0.path[-1][0]
        self.end_y = tv1.path[-1][0]
        self.end_z = 0.5*(tv0.path[-1][1] + tv1.path[-1][1])

        
        self.path_v0 = []
        self.path_v1 = []

        self.dQds_int_v0 = []
        self.dQds_int_v1 = []
        self.dQds_max_v0 = []
        self.dQds_max_v1 = []
        self.dQds_min_v0 = []
        self.dQds_min_v1 = []
        self.dQds_pv_v0  = []
        self.dQds_pv_v1  = []

        self.ds_v0 = []
        self.ds_v1 = []


    def set_view0(self, path, dqds_int, dqds_max, dqds_min, dqds_pv, ds):
        self.len_straight_v0 = math.sqrt( pow(path[0][0]-path[-1][0], 2) + pow(path[0][1]-path[-1][1],2) + pow(path[0][2]-path[-1][2],2) )     
        self.len_path_v0 = 0.

        for i in range(len(path)-1):
            self.len_path_v0 +=  math.sqrt( pow(path[i][0]-path[i+1][0], 2) + pow(path[i][1]-path[i+1][1],2)+ pow(path[i][2]-path[i+1][2],2) )
            
        #self.len_path_v0 = length

        self.path_v0     = path

        self.dQds_int_v0     = dqds_int
        self.dQds_max_v0     = dqds_max
        self.dQds_min_v0     = dqds_min
        self.dQds_pv_v0      = dqds_pv
        self.ds_v0           = ds

        self.tot_charge_int_v0 = sum(dqds_int)
        self.tot_charge_max_v0 = sum(dqds_max)
        self.tot_charge_min_v0 = sum(dqds_min)
        self.tot_charge_pv_v0  = sum(dqds_pv)


    def set_view1(self, path, dqds_int, dqds_max, dqds_min, dqds_pv, ds):
        self.len_straight_v1 = math.sqrt( pow(path[0][0]-path[-1][0], 2) + pow(path[0][1]-path[-1][1],2) + pow(path[0][2]-path[-1][2],2) )     

        self.len_path_v1 = 0.
        for i in range(len(path)-1):
            self.len_path_v1 +=  math.sqrt( pow(path[i][0]-path[i+1][0], 2) + pow(path[i][1]-path[i+1][1],2)+ pow(path[i][2]-path[i+1][2],2) )

        self.path_v1     = path

        self.dQds_int_v1     = dqds_int
        self.dQds_max_v1     = dqds_max
        self.dQds_min_v1     = dqds_min
        self.dQds_pv_v1      = dqds_pv

        self.ds_v1           = ds

        self.tot_charge_int_v1 = sum(dqds_int)
        self.tot_charge_max_v1 = sum(dqds_max)
        self.tot_charge_min_v1 = sum(dqds_min)
        self.tot_charge_pv_v1  = sum(dqds_pv)



        #self.dQds_v1     = dqds
        #self.tot_charge_v1 = sum(dqds)

    def matched(self, tv0, tv1):
        tv0.matched = evt_list[-1].nTracks3D
        tv1.matched = evt_list[-1].nTracks3D

    def set_t0_z0_corr(self, t0, z0):
        self.t0_corr = t0
        self.z0_corr = z0

    def angles(self, tv0, tv1):

        """ initial angles """
        slope_v0 = tv0.ini_slope #dx/dz
        slope_v1 = tv1.ini_slope #dy/dz
        self.ini_phi = math.degrees(math.atan2(slope_v1, slope_v0))
        self.ini_theta = math.degrees(math.atan2(math.sqrt(pow(slope_v0,2)+pow(slope_v1,2)),-1.))

        """ end angles """
        slope_v0 = tv0.end_slope #dx/dz
        slope_v1 = tv1.end_slope #dy/dz
        self.end_phi = math.degrees(math.atan2(slope_v1, slope_v0))
        self.end_theta = math.degrees(math.atan2(math.sqrt(pow(slope_v0,2)+pow(slope_v1,2)),-1.))

        
    def dump(self):
        print(" from (%.2f, %.2f, %.2f) to (%.2f, %.2f, %.2f)"%(self.ini_x, self.ini_y, self.ini_z, self.end_x, self.end_y, self.end_z))
        print(" theta, phi: [ini] %.2f ; %.2f"%(self.ini_theta, self.ini_phi), " -> [end] %.2f ; %.2f "%( self.end_theta, self.end_phi), " L = (P) %.2f / %.2f ; (S) %.2f / %.2f"%(self.len_path_v0, self.len_path_v1, self.len_straight_v0, self.len_straight_v1))
        print(" corr : %.2f cm / %.2f mus"%(self.z0_corr, self.t0_corr))
        print(" charge int V0: %.2f, V1: %.2f A= %.3f"%(self.tot_charge_int_v0, self.tot_charge_int_v1, (self.tot_charge_int_v0-self.tot_charge_int_v1)/(self.tot_charge_int_v0+self.tot_charge_int_v1)))
        print(" charge max V0: %.2f, V1: %.2f A= %.3f"%(self.tot_charge_max_v0, self.tot_charge_max_v1, (self.tot_charge_max_v0-self.tot_charge_max_v1)/(self.tot_charge_max_v0+self.tot_charge_max_v1)))
        print(" charge min V0: %.2f, V1: %.2f A= %.3f"%(self.tot_charge_min_v0, self.tot_charge_min_v1, (self.tot_charge_min_v0-self.tot_charge_min_v1)/(self.tot_charge_min_v0+self.tot_charge_min_v1)))
        print(" charge pv V0: %.2f, V1: %.2f A= %.3f"%(self.tot_charge_pv_v0, self.tot_charge_pv_v1, (self.tot_charge_pv_v0-self.tot_charge_pv_v1)/(self.tot_charge_pv_v0+self.tot_charge_pv_v1)))
