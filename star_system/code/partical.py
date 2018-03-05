'''
build partical class
'''
import math
import numpy as np
class partical(object):
    def __init__(self,loc,v,id,m=1):
        self.x,self.y = loc
        self.vx,self.vy = v
        self.id = id
        self.m = m
        return
    def get_info(self):
        loc = [self.x,self.y]
        v = [self.vx,self.vy]
        return [loc,v]
    def set_v(self,new_v):
        self.vx,self.vy = new_v
    def cal_r(self):
        # calculate radius
        self.r = math.sqrt(self.x**2 + self.y**2)
    def cal_Ek(self):
        # calculate kinetic energy
        return 0.5*self.m*(self.vx**2 + self.vy**2)
    def cal_Ep(self):
        # calculate potential energy
        self.cal_r()
        return -self.m/self.r
    def cal_E(self):
        # calculate the partical's energy
        return self.cal_Ek() + self.cal_Ep()
    def update_state_RK(self,dt):
        # update the partical state from "t0" to "t0 + dt"
        # using Runge-Kutta method
        def F(Y):
            rr,vv = Y
            rx,ry = rr
            norm_rr = math.sqrt(rx**2+ry**2)
            r_rr = vv
            r_vv = -rr/(norm_rr**3)
            return np.array([r_rr,r_vv])

        self.cal_r()
        rr = np.array([self.x,self.y])
        vv = np.array([self.vx,self.vy])
        Yn = np.array([rr,vv])
        K1 = F(Yn)
        K2 = F(Yn + 0.5*dt*K1)
        K3 = F(Yn + 0.5*dt*K2)
        K4 = F(Yn + dt*K3)
        Y_r = Yn + dt*(K1 + 2*K2 + 2*K3 + K4)/6.
        next_loc,next_v = Y_r
        self.x,self.y = next_loc
        self.vx,self.vy = next_v




