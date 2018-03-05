from partical import partical
import math
import numpy as np
import matplotlib.pyplot as plt
import random
def test_single_partical():
    loc = [1, 0]
    v = [0, 0.8]
    id = 1
    p1 = partical(loc, v, id)
    dt = 0.01
    tot_time_step = 1000
    trace_x = []
    trace_y = []
    Ep_hist = []
    Ek_hist = []
    for time_step in range(tot_time_step):
        tmp_loc, tmp_v = p1.get_info()
        tmp_x, tmp_y = tmp_loc
        tmp_vx, tmp_vy = tmp_v
        trace_x.append(tmp_x)
        trace_y.append(tmp_y)
        Ep_hist.append(p1.cal_Ep())
        Ek_hist.append(p1.cal_Ek())
        p1.update_state_RK(dt)
    Ep_hist = np.array(Ep_hist)
    Ek_hist = np.array(Ek_hist)
    E_hist = Ek_hist + Ep_hist

    plt.clf()
    plt.scatter(trace_x, trace_y)
    plt.xlabel("x")
    plt.ylabel("y")
    # plt.show()
    plt.savefig("single_partical_trace.png")
    plt.clf()
    plt.plot(np.arange(tot_time_step), Ep_hist, label="Ep")
    plt.plot(np.arange(tot_time_step), Ek_hist, label="Ek")
    plt.plot(np.arange(tot_time_step), E_hist, label="E")
    lgd = plt.legend(bbox_to_anchor=(1.04, 1), borderaxespad=0)
    plt.xlabel("iter")
    plt.ylabel("Energy")
    # plt.show()
    plt.savefig("single_partical_E.png")

class star_system(object):
    def __init__(self,number_of_partical,radius_range,percent_dist = 0,partical_init_mass = 1):
        self.number_of_partical = number_of_partical
        self.radius_range = radius_range
        self.percent_dist = percent_dist
        self.partical_init_mass = partical_init_mass
        self.partical_dic = {}
        self.trace_x_dic = {}
        self.trace_y_dic = {}
        self.Ep_hist_dic = {}
        self.Ek_hist_dic = {}
        def init_loc():
            r_in,r_out = self.radius_range
            tmp_rand = random.random()
            r = r_out*tmp_rand + r_in*(1-tmp_rand)
            tmp_rand = random.random()
            theta = tmp_rand*2*math.pi
            loc = np.array([r*math.cos(theta),r*math.sin(theta)])
            return loc

        def cal_init_v(loc):
            x,y = loc
            r = math.sqrt(x**2 + y**2)
            v = np.array([-y,x])/(r**1.5)
            return v

        for partical_id in range(self.number_of_partical):
            tmp_loc = init_loc()
            tmp_v = cal_init_v(tmp_loc)
            self.partical_dic[partical_id] = partical(tmp_loc,tmp_v,partical_id)
            self.trace_x_dic[partical_id] = []
            self.trace_y_dic[partical_id] = []
            self.Ek_hist_dic[partical_id] = []
            self.Ep_hist_dic[partical_id] = []

    def update_particals(self,dt=0.01):
        for partical_id in self.partical_dic:
            p1 = self.partical_dic[partical_id]
            p1.update_state_RK(dt)

            tmp_loc, tmp_v = p1.get_info()
            tmp_x, tmp_y = tmp_loc
            tmp_vx, tmp_vy = tmp_v
            self.trace_x_dic[partical_id].append(tmp_x)
            self.trace_y_dic[partical_id].append(tmp_y)
            self.Ep_hist_dic[partical_id].append(p1.cal_Ep())
            self.Ek_hist_dic[partical_id].append(p1.cal_Ek())

    def draw_trace(self):
        plt.clf()
        plt.title("trace")
        for partical_id in self.trace_x_dic:
            x = self.trace_x_dic[partical_id]
            y = self.trace_y_dic[partical_id]
            plt.scatter(x,y,s=3,marker=",")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.savefig("muti_partical_without_crash.png")

test_single_partical()
exit(0)

sys = star_system(10,[2,5])
tot_time_step = 5000
for tmp_step in range(tot_time_step):
    sys.update_particals()
sys.draw_trace()








