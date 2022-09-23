"""
Spyder Editor

This is a temporary script file.
"""

import math as mt
import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import networkx as nx
import math
import statistics as st
import pickle
from scipy import stats



"""
https://dsp.stackexchange.com/questions/58420/calculating-data-rate-using-bandwidth-transmission-power-noise-power-spectrum

P_rx = P_tx * G_tx * G * G_rx    

R = B * log(1+ P_rx/P_N)


SE = R/B = log(1+ P_rx/P_N)

A 0.6-m-diameter transmit antenna was incorporated on the LEO satellite 
with 40 W transmit power and effective isotropic radiated power (EIRP) of 
57 dBW. Nov 30, 2016

Your thermal noise power density is given as -174 dBm/Hz. 
This means 10**(−17.4) W/Hz which at a bandwidth of 10 MHz (10e6) becomes 10−10.4=0.04nW.
"""

def cal_rate(p_rx, bw):
#    -174 dBm/Hz = 10**(−17.4) W/Hz
    p_n = (10**(-17.4))*bw #bandwidth at 10 MHz
    return bw * mt.log10(1+(p_rx/p_n))

def spec_eff(p_rx):
    p_n = 0.04e-9
    return mt.log10(1+(p_rx/p_n))


def cal_gain_db(h, x, y, z, ox, oy, oz, fc, rician_fading, atten=0.05):
    d = sat_dist(x, y, z, ox, oy, oz)
    chi = pitch_angle(h, x, y, ox, oy)
    a_fad = angle_fading(chi)
    atm_fad = atm_fading(d, h, atten=0.05)
    g = channel(d, angle_fading=a_fad, atm_fading=atm_fad, fc=fc, rician_fading=rician_fading)
    return 10*mt.log10(g)
#    return g

def cal_gain(h, x, y, z, ox, oy, oz, fc, rician_fading, atten=0.05):
    d = sat_dist(x, y, z, ox, oy, oz)
    chi = pitch_angle(h, x, y, ox, oy)
    a_fad = angle_fading(chi)
    atm_fad = atm_fading(d, h, atten=0.05)
    a_fad=1
    g = channel(d, angle_fading=a_fad, atm_fading=atm_fad, fc=fc, rician_fading=rician_fading)
    return g
"""
     Channel Gain
     d:             distance sat-ue
     angle_fading:  pitch angle fading
     atm_fading:    atmospheric fading
     c:             speed of light
     fc:            carrier frequency
     rician_fading: Rician small-scale fading -> 20 db for satellite
"""
def channel(d, angle_fading, atm_fading, fc, rician_fading, c=3e8):
    base = (c/(4*mt.pi*d*fc))**2
    G = base*angle_fading*atm_fading*rician_fading
    return G


"""
     Channel Gain
     d:             distance sat-ue
     h: sat height
     atten: the attenuation through the clouds and rain in dB/km
     https://www.intechopen.com/chapters/53573
"""
def atm_fading(d, h, atten=0.05):
    fad = (3*d*atten)/(10*h)
    return 10**fad

"""
     Distance between sat and ue
     h:      sat height 
     x,y:    ue location
     ox, oy: beam coverage centre locaiton
     return distance in m
"""
def sat_dist(x, y, z, ox, oy, oz):
    dxs=(x-ox)**2
    dys=(y-oy)**2
    hs=(z-oz)**2
    return (hs+dxs+dys)**(1/2)

"""
     Pitch angle fading
     chi: the pitch angle
     a:   the antenna aperture efficiency
     n:   the rolloff factor for the antenna - 20 for sat
"""
def angle_fading(chi, a=1, n=20):
    const = (32*mt.log10(2)) / (2*((2*mt.acos(0.5**(1/n)))**2))
#    const = 1.09803
#    return abs(a*mt.cos(chi)**n*const)
    return a*(mt.cos(chi)**n)*const


"""
    Pitch angle
     h:      sat height
     x,y:    ue location
     ox, oy: sat locaiton
"""
def pitch_angle(h, x, y, ox, oy):
    dx=(x-ox)**2
    dy=(y-oy)**2
    num=(dx*dy)**0.5
    den=2*h
    return 2*mt.atan(num/den)

"""
    propagation delay
    t=d/v
    v is the speed of light (m/s)
    d is disntance in (m)
"""
def prop_delay(d, v=3e8):
    t=d/v
    return t
    

"""
replacement for df factor
    check single handover and omit
"""
def reduce_handover(sln_og): #add threshold of df count
    sln = sln_og.copy()
#    df_count=0
    temp=[sln[0]]
    for s in sln:
        s, t = s.split('-')
        for i in sln_og:
            st, tt = i.split('-')
            if s==st and t==tt:
                pass
            elif s==st:
                temp.append(st+'-'+tt)
            else:
                temp.append(st)

    return temp


"""
    count number of handovers
    neglect last solo satellite handover
"""
def count_handovers(sln):
    #count handovers
    handover_counter = 0
    for i in range(1, len(sln)):
        s, t = sln[i].split('-')
        sp, tp = sln[i-1].split('-')
        if sp != s:
            handover_counter = handover_counter + 1
    return handover_counter

"""
filter out virtual handover
and bounc
"""
def filter_handovers(sln_og):
    sln = sln_og.copy()
    temp=[sln[0]]
    for i in range(1, len(sln)):
        sc,tc = sln[i].split('-') 
        sp,tp = sln[i-1].split('-')
        if sc != sp:
            temp.append(sc+'-'+tc)
    return temp

def search_dictionary(dic, x,y,z):
    for name, age in dic.items():  # for name, age in dictionary.iteritems():  (for Python 2.x)
        if age == (x, y, z):
            print(name)   
#            return name

"""
set up working table from starlink data directory 
index are the satellite names
columnes are t_i of targeted seconds of each sat xyz to be loaded
path: directory of all sats
T: total time duration in sec (20 min 1200)
t_i: interval size of time
"""
def load_table(path, T, t_i):
    dir_list = os.listdir(path)
    dir_list = [x[:-4] for x in dir_list]
    ts=[]
    n_col=(T//t_i)+1
    for i in range(0,n_col):
        ts.append('t'+str(i))
    
    db = pd.DataFrame(columns=ts, index=dir_list)  
    
    for rowIndex, row in db.iterrows():
        dic = pd.read_pickle(path + str(rowIndex) + '.pkl')
        n=0
        for columnIndex, value in row.items():
            db.at[rowIndex, columnIndex]=dic.get(n)
            n=n+t_i
    return db

def load_table_hybrid(path, T, t_i, sats):
    sats = [i[0:5] for i in sats] #extract sat number only
    sats = list(set(sats)) #remove duplicates   
    ts=[]
    n_col=(T//t_i)+1
    for i in range(0,n_col):
        ts.append('t'+str(i))
    
    db = pd.DataFrame(columns=ts, index=sats)  
    
    for rowIndex, row in db.iterrows():
        dic = pd.read_pickle(path + str(rowIndex) + '.pkl')
        n=0
        for columnIndex, value in row.items():
            db.at[rowIndex, columnIndex]=dic.get(n)
            n=n+t_i
    return db

"""
value = (x, y, z, dist, rate, delay, rate_norm, delay_norm)
value[0] = x
value[1] = y
value[2] = z
value[3] = dist
value[4] = rate
value[5] = delay
value[6] = rate_norm
value[7] = delay_norm

"""
def build_table(db, h, x0, y0, z0, p_tx, bw, fc, rician_fading, atten):
    delay_total=[]
    rate_total=[]
    for rowIndex, row in db.iterrows():
        for columnIndex, value in row.items():
            x,y,z=value
            d = sat_dist(x, y, z, x0, y0, z0)
            delay = prop_delay(d*1000)
            
            g = cal_gain(h, x, y, z, x0, y0, z0, 
                           fc=fc, rician_fading=rician_fading, atten=0.05)
            p_rx = g * p_tx
            rate = cal_rate(p_rx, bw)
            if np.isnan(rate): rate=0
            
            new_value = value + (d, rate, delay)
            db.at[rowIndex, columnIndex]=new_value
            
            delay_total.append(delay)
            rate_total.append(rate)

    rate_total = [float(i)/max(rate_total) for i in rate_total]
    delay_total = [float(i)/max(delay_total) for i in delay_total]

    for rowIndex, row in db.iterrows():
        for columnIndex, value in row.items():
            rv=rate_total.pop(0)
            dv=delay_total.pop(0)
            db.at[rowIndex, columnIndex]=value+(rv,dv)

    return db

"""
build graph and add satellites within threshold distance
threshold distance to control graph size
"""
def build_graph(db, wd, wr, th=2000, df=0, df_option='no'):
    graph = nx.DiGraph()
    cols =db.columns.tolist()
    sat = db.index.tolist()
    fn='first-node'
    ln='last-node'
    #traverse columns in reverse
    for c in range(len(cols)-1, 0,-1):
        #traverse values in columns
        for v,s in zip(db[cols[c]],sat):
            w = (1-wr*v[6]) + wd*v[7] #rate at [6], delay at [7]
#            w = v[7]
            if df_option=='yes':
                w = w + df
            if v[3]>th:#distance threshold
                continue
            for sp in sat:
                from_edge = sp + "-" + cols[c-1]
                to_edge = s + "-" + cols[c]
                if c==(len(cols)-1): #last dummy node to solve graph
                    graph.add_edge(to_edge, ln, weight=0)
                if c==1 and db.loc[sp][0][3]<th: #first dummy node to solve graph
                    graph.add_edge(fn, from_edge, weight=0)
                graph.add_edge(from_edge, to_edge, weight=w)
    return graph


def solve_graph(g, mode=1):
    start_node = 'first-node'
    end_node = 'last-node'
    if mode==2:
        return [p for p in nx.all_shortest_paths(g, start_node, end_node, weight='weight', 
                                        method='dijkstra')]
    else:
#         print (nx.dijkstra_path(g, start_node, end_node, weight='weight'))
         sln = nx.dijkstra_path(g, start_node, end_node, weight='weight')
         sln=sln[1:-1]
         st, tt = sln[1].split('-')
         sln[0] = str(int(st)-1)+'-t0'
         return sln
            
def print_graph(g):
    print(list(g.edges.data()))

"""
build_table must be called before
"""
def closest_sat(db):
    min_d_st=1e4
    min_d_fn=1e4
    min_row_st=''
    min_row_fn=''
    sat = db.index.tolist()
    cols =db.columns.tolist()
    for v,s in zip(db[cols[0]],sat):
        d=v[3]
        if d<min_d_st:
            min_d_st=d
            min_row_st=s
    for v,s in zip(db[cols[-1]],sat):
        d=v[3]
        if d<min_d_fn:
            min_d_fn=d
            min_row_fn=s
    return min_row_st, cols[0], min_row_fn, cols[-1]



"""
solution data
"""
def get_data(db, sln):
    #s,t=sln[1].split('-') returns s=sat_num, t=tx -> db.loc[s][t]
    rv=[]
    dv=[]
    pv=[]
    for i in sln:
        s,t=i.split('-')
        rv.append(db.loc[s][t][4])#tuple[4] is rate
        dv.append(db.loc[s][t][5])#tuple[5] is delay   
        pv.append(db.loc[s][t][3])#tuple[3] is dist         
    return rv, dv, pv


def pull_data(path, db, sln, T, t_i, h_sat, x0, y0, z0, p_tx, bw, fc, rician_fading, atten):
    rvc=[]
    dvc=[]   
    for i in sln:
        s,t=i.split('-')
        dic = pd.read_pickle(path + s + '.pkl')
#        x,y,z = db.loc[s][t][0],db.loc[s][t][1],db.loc[s][t][2]
        index=int(t.replace('t',''))*t_i
#        print(index-(t_i//2), index, index+(t_i//2))
        for n in range(index-(t_i//2), index+(t_i//2)):
            if n<len(dic) and n<=T:
                x,y,z = dic.get(n)
            else:
#                print('out of dic range')
                break
            d = sat_dist(x, y, z, x0, y0, z0)
            delay = prop_delay(d*1000)
            g = cal_gain(h_sat, x, y, z, x0, y0, z0, fc=fc, rician_fading=rician_fading, atten=0.05)
            p_rx = g * p_tx
            rate = cal_rate(p_rx, bw)
            rvc.append(rate)
            dvc.append(delay) 
#        print(n, rate, delay)
    return rvc, dvc

def pull_sat(path, db, sat, n1, n2, T, t_i, h_sat, x0, y0, z0, p_tx, bw, fc, rician_fading, atten):
    rvc=[]
    dvc=[]   
    if 't' in sat:
        s,t=sat.split('-')

    dic = pd.read_pickle(path + s + '.pkl')
    
    for n in range(n1, n2):
        if n<len(dic) and n<=T:
            x,y,z = dic.get(n)
        else:
#                print('out of dic range')
            break
        d = sat_dist(x, y, z, x0, y0, z0)
        delay = prop_delay(d*1000)
        g = cal_gain(h_sat, x, y, z, x0, y0, z0, fc=fc, rician_fading=rician_fading, atten=0.05)
        p_rx = g * p_tx
        rate = cal_rate(p_rx, bw)
        rvc.append(rate)
        dvc.append(delay) 
#        print(n, rate, delay)
    return rvc, dvc


"""
ELEVATION ANGLE METHOD METHOD

starts at sln[0]
iterate t_i columns of sln[0] in dbe 
    determine if the threshold is reached
        search sats in column to have the lowest threshold
    repeat
record used sats, rate, delay, distance
    
fs: first satellite
th is in km


value = (x, y, z, dist, rate, delay, rate_norm, delay_norm)
value[0] = x
value[1] = y
value[2] = z
value[3] = dist
value[4] = rate
value[5] = delay
value[6] = rate_norm
value[7] = delay_norm

xvc plotting limit
tv time vairace to skip when plotting
"""
def elev_method(path, fs, th, xvc, T, t_i, h_sat, x0, y0, z0, p_tx, bw, fc, rician_fading, atten):
#    path, dbe, fs, th=2000, xvc):
    sats,rve,dve,pve=[],[],[],[]
    dir_list = os.listdir(path)
    dir_list = [x[:-4] for x in dir_list]
    sat,ts=fs.split('-')
    sats.append(sat)
    n=t_i//2
    dic = pd.read_pickle(path + sat + '.pkl')
    for i in range(0,len(xvc)):
        
        x,y,z=dic.get(n)
        dist = sat_dist(x, y, z, x0, y0, z0)
                        
        delay = prop_delay(dist*1000)
        g = cal_gain(h_sat, x, y, z, x0, y0, z0, 
                       fc=fc, rician_fading=rician_fading, atten=atten)
        p_rx = g * p_tx
        rate = cal_rate(p_rx, bw)
        
        if dist>th:
            sat=get_closest_sat_at_n(path, n, x0,y0,z0)
            sats.append(sat)
            dic = pd.read_pickle(path + sat + '.pkl')
            x,y,z=dic.get(n)
            dist = sat_dist(x, y, z, x0, y0, z0)
            delay = prop_delay(dist*1000)
            g = cal_gain(h_sat, x, y, z, x0, y0, z0, fc=fc, rician_fading=rician_fading, atten=atten)
            p_rx = g * p_tx
            rate = cal_rate(p_rx, bw)
            
        pve.append(dist)        
        rve.append(rate)
        dve.append(delay)
        n=n+1

    # append -t for handover count method
    sats = [s+'-t' for s in sats]
    
    return sats, rve, dve, pve


"""
build_table must be called before
"""
def get_closest_sat_at_n(path, n, x0,y0,z0):
    min_d=1e4
    min_sat=''
    dir_list = os.listdir(path)
    dir_list = [x[:-4] for x in dir_list]
    
    for sat in dir_list:
        dic = pd.read_pickle(path + sat + '.pkl')
        x,y,z=dic.get(n)
        dist=sat_dist(x, y, z, x0, y0, z0)
        if dist<min_d:
            min_d=dist
            min_sat=sat  
            
    return min_sat
    



"""
default is newyork at lat=40.730610, lon=-73.935242
california at lat=36.778259,lon=-119.417931
texas Lon=-100.000000, lat=31.000000
florida lon=-81.760254, lat= 27.994402
paris long=2.349014, lat=48.864716
aswan lat= 24.088938, long=32.899830
ottawa Lat=45.424721, long=-75.695000
"""
def calcPosFromLatLonRad(lat=40.730610,lon=-73.935242,radius=6378):
  
    phi   = (90-lat)*(math.pi/180);
    theta = (lon+180)*(math.pi/180);

    x = -(radius * math.sin(phi)*math.cos(theta));
    z = (radius * math.sin(phi)*math.sin(theta));
    y = (radius * math.cos(phi));
  
    return x,y,z;

"""
option = 'scatter' or 'plot' or 'step
QoS = 'delay' or 'rate'


legend location
'best'	        0
'upper right'	1
'upper left'	2
'lower left'	3
'lower right'	4
'right'     	5
'center left'	6
'center right'	7
'lower center'	8
'upper center'	9
'center'        10

"""
def plot_data(x, y, t_i, sln, nh, qos='delay', option='step', color='no'):
   
    colors = ['r', 'g', 'b', 'k', 'c', 'm', 'y']
    
    if qos == 'delay':
        tit = 'Delay'
        ylab = 'Propagation delay (msec)'
        y = [e*100 for e in y]#msec
    elif qos == 'rate':
        tit = 'Rate'
        ylab = 'Throuput (Mbps)'
        y = [e/1e6 for e in y] #Mbps
    
    if option == 'scatter':
        plt.scatter(x, y, color='green', label='graph method (%i handovers)' %nh)
        plt.legend()
    elif option == 'plot':
        if color == 'yes':
            c,sn=0,0
            for i in range(0, len(y), t_i):
                plt.plot(x[i:i+t_i], y[i:i+t_i], color=colors[c], label=sln[sn])
                sn=sn+1
                if c<len(colors)-1:
                    c=c+1
                else:
                    c=0
            plt.legend()
        else:
#            x = [e+(t_i/60/2) for e in x]
            plt.plot(x, y, color='green', label='graph method (%i handovers)' %nh)
            plt.xlim([0, 30])
            plt.legend()
    elif option == 'step':
        plt.step(x, y, color='green', where='post', label='graph method (%i handovers)' %nh)
    
    elif option == 'bar':
        x = [e+(t_i/60/2) for e in x]
        plt.bar(x[:-1], y[:-1])
        plt.xlim([0, 30])
    
    #plt.plot(x_l[:-1], delay_l[:-1], color='red', 
    #         label='elevation angle method (3 handovers)')  
    plt.axhline(st.mean(y), color = 'g', linestyle = '--')
    #plt.axhline(st.mean(delay_l[:-1]), color = 'r', linestyle = '--')

    plt.grid()
    #plt.legend(title='Parameter where:')
    plt.title(tit + '\n graph method: ' + str(nh) + ' handovers' + 
              '\n $t_i = \mu \pm \lambda=$' + str(t_i) + ' sec')
    plt.xlabel('Time (min)')
    if option=='plot' and color=='yes':
        plt.legend(bbox_to_anchor=(1.15, 1.0), loc="upper center", ncol=1)
    plt.ylabel(ylab)
#    plt.savefig('Demo/cont-rate-5-min.pdf', bbox_inches='tight')
    plt.show()
    
    
def plot_two(x, y, xe, ye, nh, nhe, t_i, qos='delay', save='no', fname=''):
    if qos == 'delay':
        tit = 'Delay'
        ylab = 'Propagation delay (msec)'
        y = [e*100 for e in y]#msec
        ye = [e*100 for e in ye]#msec

    elif qos == 'rate':
        tit = 'Data Rate'
        ylab = 'Data rate (Mbps)'
        y = [e/1e6 for e in y] #Mbps
        ye = [e/1e6 for e in ye]#msec

#'-', '--', '-.', ':'
    plt.plot(x, y, color='green', linestyle = '-', label='GM: %i handovers' %nh)
    plt.plot(xe, ye, color='red', linestyle = '-.', label='TH: %i handovers' %nhe) 
    
#    plt.axhline(st.mean(y), color = 'g', linestyle = '--')
#    plt.axhline(st.mean(ye), color = 'r', linestyle = '--')

    plt.legend( loc='lower left')
    plt.grid()
    #plt.legend(title='Parameter where:')
    plt.title('$t_i = \mu \pm \lambda=$' +str(t_i)+' sec ')
    plt.xlabel('Time (min)')
    plt.ylabel(ylab)
    if save=='yes':
        plt.savefig(fname, bbox_inches='tight')
    plt.show()
   

def pdf_cdf(data):
    # getting data of the histogram
    count, bins_count = np.histogram(data, bins='auto')      
    # finding the PDF of the histogram using count values
    pdf = count / sum(count)
    # using numpy np.cumsum to calculate the CDF
    # We can also find using the PDF values by looping and adding
    cdf = np.cumsum(pdf)
    # plotting PDF and CDF
#    plt.plot(bins_count[1:], pdf, color="red", label="PDF")
    plt.plot(bins_count[1:11], cdf[:10], label="CDF")
    plt.legend()
    return count, bins_count, pdf, cdf

def percentile(data, percentage):
#    np.percentile(rve, 10)
    return stats.scoreatpercentile(data, percentage)

def percentile_revevrse(data, value):
    return stats.percentileofscore(data,value)



def cdf_plot(path, p_wt, nth_per, include_seconds='no', save_plot='no'):
    hf=[1,2,3,4,5,6,7,8,9,10]
    seconds = [5, 15, 30, 45] #be careful, freezez computer
    rate_10nth=[]
#    n_sat, n_orb = 72, 22 #22sat/orbit; 1584 satellites
    wd, wr=0.5, 0.5
    freq =11.9e9 #ku band see reference above in comments
    ric_fad=100 #20 db
    atten = 0.05 #about 0.05 db/km for ka band in good weather
    bw=10e6 #10 MHz, total is 240 MHz
    p_tx=10 #10 db
    x0,y0,z0= calcPosFromLatLonRad(lat=45.424721, lon=-75.695000) #OTTAWA
    h_sat=550 #km
    
    for nmin in hf:
        t_i = nmin*60 # 5min (t_i is multiple of T)
        T = (30*60)+t_i # 30min, t_i to offset dummy node
        db = pd.read_pickle(p_wt+'wt-'+str(nmin)+'-min.pkl')
        g = build_graph(db,wd,wr) #th=
        sln = solve_graph(g)[1:] #ignore first dummy node, last t is not in time range
        rvc, dvc = pull_data(path, db, sln[:-1], T, t_i,h_sat, x0, y0, z0, p_tx, bw, freq, ric_fad, atten)
        xvc = list(np.linspace(0,(T-t_i)/60, len(rvc)))
        rate_10nth.append(stats.scoreatpercentile(rvc,nth_per))
        
        r_cdf=[]
        y_cdf=[]
        for i in range(0,nth_per+1):
            r_cdf.append(stats.scoreatpercentile(rvc,i))
            y_cdf.append(i)
        r_cdf = [i/1e6 for i in r_cdf]   
        plt.plot(r_cdf, y_cdf, label='GM $t_i$=%i min' %nmin)
    
    
#    if include_seconds == 'yes':
#        for sec in seconds:
#            t_i = sec # 5min (t_i is multiple of T)
#            T = (30*60)+t_i # 30min, t_i to offset dummy node
#            db = pd.read_pickle(p_wt+'wt-'+str(sec)+'-sec.pkl')
#            g = build_graph(db,wd,wr) #th=
#            sln = solve_graph(g)[1:] #ignore first dummy node, last t is not in time range
#            rvc, dvc = pull_data(path, db, sln[:-1], T, t_i,h_sat, x0, y0, z0, p_tx, bw, freq, ric_fad, atten)
#            xvc = list(np.linspace(0,(T-t_i)/60, len(rvc)))
#            rate_10nth.append(stats.scoreatpercentile(rvc,nth_per))
#            
#            r_cdf=[]
#            y_cdf=[]
#            for i in range(0,nth_per+1):
#                r_cdf.append(stats.scoreatpercentile(rvc,i))
#                y_cdf.append(i)
#            r_cdf = [i/1e6 for i in r_cdf]   
#            plt.plot(r_cdf, y_cdf, label='gm $t_i$=%i sec' %sec)
     
    #elevation method
    th=3200 #km, approx 10 degree elevation angle
    fs=sln[0]
    sats, rve, dve, pve = elev_method(path, fs, th, xvc, T, t_i, h_sat, x0, y0, z0, p_tx, bw, freq, ric_fad, atten)
   
    r_cdf=[]
    y_cdf=[]
    for i in range(0,nth_per+1):
        r_cdf.append(stats.scoreatpercentile(rvc,i))
        y_cdf.append(i)
    r_cdf = [i/1e6 for i in r_cdf]   
    plt.plot(r_cdf, y_cdf, color='red', linestyle='--', label= 'TH '+r'$\alpha=10^o$')
    
    #Hybrid method
#    nmin=5
#    t_i = nmin*60 # 5min (t_i is multiple of T)
#    T = (30*60)+t_i # 30min, t_i to offset dummy node
#    db = pd.read_pickle(p_wt+'wt-'+str(nmin)+'-min.pkl')
#    g = build_graph(db,wd,wr) #th=
#    sln = solve_graph(g)[1:] #ignore first dummy node, last t is not in time range
#    rvc, dvc = pull_data(path, db, sln[:-1], T, t_i,h_sat, x0, y0, z0, p_tx, bw, freq, ric_fad, atten)
#    xvc = list(np.linspace(0,(T-t_i)/60, len(rvc)))
#    both_methods_sln = sln+sats
#    dbh = load_table_hybrid(path, T, t_i, both_methods_sln)
#    dbh = build_table(dbh, h_sat, x0, y0, z0, p_tx, bw, freq, ric_fad, atten)
#    gh = build_graph(dbh,wd,wr) #th=
#    slnh = solve_graph(gh)[1:] 
#    rvh, dvh = pull_data(path, dbh, slnh[:-1], T, t_i,h_sat, x0, y0, z0, p_tx, bw, freq, ric_fad, atten)
#    r_cdf=[]
#    y_cdf=[]
#    for i in range(0,nth_per+1):
#        r_cdf.append(stats.scoreatpercentile(rvc,i))
#        y_cdf.append(i)
#    r_cdf = [i/1e6 for i in r_cdf]   
#    plt.plot(r_cdf, y_cdf, color='black', linestyle='loosely dashed', label='hd: th & $t_5$ ')
    
#    plt.legend()
    plt.legend(bbox_to_anchor=(1.17, 1.0), loc="upper center", ncol=1)
    plt.grid()
    #plt.legend(title='Parameter where:')
    plt.title('Data Rate CDF')
    plt.xlabel('Data rate (Mbps)')
    plt.ylabel('CDF (%)')
    if save_plot=='yes':
        plt.savefig('Demo/cdf-100.pdf', bbox_inches='tight')
    return
    
        
"""
def working_table(path, T, t_i):
    dir_list = os.listdir(path)
    dir_list = [x[:-4] for x in dir_list]
    ts=[]
    for i in range(0,(T//t_i)+1):
        ts.append('t'+str(i))
    
    db = pd.DataFrame(columns = ts, index=dir_list)  
    
    for index, row in db.iterrows():
        dic = pd.read_pickle(path + str(index) + '.pkl')
        db.at[index, 't0']=dic.get(0)
        n=t_i
        for t in range(1,len(ts)):
            db.at[index, 't'+str(t)]=dic.get(n)
            n=n+t_i
    return db
"""

"""
combine pickle files to one directory 0-1200 and 1201-2400 to 0-2400
"""
def combine_pickle(path1='', path2='', path_to=''):
    path1='SavedPositions-0000-1200/'
    path2='SavedPositions-1201-2400/'
    path_to='SavedPositions-0000-2400/'
    dir_list = os.listdir(path1)
    my_dict_final = {}  # Create an empty dictionary
    for fname in dir_list:
        with open(path1 + fname, 'rb') as f:
            my_dict_final.update(pickle.load(f))     
        with open(path2 + fname, 'rb') as f:
            my_dict_final.update(pickle.load(f))  
        with open(path_to + fname, 'wb') as handle:
            pickle.dump(my_dict_final, handle, protocol=pickle.HIGHEST_PROTOCOL)
    return


"""
14, 15 check last column special case

"""
def mean_stats(p_wt, t_i_s=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]):
    cols=['nh', 'delay-bar', 'delay-cont', 'rate-bar', 'rate-cont']
    stats = pd.DataFrame(columns=cols, index=t_i_s) 
    path = 'SavedPositions-0000-2400/'
#    n_sat, n_orb = 72, 22 #22sat/orbit; 1584 satellites
    wd, wr=0.5, 0.5
    freq =11.9e9 #see reference above in comments
    ric_fad=100 #20 db
    atten = 0.05 #about 0.05 db/km for ka band in good weather
    bw=10e6 #10 MHz, total is 240 MHz
    p_tx=10 #10 db
    x0,y0,z0= calcPosFromLatLonRad(lat=45.424721, lon=-75.695000) #OTTAWA
    h_sat=550 #km
    for t in t_i_s:
        t_i = t*60 
        T = (30*60)+t_i #30min, t_i to offset dummy node
        db = pd.read_pickle(p_wt+'wt-'+str(t)+'-min.pkl')

        g = build_graph(db,wd,wr)
        sln = solve_graph(g)[1:] #ignore first dummy node, last t is not in time range
        
        rv, dv, pv = get_data(db,sln[:-1]) 
        rvc, dvc = pull_data(path, db, sln[:-1], T, t_i,h_sat, x0, y0, z0, p_tx, bw, freq, ric_fad, atten)
        
        stats.loc[t]['delay-bar']=st.mean(dv)*100
        stats.loc[t]['rate-bar']=st.mean(rv)/1e6
        stats.loc[t]['delay-cont']=st.mean(dvc)*100
        stats.loc[t]['rate-cont']=st.mean(rvc)/1e6
        stats.loc[t]['nh']=len(sln)-1-1
    return stats
    
    
    
    