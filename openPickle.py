# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 19:30:55 2022

@author: Mohamed Hozayen

22 orbit
72 satellite per orbit
22*72 = 1584 files
each file contains 20 min (1201 rows) of xyz coordinates

satelltie altitude is 550 km
earth radius is 6378 km

"1 + orbit_num + sat_num".pkl
10101 - orbit 1 sat 1
10512 - orbit 5 sat 12
11220 - orbit 12 sat 20

https://www.techtarget.com/searchmobilecomputing/definition/downlink-and-uplink
Frequency Band	  Downlink	       Uplink
     C	       3,700-4,200 MHz	 5,925-6,425 MHz
     Ku	        11.7-12.2 GHz	  14.0-14.5 GHz
     Ka	        17.7-21.2 GHz	  27.5-31.0 GHz

Feedback WG
check with confidence legacy handover in teresstrial (surverys)
cdf of throughput 10th percentile
    outage situation, 5th percentile
    eliminate worst case

more figures: 
    cdf: handover frequnce vs rate -- saturation point
    ti size -> handover frequency
    weights, need to rerun build table

maybe put graph in title

multiple graph path, use best one
optimize reduce handover (df factor)
hybrid solution, graph and threshold

total delay = propagation + processing + queuing + transmision
handovers are based on link and loading
    capacity flavour, load balancing, satellite load distribution
    user density on the ground

"""

#import pickle
#import csv
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import pandas as pd
import eqn
import statistics as st
from scipy import stats



path = 'SavedPositions-0000-2400/'
p_wt = 'Dataframe-Tables/'

n_sat, n_orb = 72, 22 #22sat/orbit; 1584 satellites
wd, wr=0.5, 0.5
freq =11.9e9 #ku band see reference above in comments
ric_fad=100 #20 db
atten = 0.05 #about 0.05 db/km for ka band in good weather
bw=10e6 #10 MHz, total is 240 MHz
p_tx=10 #10 db
x0,y0,z0= eqn.calcPosFromLatLonRad(lat=45.424721, lon=-75.695000) #OTTAWA
h_sat=550 #km

# (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)-min
# (5, 15, 30, 45)-sec change value of t_i directly
nmin=4 # multiples of 30 specially when u plot
t_i = nmin*60 # 5min (t_i is multiple of T)
T = (30*60)+t_i # 30min, t_i to offset dummy node

#wd, wr=0.3, 0.7

"""
    GRAPH METHOD
"""
db = pd.read_pickle(p_wt+'wt-'+str(nmin)+'-min.pkl')
#db = eqn.load_table(path, T, t_i)
#db = eqn.build_table(db, h_sat, x0, y0, z0, p_tx, bw, freq, ric_fad, atten)
g = eqn.build_graph(db,wd,wr, df=0, df_option='yes') #th=
sln = eqn.solve_graph(g)[1:] #ignore first dummy node, last t is not in time range
#sln = eqn.solve_graph(g, mode=2)

#step plot 
#rv, dv, pv = eqn.get_data(db,sln) 
#xv = list(np.linspace(0, (T-t_i)/60, len(rv)))
#eqn.plot_data(xv, dv, t_i, sln, nh=len(sln)-1-1, qos='delay', option='plot')
#eqn.plot_data(xv, rv, t_i, sln, nh=len(sln)-1-1, qos='rate', option='plot')

#continuous plot
#####################   attention to 14,15 last t is deleted as not in dataset range
#                       min last solution sln[:-1]
rvc, dvc = eqn.pull_data(path, db, sln[:-1], T, t_i,h_sat, x0, y0, z0, p_tx, bw, freq, ric_fad, atten)
xvc = list(np.linspace(0,(T-t_i)/60, len(rvc)))
#eqn.plot_data(xvc, dvc, t_i, sln, nh=eqn.count_handovers(sln[:-1]), 
#              qos='delay', option='plot', color='no')
#eqn.plot_data(xvc, rvc, t_i, sln, nh=eqn.count_handovers(sln[:-1]), 
#              qos='rate', option='plot', color='no')


"""
    ELEVATION ANGLE METHOD METHOD
    
"""
#dbe = pd.read_pickle(p_wt+'wt-1-sec.pkl')
th=3200 #km, approx 10 degree elevation angle
fs=sln[0]
sats, rve, dve, pve = eqn.elev_method(path, fs, th, xvc, T, t_i, h_sat, x0, y0, z0, p_tx, bw, freq, ric_fad, atten)
#eqn.plot_data(xvc, dve, t_i, sln, nh=len(sats)-1, qos='delay', option='plot', color='no')
#eqn.plot_data(xvc, rve, t_i, sln, nh=len(sats)-1, qos='rate', option='plot', color='no')


"""
Comparing two methods
graph      xvc, dvc, rvc
elevation       dve, rve
"""
fnd='Demo/compare-delay-' +str(nmin) + '-min.pdf'
fnr='Demo/compare-rate-' + str(nmin) + '-min.pdf'
#fnd='Demo/compare-delay--'+ 'wd'+str(wd) + '--wr'+str(wr) +'-5min-.pdf'
#fnr='Demo/compare-rate--'+ 'wd'+str(wd) + '--wr'+str(wr) +'-5min.pdf'
eqn.plot_two(xvc, dvc, xvc, dve, nh=eqn.count_handovers(sln[:-1]), 
             nhe=eqn.count_handovers(sats), t_i=t_i, qos='delay', save='yes', fname=fnd)
eqn.plot_two(xvc, rvc, xvc, rve, nh=eqn.count_handovers(sln[:-1]), 
             nhe=eqn.count_handovers(sats), t_i=t_i, qos='rate', save='yes', fname=fnr)


"""
plot frequency of handovers vs rate at 10th-percentile  
"""
eqn.cdf_plot(path, p_wt, 100, save_plot='yes')
stats.scoreatpercentile(rve,20)


"""
Hydrid method

build a graph og sln+sats
"""
both_methods_sln = sln+sats
nmin=3 # multiples of 30 specially when u plot
t_i = nmin*60 # 5min (t_i is multiple of T)
T = (30*60)+t_i # 30min, t_i to offset dummy node
dbh = eqn.load_table_hybrid(path, T, t_i, both_methods_sln)
dbh = eqn.build_table(dbh, h_sat, x0, y0, z0, p_tx, bw, freq, ric_fad, atten)
gh = eqn.build_graph(dbh,wd,wr,th=2000) #th=
slnh = eqn.solve_graph(gh)[1:] #ignore first dummy node, last t is not in time range
slnhs = eqn.solve_graph(gh, mode=2) #ignore first dummy node, last t is not in time range
#slnh = eqn.filter_handovers(slnh)
rvh, dvh = eqn.pull_data(path, dbh, slnh[:-1], T, t_i,h_sat, x0, y0, z0, p_tx, bw, freq, ric_fad, atten)
xvh = list(np.linspace(0,(T-t_i)/60, len(rvh)))
eqn.plot_data(xvh, rvh, t_i, slnh, nh=eqn.count_handovers(slnh), qos='rate', option='plot', color='no')
eqn.plot_data(xvh, rvh, t_i, slnh, eqn.count_handovers(slnh), qos='rate', option='plot', color='yes')


eqn.count_handovers(slnh)
eqn.filter_handovers(slnh)






"""
Testing random
"""
n1, n2 = 0+t_i//2, 30*60
rvh_1, dvh_1 = eqn.pull_sat(path, db, sats[1], n1, n2, T, t_i,h_sat, x0, y0, z0, p_tx, bw, freq, ric_fad, atten)
rvh_0, dvh_0 = eqn.pull_sat(path, db, sats[0], n1, n2, T, t_i,h_sat, x0, y0, z0, p_tx, bw, freq, ric_fad, atten)
xvh_0 = list(np.linspace(0, n2//60, len(rvh_0)))


rvh_0 = [i/1e6 for i in rvh_0]
rvh_1 = [i/1e6 for i in rvh_1]

plt.plot(xvh_0, rvh_0)
#plt.plot(xvh_0, rvh_1) 
plt.grid()
#plt.legend(title='Parameter where:')
#plt.legend()
plt.title('Rate')
plt.xlabel('Time (min)')
plt.ylabel('Throughput (Mbps)')
#plt.savefig(paths+"/delay-" + str(interval_duration) + "sec-"+str(df)+"df.pdf")
plt.show()


for s in sats:
    n1, n2 = 0+t_i//2, 30*60
    rvh, dvh = eqn.pull_sat(path, db, s, n1, n2, T, t_i,h_sat, x0, y0, z0, p_tx, bw, freq, ric_fad, atten)
    xvh = list(np.linspace(0, n2//60, len(rvh)))
    rvh = [i/1e6 for i in rvh]
    plt.plot(xvh, rvh) 
    
plt.grid()
#plt.legend(title='Parameter where:')
#plt.legend()
plt.title('Rate')
plt.xlabel('Time (min)')
plt.ylabel('Throughput (Mbps)')
#plt.savefig(paths+"/delay-" + str(interval_duration) + "sec-"+str(df)+"df.pdf")
plt.show()

for s in (sln+sats):
    n1, n2 = 0+t_i//2, 30*60
    rvh, dvh = eqn.pull_sat(path, db, s, n1, n2, T, t_i,h_sat, x0, y0, z0, p_tx, bw, freq, ric_fad, atten)
    xvh = list(np.linspace(0, n2//60, len(rvh)))
    rvh = [i/1e6 for i in rvh]
    plt.plot(xvh, rvh) 

plt.grid()
#plt.legend(title='Parameter where:')
#plt.legend()
plt.title('Rate')
plt.xlabel('Time (min)')
plt.ylabel('Throughput (Mbps)')
#plt.savefig(paths+"/delay-" + str(interval_duration) + "sec-"+str(df)+"df.pdf")
plt.show()

  

#print(stats.scoreatpercentile(rvc,10)) #same as np.percnetile
#print(stats.percentileofscore(rve,11567960.8290049))
#
#r_cdf=[]
#x_cdf=[]
#for i in range(0,101):
#    r_cdf.append(stats.scoreatpercentile(rve,i))
#    x_cdf.append(i)
#x_cdf = [i/100 for i in x_cdf]    
#
#eqn.pdf_cdf(rve)
#count, bins_count = np.histogram(rve, bins='auto')      
## finding the PDF of the histogram using count values
#pdf = count / sum(count)
## using numpy np.cumsum to calculate the CDF
## We can also find using the PDF values by looping and adding
#cdf = np.cumsum(pdf)
## plotting PDF and CDF
#plt.plot(bins_count[1:], pdf, label="PDF")
#plt.plot(bins_count[1:], cdf, label="CDF")
#plt.legend()
#
##stats=eqn.mean_stats(p_wt,t_i_s=[1,2,3,4,5,6,7,8,9,10]).to_csv(p_wt+'stats.csv')
#stats = pd.read_csv(p_wt+'stats.csv')
#stats.rename(columns={stats.columns[0]: "t_i" }, inplace = True)
#ax=stats.plot(x='t_i', y=['delay-bar', 'delay-cont'])
#ax=stats.plot(x='t_i', y=['rate-bar', 'rate-cont'])
#ax=stats.plot.bar(x='t_i', y=['nh'])
#fig = ax.get_figure()
#fig.savefig('Demo/number-of-handovers.pdf')

#stats_temp = eqn.call_stats(t_i_s=['0.5', '1', '2', '3', '5', '6', '10'])
#stats_temp = eqn.call_stats(t_i_s=['0.5', '1', '2', '3', '5', '6', '10'])
#dic = pd.read_pickle('SavedPositions-0000-2400/12108.pkl')
#eqn.search_dictionary(dic, 811.2766800990366, 5473.36587611792, 4169.318888508697)

"""
Next steps:
    capacity
    total delay = propagation + processing + queuing + transmision

np.percentile(rve, 10)
stats.scoreatpercentile(rve,50) #same as np.percnetile
stats.percentileofscore(rve,11567960.8290049)


s0,t0=sln[0].split('-')
sln.insert(0,s0+'-t0')
sln.pop()
sl,tl=sln[-1].split('-')
sln.append(sl+'-t'+str(len(sln)))


ax=stats.plot(x='t_i', y=['delay-bar', 'delay-cont'])
fig = ax.get_figure()
fig.savefig('figure.pdf')


#    or dataframe.to_pickle("./dummy.pkl")  

#with open(p_wt + 'wt-'+str(nmin)+'-sec.pkl', 'wb') as handle:
#    pickle.dump(db, handle, protocol=pickle.HIGHEST_PROTOCOL)


    # available tables: 1:10, 15 min || 5,10,15,30,45 seconds
#    or dataframe.to_pickle("./dummy.pkl")  
    
    with open(p_wt + 'wt-'+str(nmin)+'-sec.pkl', 'wb') as handle:
        pickle.dump(db, handle, protocol=pickle.HIGHEST_PROTOCOL)


#db.to_csv(p_wt+'working-table-'+str(t_i)+'-sec.csv')

plt.plot(x, dvc, color='green', 
         label='graph method (%i handovers)' %(len(rv)-1))
#plt.plot(x_l[:-1], delay_l[:-1], color='red', 
#         label='elevation angle method (3 handovers)')  
plt.axhline(st.mean(dvc), color = 'g', linestyle = '--')
#plt.axhline(st.mean(delay_l[:-1]), color = 'r', linestyle = '--')
plt.grid()
#plt.legend(title='Parameter where:')
plt.legend()
plt.title('Delay QoS' + " - $\lambda=$" + str(t_i//60) + " min")
plt.xlabel('Duration (min)')
plt.ylabel('Propagation delay (msec)')
#plt.savefig(paths+"/delay-" + str(interval_duration) + "sec-"+str(df)+"df.pdf")
plt.show()

plt.plot(x, rvc, color='green', 
         label='graph method (%i handovers)' %(len(rv)-1))
#plt.plot(x_r[:-1], rate_l[:-1], color='red', 
#         label='elevation angle method (3 handovers)')  
plt.axhline(st.mean(rvc), color = 'g', linestyle = '--')
#plt.axhline(st.mean(rate_l[:-1]), color = 'r', linestyle = '--')
plt.grid()
#plt.legend(title='Parameter where:')
plt.legend()
plt.title('Rate QoS' + " - $\lambda=$" + str(t_i//60) + " sec")
plt.xlabel('Duration (min)')
plt.ylabel('Throughput (Mbps)')
#plt.savefig(paths+"/rate-" + str(interval_duration) + "sec-"+str(df)+"df.pdf")
plt.show()


#s,t=sln[1].split('-') returns s=sat_num, t=tx -> db.loc[s][t]

for rowIndex, row in db.iterrows(): #iterate over rows
    for columnIndex, value in row.items():
        print(value, end="\t")
    print()

#filename = 'SavedPositions-0000-1200/11220.pkl'
#df = pd.read_pickle('SavedPositions-0000-1200/10101.pkl')

X,Y,Z=[],[],[]
for k in df:
    c=df.get(k)
    X.append(c[0])
    Y.append(c[1])
    Z.append(c[2])

# Plot X,Y,Z
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_trisurf(X[:], Y[:], Z[:], color='white', edgecolors='grey', alpha=0.5)
ax.scatter(X[:], Y[:], Z[:], c='red')
plt.show()


n1=0
n2=1
earth_radius=6371
for i in range(0, len(X)-1):
    dxs=(X[n1]-X[n2])**2
    dys=(Y[n1]-Y[n2])**2
    hs=(Z[n1]-Z[n2])**2
    print((hs+dxs+dys)**(1/2))
    n1+=1
    n2+=1

# 6928.139 distance from origin
# 556.139 distance from earth
# 604.4024 distance between two adjacent satellite

dxs=(X[n1]-0)**2
dys=(Y[n1]-0)**2
hs=(Z[n1]-earth_radius)**2
print((hs+dxs+dys)**(1/2))

for name, age in df.items():  # for name, age in dictionary.iteritems():  (for Python 2.x)
    if age == (6309.727939240258, -2235.4184583579463, 1785.8783576726073):
        print(name)

for i in sln:
    s,t=i.split('-')
    x,y,z,dist,rate,delay,rate_norm,delay_norm=db.loc[s][t]
    x0,y0,z0=6378*0.707107,6378*0.707107,0
    d = eqn.sat_dist(x, y, z, x0, y0, z0)
    dd = eqn.prop_delay(d*1000)
        
    g = eqn.cal_gain(h_sat, x, y, z, x0, y0, z0, 
                         fc=freq, rician_fading=ric_fad, atten=0.05)
    p_rx = g * p_tx
    rate = eqn.cal_rate(p_rx, bw)
    print(x,y,z,d,rate,dd)
    
"""



