import numpy as np

def simulate(cue, rwd,  dt=2, T=12*1000, 
             tau1=52, sigma1=1.4, tauE1=100, tauS1=50, 
             tau2=100, sigma2=.025, tauE2=100, tauS2=50*150, 
             p=1.5, scale_cue=1, scale_rwd=2):
    
    # initialize 
    nt = T // dt
    tempWE1 = np.exp(-np.arange(0, nt)*dt/tauE1) * dt/tauE1  
    tempWS1 = np.exp(-np.arange(0, nt)*dt/tauS1) * dt/tauS1 
    tempWE2 = np.exp(-np.arange(0, nt)*dt/tauE2) * dt/tauE2 
    tempWS2 = np.exp(-np.arange(0, nt)*dt/tauS2) * dt/tauS2
    drive1, d1, s1, f1, r1 = np.zeros((4, nt)), np.zeros((4, nt)), np.zeros((4, nt)), np.zeros((4, nt)), np.zeros((4, nt))
    drive2, d2, s2, f2, r2 = np.zeros(nt), np.zeros(nt), np.zeros(nt), np.zeros(nt), np.zeros(nt)

    # simulate
    for i in range(1, nt):
        # sensory cue layer
        drive1[:, i] = (scale_cue*cue[:, i])**p
        d1[:, i] = np.sum(drive1[:, :i] * tempWE1[i-1::-1], axis=1)
        s1[:, i] = np.sum(np.abs(d1[:, :i]) * tempWS1[i-1::-1])
        f1[:, i] = d1[:, i] / (s1[:, i] + sigma1**p)
        r1[:, i] = r1[:, i-1] + dt/tau1*(-r1[:, i-1]+f1[:, i]) 
        
        # reward layer
        drive2[i] = r1[0, i]**p + (scale_rwd*rwd[i])**p
        d2[i] = np.sum(drive2[:i] * tempWE2[i-1::-1])
        s2[i] = np.sum(np.abs(d2[:i]) * tempWS2[i-1::-1]) 
        f2[i] = d2[i] / (s2[i] + sigma2**p)
        r2[i] = r2[i-1] + dt/tau2*(-r2[i-1]+f2[i]) 

    return r2