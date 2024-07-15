"""This quadfitter takes a collection of PMT hits from an event and attempts to 
reconstruct the event vertex by calculating the predicted location of the event for
sets of 4 PMTs. See Report on the Quadfitter by Ian Coulter for a theoretical
overview."""

import ROOT
import sys
import numpy as np
import random
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

filename = sys.argv[1]
f = ROOT.TFile.Open(filename)
t = f.Get("output")
m = f.Get("meta")
# speed of light in water in mm/ns
c = 225.0
# number of sets of 4 PMTs to iterate over per event
n = 3000

# takes set of 4 hit PMTs and reconstructs the position based on the hit times and locations of the PMTs
# follows math and variable naming scheme of Ian Coulter's quadfitter report
# hpmts in numpy array of form [[x1,y1,z1,t1], [x2,y2,z2,t2], [x3,y3,z3,t3], [x4,y4,z4,t4]] 
def quadrangulate(hpmts, check):
	
    # compute matricies M, N, K, G as defined in quadfitter report
    dif = np.array([hpmts[1] - hpmts[0], hpmts[2] - hpmts[0], hpmts[3] - hpmts[0]])
    M = dif[:, :3]
    N = dif[:, -1]
    K = 0.5*np.array([np.dot(p[:3], p[:3]) - np.dot(hpmts[0][:3], hpmts[0][:3]) 
                      - c*c*(np.dot(p[-1], p[-1]) - np.dot(hpmts[0][-1], hpmts[0][-1])) 
                      for p in (hpmts[1], hpmts[2], hpmts[3])])

    # discards sets of pmts that have symmetries that result in uninvertable matricies
    if(np.linalg.det(np.array(M)) == 0): return [10000,0,0,0]
    
    # compute matricies G, H as defined in quadfitter report
    M_inv = np.linalg.inv(M)
    G = np.dot(M_inv,K)
    H = np.dot(M_inv,N)

  	# compute coefficients of quadratic as defined in quadfitter report
    c0 = np.dot(hpmts[0, :3]-G,hpmts[0, :3]-G) - np.dot(c*hpmts[0, -1], c*hpmts[0, -1])
    c1 = -2*c**2*(np.dot(hpmts[0, :3]-G,H)-hpmts[0, -1])
    c2 = c**2*(c**2*np.dot(H,H)-1)

    # calculate roots of quadratic and discard the non-real / unphysical root
    times = np.roots([c2,c1,c0])
    positions = [G+c**2*time*H for time in times]

    if len(positions) == 2.0:
      if abs((np.sqrt((check[0]-positions[0][0])**2 + (check[1]-positions[0][1])**2 + (check[2]-positions[0][2])**2))/c - check[3] + times[0]) < abs((np.sqrt((check[0]-positions[1][0])**2 + (check[1]-positions[1][1])**2 + (check[2]-positions[1][2])**2))/c - check[3] + times[1]) and np.any(abs(positions[0])) > 900.0:
        position = np.append(positions[0], times[0])
      else: position = np.append(positions[1], times[1])
    else: position = np.append(positions[0], times[0])

    # discard sets of 4 hits that yeild wildly off positions
    if np.any(abs(position) > 900.0): return [10000,0,0,0]
	
    return position

# takes in array of hit pmts and selects a set of five indicies in range
def get_hpmt_indicies(nhits):
    return random.sample(range(nhits), 5)

def get_best_fit():
    reconstructed_x = []
    reconstructed_y = []
    reconstructed_z = []
	    
    # loop through all simulated events
    for ev in range(t.GetEntries()):
        x_positions = []
        y_positions = []
        z_positions = []
        times = []

        t.GetEntry(ev)
        m.GetEntry(ev)
            
        pmtids = list(t.hitPMTID)
        pmttimes = list(t.hitPMTTime)
		   
        # create array of all pmt ids and hittimes over 7ms
        times = list(pmttime for pmttime in pmttimes if pmttime < 7.0)
		    
        ids = list(pmtid for pmtid, pmttime in zip(pmtids, pmttimes) if pmttime < 7.0)

        # only attempt to quadrangulate on events with at least 10 hits (less will not be accurate enough)
        if len(ids) >= 10:
            i=0
            # quadrangulate n times for each event
            while i < n:
                pmt_hits = get_hpmt_indicies(len(ids)) # indexes of PMTs to quadrangulate over
                pmt_positions = []
                    
                # get the data for each hit pmt
                for hit in pmt_hits[:4]:
                    pmt_positions.append([list(m.pmtX)[ids[hit]], list(m.pmtY)[ids[hit]], list(m.pmtZ)[ids[hit]], times[hit]])

                hpmts = np.array(pmt_positions)
                check = np.array([list(m.pmtX)[ids[3]], list(m.pmtY)[ids[3]], list(m.pmtZ)[ids[3]], times[3]])
                    
                event_position = quadrangulate(hpmts, check)	
                        
                # add the reconstructed positions to the event cloud
                if event_position[0] != 10000 and all(np.isreal(k) for k in event_position):
                    x_positions.append(event_position[0])
                    y_positions.append(event_position[1])
                    z_positions.append(event_position[2])
                    times.append(event_position[3])
                    i+=1
        
        # average the reconstructed positions across pmt choices to find the reconstructed position for one event
        median_x = np.median(x_positions)
        median_y = np.median(y_positions)
        median_z = np.median(z_positions)
        average_x = np.mean(x_positions)
        average_y = np.mean(y_positions)
        average_z = np.mean(z_positions)

        reconstructed_x.append(median_x)
        reconstructed_y.append(median_y)
        reconstructed_z.append(median_z)
        
        #print("")
        #print("MEDIAN LOCATION: (" + str(median_x) + ", " + str(median_y) + ", " + str(median_z) +")")
        #print("AVERAGE LOCATION: (" + str(average_x) + ", " + str(average_y) + ", " + str(average_z) +")")

    avg_rec_x = np.median(reconstructed_x)
    avg_rec_y = np.median(reconstructed_y)
    avg_rec_z = np.median(reconstructed_z)

    print("MEDIAN LOCATION: (" + str(avg_rec_x) + ", " + str(avg_rec_y) + ", " + str(avg_rec_z) + ")")
    
    plt.scatter(reconstructed_x, reconstructed_y, c="blue", label="Reconstruction Cloud", s=1)
    plt.scatter(median_x, median_y, c="red", label="Median Predicted Position", s=5)
    plt.scatter(150,0,c="orange", label="Simulated Position", s=5)
    plt.legend()
    plt.title("Position in XY-Plane for Event Located at (150,0,0)")
    plt.xlim(-1000,1000)
    plt.ylim(-1000,1000)
    plt.savefig("graphs/scatterplot_150_0_0.pdf")
    plt.clf()

    plt.hist(reconstructed_x, bins=30)
    plt.xlabel("X Coordinate of Reconstructed Event")
    plt.ylabel("Number of Events")
    plt.title("Distribution of X Coordinates for Event Located at (150,0,0)")
    plt.savefig("graphs/histogram_150_0_0.pdf")

get_best_fit()
