import numpy as np

c=2.0

# hpmts in numpy array of form [[x1,y1,z1,t1],[x2,y2,z2,t2], [x3,y3,z3,t3], [x4,y4,z4,t4]] 
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

# data for c=2, event at (0,0,0,2)

hpmts = np.array([[ 0.0, -1.0, 0.0, 1/4+6],
                  [ 1.0,  0.0, 0.0, 1/4+6],
                  [-1.0,  0.0, 0.0, 1/4+6],
                  [ 0.0,  3.0, 4.0, 5/4+6]])

print(quadrangulate(hpmts, np.array([0.0,0.0,1.0,1/4+6])))
