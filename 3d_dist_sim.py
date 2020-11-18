import numpy as np


def dist_calc(x_arr, x_val):
    #x_arr is array of positions, x_val is just one position
    dist = np.sqrt((x_arr[:,0]-x_val[0])**2 + (x_arr[:,1]-x_val[1])**2 + (x_arr[:,2]-x_val[2])**2)
    return dist

def dist_calc_2d(x_arr, x_val):
    #x_arr is array of positions, x_val is just one position
    dist = np.sqrt((x_arr[:,0]-x_val[0])**2 + (x_arr[:,1]-x_val[1])**2)
    return dist

def run_sim(N, sigma=1):
    pos_vals = np.random.normal(0,sigma,(N,3))

    dist_vals_3d = []
    dist_vals_2d = []
    for num in range(N):
        dists_3d = dist_calc(pos_vals[num+1:,:], pos_vals[num,:])
        dist_vals_3d.append(dists_3d)    
        
        dists_2d = dist_calc_2d(pos_vals[num+1:,0:2], pos_vals[num,0:2])
        dist_vals_2d.append(dists_2d)
        
    full_dist_3d = np.concatenate(dist_vals_3d)
    full_dist_2d = np.concatenate(dist_vals_2d)

    mean2d = np.mean(full_dist_2d)
    mean3d = np.mean(full_dist_3d)
    
    print(f'mean, median 2d separation: {mean2d}, {np.median(full_dist_2d)}')
    print(f'mean, median 3d separation: {mean3d}, {np.median(full_dist_3d)}')
    print(f'mean3d/mean2d = {mean3d/mean2d}')

run_sim(10000)
