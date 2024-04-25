from b_ref import romberg, trap, golden, ab, scalar, gradient, A_value
import numpy as np 
import matplotlib.pyplot as plt

#everything basically as in b) just another function to minimize

r1, nhalo1 = np.load("radius1.npy"), np.load("nhalo1.npy")
r2, nhalo2 = np.load("radius2.npy"), np.load("nhalo2.npy")
r3, nhalo3 = np.load("radius3.npy"), np.load("nhalo3.npy")
r4, nhalo4 = np.load("radius4.npy"), np.load("nhalo4.npy")
r5, nhalo5 = np.load("radius5.npy"), np.load("nhalo5.npy")
r_array = np.array([r1,r2,r3,r4,r5],dtype=object)
halo_array = np.array([nhalo1,nhalo2,nhalo3,nhalo4,nhalo5])
n_sat_c_array = np.array([len(x) for x in r_array])#gesamtzahl
max_x_array = np.array([max(x) for x in r_array])
n_sat_array = n_sat_c_array/halo_array
N_array_array = np.array([np.load(f"n_array{j}.npy") for j in range(5)])

def like(a,b,c):
    """x_max_glob"""
    A = A_value(a,b,c,x_max_glob)
    tilde_array = np.zeros(74) #number of bins
    
    def integrand(x):
        return x**(a-1)*np.exp(-(x/b)**c)
    
    #print("now calculating tilde")
    for i in range(75-1):
        low = bins[i]
        high = bins[i+1]
        m = 5
        tilde_array[i] = 4*np.pi*A*n_sat_glob*(1/b)**(a-3)*romberg(integrand,low,high,m)
        
    like_value = -sum(N_array*np.log(tilde_array)-tilde_array)
    
    return like_value

def opti(lam): 
    #in our case it should be like
    return like(x_glob+lam*n_glob[i][0],y_glob+lam*n_glob[i][1],z_glob+lam*n_glob[i][2])

for j in range(5):
    x_max_glob = max_x_array[j]
    n_sat_glob = n_sat_array[j]
    halo_glob = halo_array[j]
    r_glob = r_array[j]
    #calculating N
    N_array = N_array_array[j]
    
    
    bins = np.linspace(0,x_max_glob,75)
    x_glob, y_glob, z_glob = 1.5,0.5,1 #random initial guess
    initial = np.array([x_glob, y_glob, z_glob])
    
    g = np.array([],dtype=object)
    n_glob = np.array([],dtype=object)
    gamma = np.array([])
    lam_track = np.array([])
    i = 0
    while True:
        x_glob, y_glob, z_glob = initial[0],initial[1],initial[2]
        g = np.append(g, 0)
        #calculate the direction of steepes decent for an initial point
        g_new = - gradient(like,initial)
        g[i] = g_new
        
        # if i = 0 set gamma_i = 0
        if i == 0:
            gamma = np.append(gamma,0)
            n_glob = np.append(n_glob, 0)

            #set new direction
            n_glob[0] = g[0] #n[i-1] nicht da und gamma_i = 0



        else: 
            #else calculate gamma_i
            gamma_val = scalar((g[i]-g[i-1]),g[i])/ab(g[i-1])**2 #gamma_i is a scalar, g_i is a vector 
            gamma = np.append(gamma,gamma_val)

            #set new direction
            n_val = g[i]+gamma[i]*n_glob[i-1]  
            n_glob = np.append(n_glob, 0)
            n_glob[i] = n_val


        #now do line minimization to find lambda for new x point via golden ratio search
        mini = 0 #because of graddient I would guess itstarts at zero!
        maxi = 10**-1
        
    
        mid = mini +  (maxi-mini)/2

        #provide me with the minimum
        lam_min = golden(opti,mini,mid,maxi,10**-4)
        
        lam_track = np.append(lam_track, lam_min)
        #old parameter
        xix, xiy, xiz = initial[0],initial[1],initial[2]

        #calculate new parameter
        initial = initial+lam_min*n_glob[i]

        #new parameter
        xip1x, xip1y, xip1z = initial[0],initial[1],initial[2]

        #check for convergance by comparing to its target accuracy
        check = 2*abs(like(xix, xiy, xiz)-like(xip1x, xip1y, xip1z))/abs(like(xix, xiy, xiz)+like(xip1x, xip1y, xip1z))

        print("check",check)
        print(f"currently a:{xip1x}, b:{xip1y}, c:{xip1z}")

        if check < 10**-3:
            print("reached target accuracy")
            print(f"a:{xip1x}, b:{xip1y}, c:{xip1z}")
            break 

        i += 1
        
    print("safe parameter")
    with open(f'dataset_{j}_paramsc.txt', 'w') as f:
        f.write(f'For dataset {j+1}\n')
        f.write(f'Nsat: {n_sat_glob}\n')
        f.write(f'Best-fit parameters:\n')
        f.write(f"a:{xip1x}, b:{xip1y}, c:{xip1z}")
        f.write(f'Minimum value: {opti(lam_min)}\n')
        
    np.save(f"chi2{j}c.npy",opti(lam_min))
    
    print(f"create plot {j}")
    a,b,c = xip1x, xip1y, xip1z
    tilde_array = np.zeros(74) #number of bins
    A = A_value(a,b,c,5)

    def integrand(x):
        return x**(a-1)*np.exp(-(x/b)**c)

    #print("now calculating tilde")
    for i in range(75-1):
        low = bins[i]
        high = bins[i+1]
        m = 5
        tilde_array[i] = 4*np.pi*A*n_sat_glob*(1/b)**(a-3)*romberg(integrand,low,high,m)
        
    np.save(f"tilde{j}c.npy",tilde_array)
        
    starts = bins[:-1] + (bins[1:]-bins[:-1])/2
    plt.bar(starts, tilde_array, width=0.015, align ="center",label = "best fit data")
    plt.bar(starts, N_array, width=0.015, align ="center", label = "data")
    #bisher kein log log plot, since 
    plt.xlabel("x=r/r_vir")
    plt.ylabel("N(x)")
    if check == 0:
        plt.title(f'Dataset {j+1}, c), error in golden ratio')
    else: 
        plt.title(f'Dataset {j+1}, c)')
    plt.legend()
    plt.savefig(f'dataset_{j+1}_plotc.jpg')
    plt.close()
    """problems with golden ratio search and check point! - because grs provides twice the same value-> ckeck = 0 -> end"""
        
    
