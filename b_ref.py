import numpy as np
import matplotlib.pyplot as plt

#import datasets
def readfile(filename):
    f = open(filename, 'r')
    data = f.readlines()[3:] #Skip first 3 lines 
    nhalo = int(data[0]) #number of halos
    radius = []
    
    for line in data[1:]:
        if line[:-1]!='#':
            radius.append(float(line.split()[0]))
    
    radius = np.array(radius, dtype=float)    
    f.close()
    return radius, nhalo



# for integration
#trapezoid rule
def trap(func,low,high,N): #a=mini, b = maxi
    h = (high-low)/N
    x = np.linspace(low,high,N+1)
    result = (func(low)+func(high))*0.5
    for element in x[1:N-1]: 
        result += func(element)
    result = result*h
    return result 

#romberg integration - from the assignment before
def romberg(func, low, high, m): #a=mini, b = maxi
    h = high - low
    r = np.zeros(m)
    r[0] = trap(func,low,high,1)
    
    #step 4
    Np = 1
    for i in range(1,m): 
        #step 5
        delta = h 
        h = 0.5*h
        x = low + h #a = mini
        for _ in range(Np): 
            r[i] = r[i] + func(x)
            x = x + delta

        #step 6
        r[i] = 0.5*(r[i-1]+ delta*r[i])
        #then double Np
        Np *=2
        
    #step 8
    Np = 1
    for i in range(1,m):
        #step 9
        Np *= 4
        for j in range(m-i):
            r[j] = (Np*r[j+1]-r[j])/(Np-1)
            
    return r[0]


#for linear minimalization
def golden(f,mini,mid,maxi,acc):
    # for test
    x_array = np.array([[None],[None],[None]])
    count = 0
    while True: 
        
        # step 1
        interval_1 = abs(maxi - mid)
        interval_2 = abs(mid - mini)

        if interval_1 > interval_2:
            larger_interval = (mid, maxi)
            x = maxi

        else:
            larger_interval = (mini, mid)
            x = mini
        
        #also for later! - to check whether this loop is stuck
        x_array[count%3] = x
        if (x_array[0]==x_array[1]) and (x_array[0]==x_array[2]):
            value = x
            break
        
        
        w = 0.38197 #golden ratio
        d = mid + (x-mid)*w
        
        # step 2
        if abs(maxi - mini) < acc:
            if f(d) < f(mid):
                #print(f"d it is!. d = {d}")
                #print("a",mini)
                #print("b",mid)
                #print("c",maxi)
                value = d
                break
            else: 
                #print(f"b it is!. b = {mid}")
                #print("a",mini)
                #print("c",maxi)
                #print("d",d)
                value = mid
                break

        # step 3
        #tighten towards d
        if f(d) < f(mid): 
            if mid < d < maxi: 
                mini,mid = mid,d

            else:
                maxi, mid = mid ,d

        #step 4
        #tighten towads b
        if f(d) > f(mid): 
            if mid < d < maxi:
                maxi = d

            else:
                mini = d
        count += 1
    return value

#define absolute value of a vector
def ab(vec):
    return np.sqrt(sum(vec**2))

#scalar product of a vector
def scalar(vec1,vec2):
    return sum(vec1*vec2)

#gradient of a scalar field 
def gradient(function,vec): # vec_x = x,y,z
    x,y,z = vec[0],vec[1],vec[2]
    array = np.zeros(3)
    h = 0.0001 #might causes issues for x,y or z <<h 
    #central method 
    array[0] = (function(x+h,y,z)-function(x-h,y,z))/2/h
    array[1] = (function(x,y+h,z)-function(x,y-h,z))/2/h
    array[2] = (function(x,y,z+h)-function(x,y,z-h))/2/h
    return array

#for linear minimization
def opti(lam): 
    return chi2_eq(x_glob+lam*n_glob[i][0],y_glob+lam*n_glob[i][1],z_glob+lam*n_glob[i][2]) #given calculated values 

#calculate A for a new set of a,b,c
def A_value(a,b,c,x_max):
    if (a < 1) or (c < 0):
        print("Warning") #important for minimization for specific case!
        
    m = 10 #seemed (tried out) to be a good value
    mini = 0
    maxi = x_max
    def integrad_a(x):
        return x**(a-1)*np.exp(-(x/b)**c)
    
    k = 4*np.pi*(1/b)**(a-3)*romberg(integrad_a,mini,maxi,m)
    return 1/k # returns A

#chi^2 function
def chi2_eq(a,b,c):
    #N_array is given
    """x_max_glob"""
    A = A_value(a,b,c,x_max_glob)
    #tilde_i
    """x_max_glob,N_sat_glob"""
    tilde_array = np.zeros(74) #number of bins
    
    def integrand(x):
        return x**(a-1)*np.exp(-(x/b)**c)
    
    #print("now calculating tilde")
    for i in range(75-1):
        low = bins[i]
        high = bins[i+1]
        m = 5
        tilde_array[i] = 4*np.pi*A*n_sat_glob*(1/b)**(a-3)*romberg(integrand,low,high,m)
    
    chi_value = sum((N_array-tilde_array)**2/tilde_array**2)
    
    return chi_value 

for i in range(1,6):
    radius, nhalo = readfile(f'satgals_m1{i}.txt')
    np.save(f"radius{i}.npy", radius)
    np.save(f"nhalo{i}.npy", nhalo)

if __name__ == "__main__":
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

    #the conjugated gradient method for the different datasets
    for j in range(5):
        x_max_glob = max_x_array[j]
        n_sat_glob = n_sat_array[j]
        halo_glob = halo_array[j]
        r_glob = r_array[j]
        
        
        bins = np.linspace(0,x_max_glob,75)
        x_glob, y_glob, z_glob = 1.5,0.5,1 #random initial guess
        initial = np.array([x_glob, y_glob, z_glob])
        
        g = np.array([],dtype=object)
        n_glob = np.array([],dtype=object)
        gamma = np.array([])
        lam_track = np.array([])
        
        #calculating N
        
        N_array = np.zeros(len(bins)-1)
        for i in range(len(bins)-1):
            #number of galaxies in bin
            N_array[i] = sum(r_glob* (bins[i] < r_glob)*(r_glob < bins[i+1]))  #returns True/False = 1/0
        
        #mean number of galaxy per bin per halo
        N_array = N_array/halo_glob 
            
        np.save(f"n_array{j}.npy",N_array)

        #start of the actual method
        maxi = 0.01
        i = 0
        while True:
            x_glob, y_glob, z_glob = initial[0],initial[1],initial[2]
            g = np.append(g, 0)
            #calculate the direction of steepes decent for an initial point
            g_new = - gradient(chi2_eq,initial)
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
            mini = 0 
            #to reduce mistakes by golden ratio search
            if i>2:
                lam_min_nacher = lam_track[i-1]
                lam_min_vorher = lam_track[i-2]

                if lam_min_vorher == lam_min_nacher:
                    maxi = maxi/2
            

            mid = mini +  (maxi-mini)/2

            #provide me with the minimum
            lam_min = golden(opti,mini,mid,maxi,10**-4)
            
            lam_track = np.append(lam_track, lam_min)
            xix, xiy, xiz = initial[0],initial[1],initial[2]

            #calculate new parameter
            initial = initial+lam_min*n_glob[i]

            #new parameter
            xip1x, xip1y, xip1z = initial[0],initial[1],initial[2]

            #check for convergance by comparing to its target accuracy
            check = 2*abs(chi2_eq(xix, xiy, xiz)-chi2_eq(xip1x, xip1y, xip1z))/abs(chi2_eq(xix, xiy, xiz)+chi2_eq(xip1x, xip1y, xip1z))

            print("check",check)
            print(f"currently a:{xip1x}, b:{xip1y}, c:{xip1z}")

            if check < 10**-3:
                print("reached target accuracy")
                print(f"a:{xip1x}, b:{xip1y}, c:{xip1z}")
                break 

            i += 1
            
        print("safe parameter")
        with open(f'dataset_{j}_paramsb.txt', 'w') as f:
            f.write(f'For dataset {j+1}\n')
            f.write(f'Nsat: {n_sat_glob}\n')
            f.write(f'Best-fit parameters:\n')
            f.write(f"a:{xip1x}, b:{xip1y}, c:{xip1z}")
            f.write(f'Minimum value chi2: {opti(lam_min)}\n')
        np.save(f"chi2{j}b.npy",opti(lam_min))
            
        print("create plot")
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
        
        np.save(f"tilde{j}b.npy",tilde_array)
            
        starts = bins[:-1] + (bins[1:]-bins[:-1])/2
        plt.bar(starts, tilde_array, width=0.015, align ="center",label = "best fit data")
        plt.bar(starts, N_array, width=0.015, align ="center", label = "data")
        #bisher kein log log plot, since 
        plt.xlabel("x=r/r_vir")
        plt.ylabel("N(x)")
        if check == 0:
            plt.title(f'Dataset {j+1}, b), error in golden ratio')
        else: 
            plt.title(f'Dataset {j+1}, b)')
        plt.legend()
        plt.savefig(f'dataset_{j+1}_plotb.jpg')
        plt.close()
        """problems with golden ratio search and check point! - because grs provides twice the same value-> ckeck = 0 -> end"""
            
    
