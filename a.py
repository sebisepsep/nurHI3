import numpy as np
import matplotlib.pyplot as plt
#write maximazation algorithm


def golden(f,mini,mid,maxi,acc):
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

        w = 0.38197 #golden ratio
        d = mid + (x-mid)*w

        # step 2
        if abs(maxi - mini) < acc:
            if f(d) < f(mid):
                print(f"d it is!. d = {d}")
                value = d
                break
            else: 
                print(f"b it is!. b = {mid}")
                value = mid
                break

        # step 3
        #tighten towards d
        if f(d) < f(mid): 
            if mid < d < maxi: #b and c must be ordered :)
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
    return value



def main():
    #given parameters
    a, b, c, x_max, Nsat, A = 2.4, 0.25, 1.6, 5, 100, 256/(5*np.pi**(3/2))
    #given function
    def N(x): 
        return -4*np.pi*A*Nsat*x**(a-1)*(1/b)**(a-3)*np.exp(-(x/b)**c) 
    mini,mid,maxi = 0,2,5
    acc = 10**-6
    xx = golden(N,mini,mid,maxi,acc)
    #saving data
    with open("1a.txt", 'w') as file:
        file.write(f"The minimum is at x={xx},N(x) = {-N(xx)}")
        
if __name__ == "__main__":
    main()
