# PROBLEM SET 1 QUESTION 2

import numpy as np 
import matplotlib.pyplot as plt 

plt.rcParams.update({"text.usetex" : True})

integrate = True
if integrate:
    def simpsons(
            func: callable, 
            bounds: list[int|float], 
            num_slices: int|float
            ) -> list[int|float]:
        """
        Compute the single-variabel integral of a function using quadratic curves
        """
        if len(bounds) != 2:
            raise ValueError('Input bounds must be array of length 2!')
        
        if num_slices%2 != 0:
            raise ValueError("Input variable 'num_slices' must be an even value!")

        lower = bounds[0]
        upper = bounds[1]
        h = (upper-lower) / num_slices #width of each slice 
        A = np.zeros(int(num_slices/2)) #blank array for output length
        for k in range(int(num_slices/2)):
            A[k] = (1/3)*h*(func(lower + 2*k*h, T, alpha) + 4*func(lower + h*(1+2*k), T, alpha ) + func(lower + 2*h*(1+k), T, alpha)) #sum of lefthand, center, and righthand
        
        value =  np.sum(A)
        return value






    T = 1/52 #1 week; in years
    alpha = 15 * ((3600*24*365)**2)*1.057e-16 #in light years / years^2



    def dx_dir(tau, T, alpha):
        C = (alpha * T) / (2*np.pi)
        pre_sinh = C*(1-np.cos((2*np.pi*tau)/T))
        return np.sinh(pre_sinh)

    def dt_dir(tau, T, alpha):
        C = (alpha * T) / (2*np.pi)
        pre_cosh = C*(1-np.cos((2*np.pi*tau)/T))
        return np.cosh(pre_cosh)



    N=100
    tau = np.linspace(0, T, 6048)

    #perform numerical integration on dt and dx from 0 to tau

    def integrate_timestep(tau):
        start = 0
        stop = tau
        x = simpsons(dx_dir, (start, stop), N)
        t = simpsons(dt_dir, (start, stop), N)
        return x, t

    x_empty = np.zeros(len(tau))
    t_empty = np.zeros(len(tau))
    for k in range(len(tau)):
        x_empty[k] = integrate_timestep(tau[k])[0]
        t_empty[k] = integrate_timestep(tau[k])[1]

    x = x_empty
    t = t_empty



    plt.plot(x, t, label=r'Spacetime Coordinates ($\tau$ = 1 week)', lw=2, color='b')
    plt.title('Worldline of Twin Paradox Astronaut', fontsize=16)
    plt.xlabel('$x^1$ [light years]', fontsize=14)
    plt.ylabel('$x^0$ [years]', fontsize=14)
    plt.legend(loc='best', fontsize=14)
    plt.grid(alpha = 0.4)
    plt.show()


    print('Time Dilation for two-way trip={} seconds'.format(2*(np.max(t) - 1/52)*365*24*3600))
    # 0.01923091257063526 years
    # of time in lab frame 
    # time in proper frame  = 1/52 years 








homebody = np.linspace(0, 2, 400)
astronaut_to = np.linspace(0, 1, 400)
astronaut_back = np.linspace(1, 0, 400)
time = np.linspace(0, 1, 400)

xs1 = (0, 0, 0, 0,)
ys1 = (0.125, 0.25, 0.375, 0.5)
xs2 = (0, 0, 0, 0)
ys2 = (1.5, 1.625, 1.75, 1.875)

plt.plot(np.zeros(len(homebody)), homebody, color='blue', lw=2, label='Homebody')
plt.plot(time, astronaut_to, color='red', lw=2, label='Astronaut')
plt.plot(time, astronaut_back + 1, color='red', lw=2)
for k in range(len(xs1)):
    plt.plot((xs1[k], 1), (ys1[k],1), ls='--', color='orange', alpha = 0.3)
    plt.plot((xs2[k], 1), (ys2[k],1), ls='--', color='orange', alpha = 0.3)
plt.grid(alpha=0.4)
plt.xlabel('Distance [Arbitrary Units]', fontsize=14)
plt.ylabel('Time [years]', fontsize=14)
plt.title('Worldline of Twins Paradox', fontsize=16)
plt.legend(loc='best', fontsize=14)
plt.show()





# PROBLEM SET 2 QUESTIONS 3 AND 4 

import numpy as np
import sympy as sp
    #https://docs.sympy.org/latest/modules/tensor/array.html

if __name__ == "__main__":
        #create the index list for the indices
    indexlist = (0, 1, 2)
        #define variables     
    L, t, theta, phi= sp.symbols('L t theta phi')
        #define tensor indice: 0, 1, 2
    mu, nu, ka, si, la, rho = sp.symbols('mu nu kappa sigma lambda rho')

        #create metric (1, -1, -1)
    eta = np.array([[1, 0, 0],          #one time dimension, 2 space dimensions
                    [0, -1, 0],
                    [0, 0, -1]], dtype=object)

        #write xupup as a 3-vector in matrix representation
    xup = np.array(([t, theta, phi]))
        #write gdndn as a 3x3-tensor in matrix representation
    gdndn = np.array([[L**2, 0, 0], 
                    [0, -(L**2)*(sp.cosh(t)**2),0], 
                    [0, 0, -(L**2)*(sp.sin(theta)**2)*(sp.cosh(t)**2)]], dtype = object)
        #compute gupup by taking the inverse of gdndn
    gupup= sp.Matrix(gdndn).inv()
        #re-establish gupup as a numpy array to perform operations on it
    gupup = np.array(gupup, dtype = object)




        #contract the tensor product along the 2nd index of gdndn and 1st index of gupup
    # delta = sp.tensorcontraction(sp.tensorproduct(gdndn, gupup), (2, 1))  #this was one way to do it according to the sympy wiki, don't understand how it works though
        #do it another way: 
        #create an array to perform a 1-dimensional matrix contraction (via dot product)
    vdotfill = (1, 1, 1)
        # create empty 3x3 array to append
    delta = np.empty((3,3), dtype=object)
        #create empty 1x3 array to perform contraction over 1 of the 3 indices
    element = np.empty(3, dtype=object)
        #iterate through each index
    for la in indexlist:
        for nu in indexlist:
            for mu in indexlist:    #append the concraction array with elements of gupup / gdndn
                element[mu] = gdndn[la][mu]*gupup[mu][nu]
            delta[la][nu] = sp.simplify(element.dot(vdotfill))   #contract over the mu index by summing via a dot product

    print('kronecker delta')
    print(delta)    #yields delta-function
    print('------') #print a space between answers





        #do something similar
        #create empty 3x3x3 array filled with any kind of object (we will be working with sympy matrices as tensor matrix representations)
    Chrisupdndn = np.empty((3,3,3), dtype=object)
        #create an empty 1x3 array, again for any object, since we will need to perform a contraction over 1 of the 4 indices
    v = np.empty((3), dtype=object)
        #iterate through each index 
    for mu in indexlist:
        for nu in indexlist:
            for la in indexlist:
                for si in indexlist: #compute individual v term mu/nu/lambda/sigma  (0.5 * gupup * d.gdndn + d.gdndn - d.gdndn) along each specified index axis
                    v[si] = 0.5*(gupup[mu][si])*(sp.diff(gdndn[si][la], xup[nu]) + sp.diff(gdndn[nu][si], xup[la]) - sp.diff(gdndn[nu][la], xup[si]))
                Chrisupdndn[mu][nu][la] = sp.simplify(v.dot(vdotfill))  #perform a dot product using vdotfill array to sum the elements of the array (I had no other way to contract the 3-vector)
                                                                        # append the Chrisupdndn array in mu/nu/lambda (contract along sigma for each term) 

    print('christoffel symbol')
    print(Chrisupdndn)  #yields connection coefficients
    print('------') #space



        #create empty 3x3x3x3 matrix representation, taking in any object as dtype
    Riemupdndndn = np.empty((3, 3, 3, 3), dtype=object)
        #create empty 3x1 array to take the dot product with and contract along the lambda axis
    entry = np.empty(3, dtype = object)
        #iterate over every index
    for rho in indexlist:
        for si in indexlist:
            for mu in indexlist:
                for nu in indexlist:
                    for la in indexlist:    #contract ONLY over the lambda index (not the derivatives - only the christoffel contractions)
                        entry[la] = (Chrisupdndn[la][si][nu]*Chrisupdndn[rho][la][mu]) - (Chrisupdndn[la][si][mu]*Chrisupdndn[rho][la][nu])
                    #add the derivative terms to the christoffel contraction terms - DO NOT contract over derivatives which are not contracted
                    #iterate over riemaupdndn and fill with values 
                    Riemupdndndn[rho][si][mu][nu] = sp.simplify(sp.trigsimp(sp.expand(entry.dot(vdotfill)))) + sp.diff(Chrisupdndn[rho][si][nu], xup[mu]) - sp.diff(Chrisupdndn[rho][si][mu], xup[nu])  

    print('riemann tensor')
    print(Riemupdndndn)
    print('------') #space 



        #define Ricci tensor as empty array, then fill it 
    Riccdndn = np.empty((3,3), dtype = object)
        #create empty temporary array to contract lambda indices over, via dot product
    temp = np.empty(3, dtype = object)
        #iterate over indices
    for mu in indexlist:
        for nu in indexlist:
            for la in indexlist:        #contract over la index via dot product of la temp array
                temp[la] = Riemupdndndn[la][mu][nu][la]
            Riccdndn[mu][nu] = sp.simplify(temp.dot(vdotfill))  #contract
        
    print('ricci tensor')
    print(Riccdndn)
    print('------') #space


        #create a loop using the metric to raise a ricci index 
        #empty ricci array prior to contracting mu
    RiccContract = np.empty(3, dtype=object)
        #create temp to contract over nu when raising an index
    temp = np.empty(3, dtype = object)
    for mu in indexlist:
        for nu in indexlist:
            temp[nu] = Riccdndn[mu][nu] * gupup[mu][nu]     #determine the  nu-th term
        RiccContract[mu] = sp.simplify(temp.dot(vdotfill))  #contract over nu 
    RiccScal = sp.simplify(RiccContract.dot(vdotfill))  #now contract over mu 


    print('ricci scalar')
    print(RiccScal)
    print('------') #space




        #specify array for general output of the (reduced) field equations
    output = np.empty((3,3), dtype = object)        
        #specify array to determine the constant \Alpha
    AlphaConstant = np.empty((3,3), dtype = object)
        #specify a temporary array to contract over the nu index
    contr = np.empty(3, dtype = object)
        #iterate over indices
    for mu in indexlist:
        for la in indexlist:
            for nu in indexlist:        #determine the (mu, nu) elements of the output array by programming the field equation in 
                output[mu][nu] = Riccdndn[mu][nu] - 0.5 * gdndn[mu][nu] * RiccScal + (1/L**2)*gdndn[mu][nu]
                                        #contract over the nu index ONlY (do not add other terms like the RiccScal term)
                contr[nu] = - Riccdndn[mu][nu]*gupup[nu][la]
                                        #do dot product to get contraction; then add the 0.5R term on the back to get the proper output
            AlphaConstant[mu][la] = sp.simplify(contr.dot(vdotfill)) + 0.5*RiccScal*delta[mu][la]   

        #print values
    print('field equation test (get zeros) & find Alpha')
    print(output)
    print(AlphaConstant)
    print('------') #space




#for Q1 calculating inverse
# A = np.array([[sp.sin(theta)*sp.cos(phi), L*sp.cos(theta)*sp.cos(phi), -L*sp.sin(theta)*sp.sin(phi)], 
#               [sp.sin(theta)*sp.sin(phi), L*sp.cos(theta)*sp.sin(phi), L*sp.sin(theta)*sp.cos(phi)],
#               [sp.cos(theta), -L*sp.sin(theta), 0]], dtype=object)



# A_inv = sp.Matrix(A).inv()
# A_inv = np.array(A_inv, dtype=object)

# print(sp.simplify(A_inv))
