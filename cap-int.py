import numpy as np
from scipy.special import comb
from scipy.special import gamma
from scipy.special import gammainc


def complexint(a,lmn1,A,b,lmn2,B):

    l1,m1,n1 = lmn1    # ANGULAR MOMENTUM SHELL ON ATOM A

    l2,m2,n2 = lmn2    # ANGULAR MOMENTUM SHELL ON ATOM B

#A=[0.0, 0.0, 0.0]

#B=[0.0, 0.0, 1.42]

#l1 = 1
#m1 = 1
#n1 = 1
#l2 = 1
#m2 = 1
#n2 = 1

#a = 3.425250914
#b = 0.6239137298

alpha = a + b

n = 2              # QUADRATIC CAP

 # CAP BOX PARAMETERS-----------------------
c1 = 0.00000

c2 = 0.00000

c3 = 0.00000
 # -----------------------------------------

AB1 = ((a*A[0]) + (b*B[0]))/alpha

AB2 = ((a*A[1]) + (b*B[1]))/alpha

AB3 = ((a*A[2]) + (b*B[2]))/alpha


# CONSTRUCTION OF CAP ALONG X-DIRECTION    

cap_xi_X = 0.0

for rho in range(0,l1):
    
    for sigma in range(0,l2):
        
        I1 = 0.0

        k = rho+sigma
                   
        for tau in range(0,n):
            
            incgamma1 = gamma((k+tau+1)*0.5)*gammainc((k+tau+1)*0.5,alpha*(c1+AB1)*(c1+AB1))
            incgamma2 = gamma((k+tau+1)*0.5)*gammainc((k+tau+1)*0.5,alpha*(c1-AB1)*(c1-AB1))
            
            if (c1+AB1) >= 0:
                
                xi = 1
            
            else:
                           
                xi = np.power(-1,k+tau+1)

            delta1 = gamma((k+tau+1)*0.5)-(incgamma1*xi)

            if (c1-AB1) >= 0:
                           
                xi = 1
                       
            else:
                           
                xi = np.power(-1,k+tau+1)

            delta2 = gamma((k+tau+1)*0.5)-(incgamma2*xi)

            I1 += np.power(-1,n-tau)*comb(n,tau)*np.power(alpha,-0.5*(k+tau+1))*(np.power(-1,k)*np.power((c1+AB1),n-tau)*delta1 + np.power((c1-AB1),n-tau)*delta2)

        I1 = 0.5*I1

        cap_xi_X += comb(l1,rho)*comb(l2,sigma)*np.power((AB1-A[0]),l1-rho)*np.power((AB1-B[0]),l2-sigma)*I1

# CONSTRUCTION OF CAP ALONG Y-DIRECTION    

cap_xi_Y = 0.0

for rho in range(0,m1):
    
    for sigma in range(0,m2):
        
        I1 = 0.0

        k = rho+sigma
                   
        for tau in range(0,n):
            
            incgamma1 = gamma((k+tau+1)*0.5)*gammainc((k+tau+1)*0.5,alpha*(c2+AB2)*(c2+AB2))
            incgamma2 = gamma((k+tau+1)*0.5)*gammainc((k+tau+1)*0.5,alpha*(c2-AB2)*(c2-AB2))
            
            if (c2+AB2) >= 0:
                
                xi = 1
            
            else:
                           
                xi = np.power(-1,k+tau+1)

            delta1 = gamma((k+tau+1)*0.5)-(incgamma1*xi)

            if (c2-AB2) >= 0:
                           
                xi = 1
                       
            else:
                           
                xi = np.power(-1,k+tau+1)

            delta2 = gamma((k+tau+1)*0.5)-(incgamma2*xi)

            I1 += np.power(-1,n-tau)*comb(n,tau)*np.power(alpha,-0.5*(k+tau+1))*(np.power(-1,k)*np.power((c2+AB2),n-tau)*delta1 + np.power((c2-AB2),n-tau)*delta2)

        I1 = 0.5*I1

        cap_xi_Y += comb(m1,rho)*comb(m2,sigma)*np.power((AB2-A[1]),m1-rho)*np.power((AB2-B[1]),m2-sigma)*I1

# CONSTRUCTION OF CAP ALONG Z-DIRECTION    

cap_xi_Z = 0.0

for rho in range(0,n1):
    
    for sigma in range(0,n2):
        
        I1 = 0.0

        k = rho+sigma
                   
        for tau in range(0,n):
            
            incgamma1 = gamma((k+tau+1)*0.5)*gammainc((k+tau+1)*0.5,alpha*(c3+AB3)*(c3+AB3))
            incgamma2 = gamma((k+tau+1)*0.5)*gammainc((k+tau+1)*0.5,alpha*(c3-AB3)*(c3-AB3))
            
            if (c3+AB3) >= 0:
                
                xi = 1
            
            else:
                           
                xi = np.power(-1,k+tau+1)

            delta1 = gamma((k+tau+1)*0.5)-(incgamma1*xi)

            if (c3-AB3) >= 0:
                           
                xi = 1
                       
            else:
                           
                xi = np.power(-1,k+tau+1)

            delta2 = gamma((k+tau+1)*0.5)-(incgamma2*xi)

            I1 += np.power(-1,n-tau)*comb(n,tau)*np.power(alpha,-0.5*(k+tau+1))*(np.power(-1,k)*np.power((c3+AB3),n-tau)*delta1 + np.power((c3-AB3),n-tau)*delta2)

        I1 = 0.5*I1

        cap_xi_Z += comb(n1,rho)*comb(n2,sigma)*np.power((AB3-A[2]),n1-rho)*np.power((AB3-B[2]),n2-sigma)*I1
