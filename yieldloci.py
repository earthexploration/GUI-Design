# coding -*- uft-8 -*-
#
import sys
import scipy.optimize as opt
import numpy as np
import math
import pyqtgraph as pg
from PyQt4 import QtCore,QtGui

# Plane Stress
class YL:
    def _init__(self):
        self.yieldcrits = {'Vonmises':self.Vonmises,
                           'Hill48':self.Hill48,
                           'Yld2000':self.Yld2000_plot,
                           'CPB':self.CPB_plot
                           }
        self.Vonmises_params = None
        self.Hill48_params = None
        self.Yld2000_params = None
        self.CPB_params = None
##        yieldparams = {'Vonmises':self.Vonmises_params,
##                       'Hill48':self.Hill48_params,
##                       'Yld2000':self.Yld2000_params,
##                       'CPB':self.CPB_params,
##                      }
        
    def params(self,name):
        if name == "Vonmises":
            return self.Vonmises_params
        if name == "Hill48":
            return self.Hill48_params
        if name == "Yld2000":
            return self.Yld2000_params
        if name == "CPB":
            return self.CPB_params

    def Vonmises(self,ay):
        # a1^2-a1*a2+a2^2=ay^2
        # (R,theta) <-> (a1,a2)
        # R^2*[1-cos(theta)*sin(theta)]=ay^2
        point = {'X':[],'Y':[]}
        for theta in np.linspace(0,2*np.pi):
            tmp = 1 - np.cos(theta)*np.sin(theta)
            R = np.absolute(ay)/np.sqrt(tmp)
            x = R*np.cos(theta)
            y = R*np.sin(theta)
            point['X'].append(x)
            point['Y'].append(y)
        self.Vonmises_params = [ay]
        print "Parameters:",ay
        return point

    def Hill48(self,a0,a90,r0,r90):
        # (G+H)*axx^2-2H*axx*ayy+(H+F)*ayy^2+2N*axy^2=1
        # r0 = H/G; r90 = H/F; r45 = H/(F+G)-1/2
        # a0/a90 = sqrt(r0*(1+r90)/r90/(1+r0))
        # a1^2-2*r0/(1+r0)*a1*a2+r0*(1+r90)/r90/(1+r0)*a2^2=
        # a0^2;r0*(1+r90)/r90/(1+r0)*a90^2
        point = {'X':[],'Y':[]}
        if a90 == None and a0 is not None:
            GH=1/a0**2
            H=r0/(1+r0)/a0**2
            HF=r0*(1+r90)/r90/(1+r0)/a0**2
            G=GH-H
            F=HF-H
        if a0 == None and a90 is not None:
            a0_2 = r0*(1+r90)/r90/(1+r0)*a90**2
            GH=1/a0_2
            H=r0/(1+r0)/a0_2
            HF=r0*(1+r90)/r90/(1+r0)/a0_2
            G=GH-H
            F=HF-H
        self.Hill48_params = [G,H,F]
        print "Parameters:",G,H,F
        for theta in np.linspace(0,2*np.pi):
            tmp = GH*np.cos(theta)**2-2*H*np.cos(theta)*np.sin(theta)+HF*np.sin(theta)**2
            R = 1./np.sqrt(tmp)
            x = R*np.cos(theta)
            y = R*np.sin(theta)
            point['X'].append(x)
            point['Y'].append(y)
        return point

    def Yld2000_a(self,Ytheta,theta,a,alpha):
        if theta=='b':
            axx = Ytheta
            ayy = Ytheta
            axy = 0
        else:
            theta = theta/180*math.pi
            axx = Ytheta*np.cos(theta)**2
            ayy = Ytheta*np.sin(theta)**2
            axy = Ytheta*np.cos(theta)*np.sin(theta)
        # a=8 for fcc 6 for bcc
        (alpha1,alpha2,alpha3,alpha4,alpha5,alpha6,alpha7,alpha8)=alpha
        # L',L''
        L1_11 = 2.0/3*alpha1
        L1_12 = -1.0/3*alpha1
        L1_21 = -1.0/3*alpha2
        L1_22 = 2.0/3*alpha2
        L1_66 = alpha7
        L2_11 = 1.0/9*(-2*alpha3+2*alpha4+8*alpha5+-2*alpha6)
        L2_12 = 1.0/9*(1*alpha3+-4*alpha4+-4*alpha5+4*alpha6)
        L2_21 = 1.0/9*(4*alpha3+-4*alpha4+-4*alpha5+1*alpha6)
        L2_22 = 1.0/9*(-2*alpha3+8*alpha4+2*alpha5+-2*alpha6)
        L2_66 = alpha8
        # X',X''
        X1_xx = L1_11*axx+L1_12*ayy
        X1_yy = L1_21*axx+L1_22*ayy
        X1_xy = L1_66*axy
        X2_xx = L2_11*axx+L2_12*ayy
        X2_yy = L2_21*axx+L2_22*ayy
        X2_xy = L2_66*axy
        # X1',X2'
        X1_11 = 1.0/2*(X1_xx+X1_yy+np.sqrt((X1_xx-X1_yy)**2+4*X1_xy**2))
        X1_22 = 1.0/2*(X1_xx+X1_yy-np.sqrt((X1_xx-X1_yy)**2+4*X1_xy**2))
        # X1'',X2''
        X2_11 = 1.0/2*(X2_xx+X2_yy+np.sqrt((X2_xx-X2_yy)**2+4*X2_xy**2))
        X2_22 = 1.0/2*(X2_xx+X2_yy-np.sqrt((X2_xx-X2_yy)**2+4*X2_xy**2))
        # phi',phi'',phi
        Phi1 = abs(X1_11-X1_22)**a
        Phi2 = abs(2*X2_22+X2_11)**a+abs(2*X2_11+X2_22)**a
        Phi = Phi1 + Phi2
##        print 'a_eff:',(Phi/2.)**(1./a)
        return (Phi/2.0)**(1.0/a)

    def Yld2000_r(self,Ytheta,theta,a,alpha):
        if theta=='b':
            axx = Ytheta
            ayy = Ytheta
            axy = 0
        else:
            theta = theta/180*math.pi
            axx = Ytheta*np.cos(theta)**2
            ayy = Ytheta*np.sin(theta)**2
            axy = Ytheta*np.cos(theta)*np.sin(theta)
##            print 'axx,ayy,axy,theta:',axx,ayy,axy,theta/math.pi*180
        (alpha1,alpha2,alpha3,alpha4,alpha5,alpha6,alpha7,alpha8)=alpha
        # L',L''
        L1_11 = 2.0/3*alpha1
        L1_12 = -1.0/3*alpha1
        L1_21 = -1.0/3*alpha2
        L1_22 = 2.0/3*alpha2
        L1_66 = alpha7
        L2_11 = 1.0/9*(-2*alpha3+2*alpha4+8*alpha5+-2*alpha6)
        L2_12 = 1.0/9*(1*alpha3+-4*alpha4+-4*alpha5+4*alpha6)
        L2_21 = 1.0/9*(4*alpha3+-4*alpha4+-4*alpha5+1*alpha6)
        L2_22 = 1.0/9*(-2*alpha3+8*alpha4+2*alpha5+-2*alpha6)
        L2_66 = alpha8
        # X',X''
        X1_xx = L1_11*axx+L1_12*ayy
        X1_yy = L1_21*axx+L1_22*ayy
        X1_xy = L1_66*axy
        X2_xx = L2_11*axx+L2_12*ayy
        X2_yy = L2_21*axx+L2_22*ayy
        X2_xy = L2_66*axy
        # X1',X2'
        X1_11 = 1.0/2*(X1_xx+X1_yy+np.sqrt((X1_xx-X1_yy)**2+4*X1_xy**2))
        X1_22 = 1.0/2*(X1_xx+X1_yy-np.sqrt((X1_xx-X1_yy)**2+4*X1_xy**2))
        # X1'',X2''
        X2_11 = 1.0/2*(X2_xx+X2_yy+np.sqrt((X2_xx-X2_yy)**2+4*X2_xy**2))
        X2_22 = 1.0/2*(X2_xx+X2_yy-np.sqrt((X2_xx-X2_yy)**2+4*X2_xy**2))
        # delta'
        delta1 = (X1_xx-X1_yy)**2+4*X1_xy**2
##        print 'delta1:',delta1
        dPhi1_dX1_11 = a*(X1_11-X1_22)**(a-1)
        dPhi1_dX1_22 = -a*(X1_11-X1_22)**(a-1)
        if abs(delta1)>1e-10:
            dX1_11_dX1_xx = 1.0/2*(1+(X1_xx-X1_yy)/np.sqrt(delta1))
            dX1_11_dX1_yy = 1.0/2*(1-(X1_xx-X1_yy)/np.sqrt(delta1))
            dX1_11_dX1_xy = 2.0*X1_xy/np.sqrt(delta1)
            dX1_22_dX1_xx = 1.0/2*(1-(X1_xx-X1_yy)/np.sqrt(delta1))
            dX1_22_dX1_yy = 1.0/2*(1+(X1_xx-X1_yy)/np.sqrt(delta1))
            dX1_22_dX1_xy = -2.0*X1_xy/np.sqrt(delta1)
            
            dPhi1_dX1_xx = dPhi1_dX1_11*dX1_11_dX1_xx + dPhi1_dX1_22*dX1_22_dX1_xx
            dPhi1_dX1_yy = dPhi1_dX1_11*dX1_11_dX1_yy + dPhi1_dX1_22*dX1_22_dX1_yy
            dPhi1_dX1_xy = dPhi1_dX1_11*dX1_11_dX1_xy + dPhi1_dX1_22*dX1_22_dX1_xy
            
        else:
            dPhi1_dX1_xx = 0
            dPhi1_dX1_11 = 0
            dPhi1_dX1_yy = 0
            dPhi1_dX1_22 = 0
            dPhi1_dX1_xy = 0
##        print "dPhi1_dX1_xx:",dPhi1_dX1_xx
##        print "dPhi1_dX1_yy:",dPhi1_dX1_yy
##        print "dPhi1_dX1_xy:",dPhi1_dX1_xy
        # delta''
        delta2 = (X2_xx-X2_yy)**2+4*X2_xy**2
##        print 'delta2:',delta2
        dPhi2_dX2_11 = a*abs(2*X2_22+X2_11)**(a-1)*np.sign(2*X2_22+X2_11)+\
                               2*a*abs(2*X2_11+X2_22)**(a-1)*np.sign(2*X2_11+X2_22)
        dPhi2_dX2_22 = 2*a*abs(2*X2_22+X2_11)**(a-1)*np.sign(2*X2_22+X2_11)+\
                           a*abs(2*X2_11+X2_22)**(a-1)*np.sign(2*X2_11+X2_22)
        if abs(delta2)>1e-10:
            dX2_11_dX2_xx = 1.0/2*(1+(X2_xx-X2_yy)/np.sqrt(delta2))
            dX2_11_dX2_yy = 1.0/2*(1-(X2_xx-X2_yy)/np.sqrt(delta2))
            dX2_11_dX2_xy = 2.0*X2_xy/np.sqrt(delta2)
            dX2_22_dX2_xx = 1.0/2*(1-(X2_xx-X2_yy)/np.sqrt(delta2))
            dX2_22_dX2_yy = 1.0/2*(1+(X2_xx-X2_yy)/np.sqrt(delta2))
            dX2_22_dX2_xy = -2.0*X2_xy/np.sqrt(delta2)
            
            dPhi2_dX2_xx = dPhi2_dX2_11*dX2_11_dX2_xx + dPhi2_dX2_22*dX2_22_dX2_xx
            dPhi2_dX2_yy = dPhi2_dX2_11*dX2_11_dX2_yy + dPhi2_dX2_22*dX2_22_dX2_yy
            dPhi2_dX2_xy = dPhi2_dX2_11*dX2_11_dX2_xy + dPhi2_dX2_22*dX2_22_dX2_xy
            
        else:
            dPhi2_dX2_xx = dPhi2_dX2_11
            dPhi2_dX2_yy = dPhi2_dX2_22
            dPhi2_dX2_xy = 0
##        print "dPhi2_dX2_xx:",dPhi2_dX2_xx
##        print "dPhi2_dX2_yy:",dPhi2_dX2_yy
##        print "dPhi2_dX2_xy:",dPhi2_dX2_xy
        dPhi_daxx = dPhi1_dX1_xx*L1_11+dPhi1_dX1_yy*L1_21+dPhi2_dX2_xx*L2_11+dPhi2_dX2_yy*L2_21
        dPhi_dayy = dPhi1_dX1_xx*L1_12+dPhi1_dX1_yy*L1_22+dPhi2_dX2_xx*L2_12+dPhi2_dX2_yy*L2_22
        dPhi_daxy = dPhi1_dX1_xy*L1_66+dPhi2_dX2_xy*L2_66
##        print 'dPhi_daxx:',dPhi_daxx
##        print 'dPhi_dayy:',dPhi_dayy
##        print 'dPhi_daxy:',dPhi_daxy
        if theta=='b':
            r=dPhi_dayy/dPhi_daxx
        else:
            r=-(np.sin(theta)**2*dPhi_daxx-np.sin(2*theta)*dPhi_daxy+np.cos(theta)**2*dPhi_dayy)/\
               (dPhi_daxx+dPhi_dayy)
        return r

    def Yld2000_err(self,alpha,a0,a45,a90,ab,r0,r45,r90,rb,a):
        [alpha1,alpha2,alpha3,alpha4,alpha5,alpha6,alpha7,alpha8]=alpha
        Yref = self.Yld2000_a(a0,0,a,alpha)
        err_a0  = (Yref - a0)/Yref
        err_a45 = (self.Yld2000_a(a45,45,a,alpha) - Yref)/Yref
        err_a90 = (self.Yld2000_a(a90,90,a,alpha) - Yref)/Yref
        err_ab = (self.Yld2000_a(ab,'b',a,alpha) - Yref)/Yref
        
        err_r0 = (self.Yld2000_r(a0,0,a,alpha) - r0)/r0
        err_r45 = (self.Yld2000_r(a45,45,a,alpha) - r45)/r45
        err_r90 = (self.Yld2000_r(a90,90,a,alpha) - r90)/r90
        err_rb = (self.Yld2000_r(ab,'b',a,alpha) - rb)/rb

##        error = (err_ab/ab)**2+(err_a0/a0)**2+(err_a45/a45)**2+(err_a90/a90)**2+\
##                (err_rb/rb)**2+(err_r0/r0)**2+(err_r45/r45)**2+(err_r90/r90)**2
        #print "error:",error
        return [err_a0,err_a45,err_a90,err_ab,err_r0,err_r45,err_r90,err_rb]

    def Yld2000_alpha(self,a0,a45,a90,ab,r0,r45,r90,rb,a):
        Init=[1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0]
        Alpha = opt.leastsq(self.Yld2000_err,Init,args=(a0,a45,a90,ab,r0,r45,r90,rb,a))
        print "Parameters:", Alpha
        return Alpha[0]
    def Yld2000_point(self,axx,ayy,axy,a,alpha):
        # a=8 for fcc 6 for bcc
        [alpha1,alpha2,alpha3,alpha4,alpha5,alpha6,alpha7,alpha8]=alpha
        # L',L''
        L1_11 = 2.0/3*alpha1
        L1_12 = -1.0/3*alpha1
        L1_21 = -1.0/3*alpha2
        L1_22 = 2.0/3*alpha2
        L1_66 = alpha7
        L2_11 = 1.0/9*(-2*alpha3+2*alpha4+8*alpha5+-2*alpha6)
        L2_12 = 1.0/9*(1*alpha3+-4*alpha4+-4*alpha5+4*alpha6)
        L2_21 = 1.0/9*(4*alpha3+-4*alpha4+-4*alpha5+1*alpha6)
        L2_22 = 1.0/9*(-2*alpha3+8*alpha4+2*alpha5+-2*alpha6)
        L2_66 = alpha8
        # X',X''
        X1_xx = L1_11*axx+L1_12*ayy
        X1_yy = L1_21*axx+L1_22*ayy
        X1_xy = L1_66*axy
        X2_xx = L2_11*axx+L2_12*ayy
        X2_yy = L2_21*axx+L2_22*ayy
        X2_xy = L2_66*axy
        # X1',X2'
        X1_11 = 1.0/2*(X1_xx+X1_yy+np.sqrt((X1_xx-X1_yy)**2+4*X1_xy**2))
        X1_22 = 1.0/2*(X1_xx+X1_yy-np.sqrt((X1_xx-X1_yy)**2+4*X1_xy**2))
        # X1'',X2''
        X2_11 = 1.0/2*(X2_xx+X2_yy+np.sqrt((X2_xx-X2_yy)**2+4*X2_xy**2))
        X2_22 = 1.0/2*(X2_xx+X2_yy-np.sqrt((X2_xx-X2_yy)**2+4*X2_xy**2))
        # phi',phi'',phi
        Phi1 = abs(X1_11-X1_22)**a
        Phi2 = abs(2*X2_22+X2_11)**a+abs(2*X2_11+X2_22)**a
        Phi = Phi1 + Phi2
##        print 'a_eff:',(Phi/2.)**(1./a)
        return (Phi/2.0)**(1.0/a)
    
    def Yld2000_plot(self,a0,a45,a90,ab,r0,r45,r90,rb,a):
        alpha = self.Yld2000_alpha(a0,a45,a90,ab,r0,r45,r90,rb,a)
        self.Yld2000_params = alpha
        point={'X':[],'Y':[]}
        Yref = self.Yld2000_point(a0,0,0,a,alpha)
        for theta in np.linspace(0,2*np.pi,300):
            g = lambda Y: self.Yld2000_point(Y*np.cos(theta),Y*np.sin(theta),0,a,alpha)-Yref
            Init = 1.0
            arc = opt.fsolve(g,Init)
            R = arc[0]
            point['X'].append(R*np.cos(theta))
            point['Y'].append(R*np.sin(theta))
        return point

    def CPB_F(self,axx,ayy,azz,ayz,axz,axy,a,Ck):
        # (|Sigma_1|-k*Sigma_1)**a+(|Sigma_2|-k*Sigma_2)**a+(|Sigma_3|-k*Sigma_3)**a=F
        # k in [-1,1];a >= 1
        [C11,C22,C33,C13,C12,C23,k] = Ck
        C66 = 1.0
        C44 = 1.0
        C55 = 1.0

        Sigma_xx = (2./3*C11 - 1./3*C12 - 1./3*C13)*axx + (-1./3*C11 + 2./3*C12 - 1./3*C13)*ayy + (-1./3*C11 - 1./3*C12 + 2./3*C13)*azz
        Sigma_yy = (2./3*C12 - 1./3*C22 - 1./3*C23)*axx + (-1./3*C12 + 2./3*C22 - 1./3*C23)*ayy + (-1./3*C12 - 1./3*C22 + 2./3*C23)*azz
        Sigma_zz = (2./3*C13 - 1./3*C23 - 1./3*C33)*axx + (-1./3*C13 + 2./3*C23 - 1./3*C33)*ayy + (-1./3*C13 - 1./3*C23 + 2./3*C33)*azz
        Sigma_xy = C66*axy
        Sigma_yz = C44*ayz
        Sigma_xz = C55*axz
        if azz==0 and ayz==0 and axz==0:
            Sigma_1 = 1.0/2*(Sigma_xx + Sigma_yy + np.sqrt((Sigma_xx - Sigma_yy)**2 + 4*Sigma_xy**2))
            Sigma_2 = 1.0/2*(Sigma_xx + Sigma_yy - np.sqrt((Sigma_xx - Sigma_yy)**2 + 4*Sigma_xy**2))
            Sigma_3 = Sigma_zz
            
        else:
            array = np.asarray([[Sigma_xx,Sigma_xy,Sigma_xz],[Sigma_xy,Sigma_yy,Sigma_yz],[Sigma_xz,Sigma_yz,Sigma_zz]])
            w,v = np.linalg.eig(array)
            Sigma = [ i for i in np.array(w).tolist()]
            print Sigma
            Sigma_1 = Sigma[0]
            Sigma_2 = Sigma[1]
            Sigma_3 = Sigma[2]
 
        F = (abs(Sigma_1)-k*Sigma_1)**a+(abs(Sigma_2)-k*Sigma_2)**a+(abs(Sigma_3)-k*Sigma_3)**a
        #print 'F:',F      
        return F
    def CPB_a0(self,F,a,Ck):
        [C11,C22,C33,C13,C12,C23,k] = Ck
  
        Phi_1 = (2./3*C11 - 1./3*C12 - 1./3*C13)
        Phi_2 = (2./3*C12 - 1./3*C22 - 1./3*C23)
        Phi_3 = (2./3*C13 - 1./3*C23 - 1./3*C33)
        a0_T = (F/((np.absolute(Phi_1) - k*Phi_1)**a + (np.absolute(Phi_2) - k*Phi_2)**a + (np.absolute(Phi_3) - k*Phi_3)**a))**(1./a)
        a0_C = (F/((np.absolute(Phi_1) + k*Phi_1)**a + (np.absolute(Phi_2) + k*Phi_2)**a + (np.absolute(Phi_3) + k*Phi_3)**a))**(1./a)

        return [a0_T,a0_C]
    def CPB_a90(self,F,a,Ck):
        [C11,C22,C33,C13,C12,C23,k] = Ck
        
        Psi_1 = -1./3*C11 + 2./3*C12 - 1./3*C13
        Psi_2 = -1./3*C12 + 2./3*C22 - 1./3*C23
        Psi_3 = -1./3*C13 + 2./3*C23 - 1./3*C33
        a90_T = (F/((np.absolute(Psi_1) - k*Psi_1)**a + (np.absolute(Psi_2) - k*Psi_2)**a + (np.absolute(Psi_3) - k*Psi_3)**a))**(1./a)
        a90_C = (F/((np.absolute(Psi_1) + k*Psi_1)**a + (np.absolute(Psi_2) + k*Psi_2)**a + (np.absolute(Psi_3) + k*Psi_3)**a))**(1./a)
        return [a90_T,a90_C]
    
    def CPB_aN(self,F,a,Ck):
        [C11,C22,C33,C13,C12,C23,k] = Ck

        Zeta_1 = -1./3*C11 - 1./3*C12 + 2./3*C13
        Zeta_2 = -1./3*C12 - 1./3*C22 + 2./3*C23
        Zeta_3 = -1./3*C13 - 1./3*C23 + 2./3*C33
        aN_T = (F/((np.absolute(Zeta_1) - k*Zeta_1)**a + (np.absolute(Zeta_2) - k*Zeta_2)**a + (np.absolute(Zeta_3) - k*Zeta_3)**a))**(1./a)
        aN_C = (F/((np.absolute(Zeta_1) + k*Zeta_1)**a + (np.absolute(Zeta_2) + k*Zeta_2)**a + (np.absolute(Zeta_3) + k*Zeta_3)**a))**(1./a)

        return [aN_T,aN_C]
    
    def CPB_ab(self,F,a,Ck):
        [C11,C22,C33,C13,C12,C23,k] = Ck
        
        Omega_1 = 1./3*C11 + 1./3*C12 - 2./3*C13
        Omega_2 = 1./3*C12 + 1./3*C22 - 2./3*C23
        Omega_3 = 1./3*C13 + 1./3*C23 - 2./3*C33
        ab_T = (F/((np.absolute(Omega_1) - k*Omega_1)**a + (np.absolute(Omega_2) - k*Omega_2)**a + (np.absolute(Omega_3) - k*Omega_3)**a))**(1./a)
        ab_C = (F/((np.absolute(Omega_1) + k*Omega_1)**a + (np.absolute(Omega_2) + k*Omega_2)**a + (np.absolute(Omega_3) + k*Omega_3)**a))**(1./a)

        return [ab_T,ab_C]

    def CPB_tau0(self,F,a,Ck):
        [C11,C22,C33,C13,C12,C23,k] = Ck
        C66 = 1.0
        tau0 = (F/((np.absolute(C66) + k*C66)**a + (np.absolute(C66) - k*C66)**a))**(1.0/a)
        return tau0
    
    def CPB_err(self,Ck,a0_T,a0_C,a90_T,a90_C,aN_C,ab_T,ab_C,a):
        F = self.CPB_F(a0_T,0,0,0,0,0,a,Ck)
        [a0T,a0C] = self.CPB_a0(F,a,Ck)
        err_a0_T = (a0T - a0_T)/a0_T

        err_a0_C = (a0C - a0_C)/a0_C
        
        [a90T,a90C] = self.CPB_a90(F,a,Ck)
        err_a90_T = (a90T - a90_T)/a90_T
        err_a90_C = (a90C - a90_C)/a90_C
        
        aNC = self.CPB_aN(F,a,Ck)[1]
        err_aN_C = (aNC - aN_C)/aN_C
        
        [abT,abC] = self.CPB_ab(F,a,Ck)
        err_ab_T = (abT - ab_T)/ab_T
        err_ab_C = (abC - ab_C)/ab_C
        
        error = [err_a0_T,err_a0_C,err_a90_T,err_a90_C,err_aN_C,err_ab_T,err_ab_C]
        #print "error:",np.sqrt(sum(np.square(error))/len(error))
        #print error
        return error

    def CPB_Ck(self,a0_T,a0_C,a90_T,a90_C,aN_C,ab_T,ab_C,a):
        Init = [1.0,1.0,1.0,0.0,0.0,0.0,0.0]
        Ck = opt.leastsq(self.CPB_err,Init,args=(a0_T,a0_C,a90_T,a90_C,aN_C,ab_T,ab_C,a))
        print "Parameters:",Ck[0]
        return Ck[0]

    def CPB_plot(self,a0_T,a0_C,a90_T,a90_C,aN_C,ab_T,ab_C,a):
        Ck = self.CPB_Ck(a0_T,a0_C,a90_T,a90_C,aN_C,ab_T,ab_C,a)
        self.CPB_params = Ck
        point={'X':[],'Y':[]}
        #Ck = [1.0,1.0,1.0,0.0,0.0,0.0,-0.2]
        F = self.CPB_F(a0_T,0,0,0,0,0,a,Ck)
        for theta in np.linspace(0,2*np.pi,600):
            g = lambda Y: self.CPB_F(Y*np.cos(theta),Y*np.sin(theta),0.,0.,0.,0.,a,Ck) - F
            Init = 1.0
            arc = opt.fsolve(g,Init)
            R = arc[0]
            point['X'].append(R*np.cos(theta))
            point['Y'].append(R*np.sin(theta))
        return point

##def main():
##    app = QtGui.QApplication(sys.argv)
##    
##    yl = YL()
    #point = yl.CPB_plot(1.0,1.0,1.0,1.0,1.0,1.0,1.0,2.0)

##    point = yl.Yld2000_plot(1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,8.0)
##    pw = pg.plot()
##    pw.showGrid(True,True,alpha=1.0)
##    
##    pw.plot(point['X'],point['Y'])
##    
##    
##    app.exec_()
    
##if __name__ == "__main__":
##    main()

        
