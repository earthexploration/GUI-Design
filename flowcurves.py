# !/usr/bin/python
# coding: UTF-8
#
import sys
import numpy as np
import scipy.optimize as opt   
import math

from reader import File
import pyqtgraph as pg
from PyQt4 import QtGui,QtCore

class FC:
    def __init__(self,strain,stress,temp=20,rate=0.5,ue=False):
        self.xdata = strain
        self.ydata = stress
        self.models = {'Swift':[self.Swift, self.SwiftObjFunc],
                       'Voce': [self.Voce, self.VoceObjFunc],
                       'Gosh': [self.Gosh, self.GoshObjFunc],
                       'Hockett_Sherby':[self.Hockett_Sherby, self.Hockett_SherbyObjFunc],
                       'SwiftVoce': [self.SwiftVoce, self.SwiftVoceObjFunc],
                       'Modified_Zener_Hollomon': [self.Modified_Zener_Hollomon, self.Modified_Zener_HollomonObjFunc],
                      }
        if ue == False:
            self.w = 0.0
        else:
            self.w = 1e3
            n = stress.index(max(stress))
            self.ue_strain = strain[n]
            self.ue_stress = stress[n]
        self.temp = temp
        self.rate = rate
        #print self.temp,self.rate
    def params(self,name,bounds):
        if name in self.models.keys():
            if name == u'Modified_Zener_Hollomon':
                ans = opt.differential_evolution(self.models[name][1],bounds,args=(self.temp,self.rate,self.xdata,self.ydata,),
                                                 strategy='best1bin',disp=True)
                error = np.sqrt(sum(np.square(self.models[name][0](ans.x,self.temp,self.rate,self.xdata)-self.ydata))/len(self.xdata))
            else:
                ans = opt.differential_evolution(self.models[name][1],bounds,args=(self.xdata,self.ydata,),disp=True)
                error = np.sqrt(sum(np.square(self.models[name][0](ans.x,self.xdata)-self.ydata))/len(self.xdata))
            print 'Model:',name
            print 'Params:',ans.x
            print 'Error:',error
            return [ans.x,error]
        else:
            return None
    def values(self,name,params,xdata):
        if name in self.models.keys():
            if name == u'Modified_Zener_Hollomon':
                values = self.models[name][0](params,self.temp,self.rate,xdata)
                return values
            else:
                values = self.models[name][0](params,xdata)
                return values
            #print 'Model:',name
            #print '(x,y):',xdata,values
        else:
            return None
    
    # models 
    def Swift(self,params,xdata):
        '''stress(strain)=q*(e0+strain)**n'''
        '''[1500,2e-3,0.1]'''
        xdata = np.array(xdata)
        return params[0]*(params[1]+xdata)**(params[2])
    def SwiftPrime(self,params,xdata):
        '''stress'(strain)=q*n*(e0+strain)**(n-1)'''
        xdata = np.array(xdata)
        return params[0]*params[2]*(params[1]+xdata)**(params[2]-1)
    def SwiftJac(self,params,xdata,ydata):
        dfd0 = (params[1]+xdata)**(params[2])
        dfd1 = params[0]*params[2]*(params[1]+xdata)**(params[2]-1)
        dfd2 = params[0]*(params[1]+xdata)**(params[2])*np.log(params[1]+xdata)
        return np.array([dfd0,dfd1,dfd2]).T
    def SwiftObjFunc(self,params,xdata,ydata):
        xdata = np.array(xdata)
        ydata = np.array(ydata)
        addition = self.w*(np.absolute(self.SwiftPrime(params,self.ue_strain)-self.Swift(params,self.ue_strain)))
        return sum(np.absolute(self.Swift(params,xdata)-ydata))+ addition
    def SwiftFit(self,xdata,ydata):
        cons = ({'type':'eq','fun': lambda x: self.SwiftPrime(x,self.ue_strain) - self.ue_stress},
                )
        bnds = ((0,1e4),(0,1e-1),(0,1))
        init = [ydata[0],1e-2,0.1]
        res = opt.minimize(self.SwiftObjFunc,init,args=(xdata,ydata),
                           bounds=bnds,constraints=cons)
        print res
        return res
    
    def Voce(self,params,xdata):
        '''stress(strain) = a + b*(1-exp(-c*strain))'''
        ''' [715,400,38]'''
        xdata = np.array(xdata)
        return params[0]+params[1]*(1-np.exp(-params[2]*xdata))
    def VoceJac(self,params,xdata,ydata):
        xdata = np.array(xdata)
        dfd0 = np.ones(len(xdata))
        dfd1 = (1-np.exp(-params[2]*xdata))
        dfd2 = params[1]*xdata*np.exp(-params[2]*xdata)
        return np.array([dfd0,dfd1,dfd2]).TSwiftVoceObjFunc
    def VocePrime(self,params,xdata):
        '''stress'(strain) = b*c*exp(-c*strain)'''
        xdata = np.array(xdata)
        return params[1]*params[2]*np.exp(-params[2]*xdata)
    def VoceObjFunc(self,params,xdata,ydata):
        xdata = np.array(xdata)
        ydata = np.array(ydata)
        addition = self.w*(np.absolute(self.VocePrime(params,self.ue_strain)-self.Voce(params,self.ue_strain)))
        return sum(np.absolute(self.Voce(params,xdata)-ydata))+addition

    def Gosh(self,params,xdata):
        '''stress(strain)=q*(e0+strain)**n-p'''
        xdata = np.array(xdata)
        return params[0]*(params[1]+xdata)**params[2]-params[3]
    def GoshJac(self,params,xdata,ydata):
        dfd0 = (params[1]+xdata)**params[2]
        dfd1 = params[0]*params[2]*(params[1]+xdata)**(params[2]-1)
        dfd2 = params[0]*(params[1]+xdata)**(params[2])*np.log(params[1]+xdata)
        dfd3 = -1.0*np.ones(len(xdata))
        return np.array([dfd0,dfd1,dfd2,dfd3]).T
    def GoshPrime(self,params,xdata):
        '''stress'(strain)=q*n*(e0+strain)**(n-1)'''
        xdata = np.array(xdata)
        return params[0]*params[2]*(params[1]+xdata)**(params[2]-1)
    def GoshObjFunc(self,params,xdata,ydata):
        xdata = np.array(xdata)
        ydata = np.array(ydata)
        addition = self.w*(np.absolute(self.GoshPrime(params,self.ue_strain)-self.Gosh(params,self.ue_strain))) 
        return sum(np.absolute(self.Gosh(params,xdata)-ydata))+addition
    def GoshFit(self):
        bounds=[(0,1e6),(0,1e-2),(0,1),(0,1e6)]
        res = opt.differential_evolution(self.GoshObjFunc,bounds,args=(self.xdata,self.ydata))
        #print 'Error:',sum(np.square(self.GoshObjFunc(res[0],self.xdata,self.ydata)))
        return res
    
    def Hockett_Sherby(self,params,xdata):
        '''stress(strain) = a+b*(1-exp(-c*strain**n))'''
        xdata = np.array(xdata)
        return params[0]+params[1]*(1-np.exp(-params[2]*xdata**(params[3])))
    def Hockett_SherbyPrime(self,params,xdata):
        '''stress'(strain) = b*c*n*strain**(n-1)*exp(-c*strain**n)'''
        xdata = np.array(xdata)
        return params[1]*params[2]*params[3]*xdata**(params[3]-1)*np.exp(-params[2]*xdata**(params[3]))
    def Hockett_SherbyObjFunc(self,params,xdata,ydata):
        xdata = np.array(xdata)
        ydata = np.array(ydata)
        addition = self.w*(np.absolute(self.Hockett_SherbyPrime(params,self.ue_strain)-self.Hockett_Sherby(params,self.ue_strain)))
        return sum(np.absolute(self.Hockett_Sherby(params,xdata)-ydata))+addition

    def SwiftVoce(self,params,xdata):
        '''stress(strain) = c*Swift+(1-c)*Voce'''
        xdata = np.array(xdata)
        return params[0]*(params[1]+xdata)**(params[2])*params[6]+\
               params[3]+params[4]*(1-np.exp(-params[5]*xdata))*(1-params[6])
    def SwiftVocePrime(self,params,xdata):
        xdata = np.array(xdata)
        return params[0]*params[2]*(params[1]+xdata)**(params[2]-1)*params[6]+\
               params[3]*params[4]*np.exp(-params[5]*xdata)*(1-params[6])
    def SwiftVoceObjFunc(self,params,xdata,ydata):
        xdata = np.array(xdata)
        ydata = np.array(ydata)
        addition = self.w*(np.absolute(self.SwiftVocePrime(params,self.ue_strain)-self.SwiftVoce(params,self.ue_strain)))
        return sum(np.absolute(self.SwiftVoce(params,xdata)-ydata))+addition

    def Modified_Zener_Hollomon(self,params,xdata,ydata,zdata):
        '''stress(strain,strainrate,temperature)=A*exp(Q/(R*T))*strainrate**m*
            (1+a*exp(-c*(strain-e)**2)*(1-b*exp(-N*strain**n)'''
        xdata = np.array(xdata) #temperature
        ydata = np.array(ydata) #strainrate
        zdata = np.array(zdata) #strain
        R = 8.314
        First = params[0]*np.exp(params[1]/(R*xdata))*ydata**(params[2])
        Second = 1+params[3]*np.exp(-params[4]*(zdata-params[5])**2)
        Third =  1-params[6]*np.exp(-params[7]*zdata**(params[8]))        
        return First*Second*Third
                  
    def Modified_Zener_HollomonPrime(self,params,xdata,ydata,zdata):
        pass
    def Modified_Zener_HollomonObjFunc(self,params,xdata,ydata,zdata,kdata):
        xdata = np.array(xdata)
        ydata = np.array(ydata)
        zdata = np.array(zdata)
        kdata = np.array(kdata)
        #return sum(np.absolute(self.Modified_Zener_Hollomon(params,xdata,ydata,zdata)-kdata))
        return self.Modified_Zener_Hollomon(params,xdata,ydata,zdata)-kdata
    def Modified_Zener_HollomonFit(self,xdata,ydata,zdata,kdata):
        # Modified_Zener_Hollomon Initial Params
        Init = [1,1,1,1,1,1,1,1,1]
        paramsFit = opt.leastsq(self.Modified_Zener_HollomonObjFunc,Init,args=(xdata,ydata,zdata,kdata))
        error = np.sqrt(sum(np.square(self.Modified_Zener_Hollomon(paramsFit[0],xdata,ydata,zdata)-kdata))/len(zdata))
        print 'Model:',"Modified_Zener_Hollomon"
        print 'Params:',paramsFit[0]
        print "Error:",error
        return [paramsFit[0],error]


##def main():
##    app = QtGui.QApplication(sys.argv)
##    f = reader.File('F:\pyform\pyform0\dc05_lank_l2.TRA')
##    f = File('F:\pyform\StahlRohdaten\St_WR_100_0_5_005_1.asc')
##    f = File('F:\pyform\Hui\St_WR_400_0_5_017.asc')
##    f = File('F:\pyform\Hui\data_0.asc')
##    fc = FC(f.strainT_exp(),f.stressT_exp(),temp=f.temp_exp(),rate=f.rate_exp())
##    res = fc.Modified_Zener_HollomonFit(f.temp_exp(),f.rate_exp(),f.strainT_exp(),f.stressT_exp())
##    Y = fc.values('Modified_Zener_Hollomon',res[0],f.strainT_exp())
    #ans = fc.params('Modified_Zener_Hollomon',[(-1e1,1e1),(-1e0,1e4),(-1e1,1e1),(-1e1,1e1),(-1e1,1e1),(0,1e2),(0,1e2),(0,1e4),(0,1e1)])
    #Y = fc.values('Modified_Zener_Hollomon',ans[0],f.strainT())
##    Y_prime = fc.SwiftPrime(ans[0],f.strainT())
##    pw = pg.plot()
##    pw.plot(f.strainT_exp(),Y,pen='r')
##    pw.plot(f.strainT_exp(),f.stressT_exp())
##    pw.plot(f.strainT(),fc.Swift(res.x,f.strainT()),pen='g')
##    pw.plot(f.strainT(),Y_prime)
##    ans = fc.params('Voce',[(0,1e4),(0,1e4),(0,1e3)])
##    fc.values('Voce',ans[0],2e-3)
##    ans = fc.params('Gosh',[(0,1e4),(0,1e-2),(0,1),(0,1e4)])
##    fc.values('Gosh',ans[0],2e-3)
##    ans = fc.params('Hockett_Sherby',[(0,1e4),(0,1e4),(0,1e3),(0,1)])
##    fc.values('Hockett_Sherby',ans[0],2e-3)
##    ans = fc.params('SwiftVoce',[(0,1e4),(0,1e-2),(0,1),(0,1e4),(0,1e4),(0,1e3),(0,1)])
##    fc.values('SwiftVoce',ans[0],2e-3)
    
##
##
##
##
##    app.exec_()
##    
##if __name__=='__main__':
##    main()
