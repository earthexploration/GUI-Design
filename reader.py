#!/usr/bin/python
# -*- coding: UTF-8 -*-
# read .TRA
# $env 
# $dat/exp

# ??? Unit
import numpy as np
import codecs
import csv
import math

class File:
    def __init__(self,filename):
        col = 40
        step = 5
        row = 0
        dol = 0

        env = {}
        exp = {}
        #dat = {}
        dat_s = {} # string
        dat_f = {} # float
        env_n = 0
        dat_n = 0
        
        
        self.exp = None
        self.env = None
        self.exp_n = None
        self.env_n = None
        #self.dat = None
 
        with codecs.open(filename,'rb',encoding='iso-8859-1') as f:
            for line in f.readlines():
                row = row + 1
                if line[0] == '$':
                    dol = dol + 1
                    continue
                # dol=1/dol%2!=0
                if dol==1:
                    name = line[0:col].strip()
                    value = line[col:].split(' ')[0].strip()
                    env.setdefault(name,value)
                    env_n = env_n + 1
                if dol%2 == 0:
                    if ';' in line:
                        elements = line.split(';')
                    elif ',' in line:
                        elements = line.split(',')
                    else:
                        elements = line.split()
                    #print elements
                    if len(elements) >= 2: # at least (x,y)
                        #print elements
                        try:
                            for i in range(len(elements)):
                                    dat_f.setdefault(i,[]).append(float(elements[i].strip()))
                            dat_n = dat_n + 1
                        except ValueError:
                            p = 0
                            el = []
                            for i in elements:
                                if u'(' in i:
                                    p = p + 1
                                if u')' in i:
                                    p = p - 1
                                    continue
                                if p%2 == 0:
                                    el.append(i)
                            elements = el
                            #print elements
                            for i in range(len(elements)): 
                                dat_s.setdefault(i,[]).append(elements[i].strip())
                        except IOError:
                            print "Fail to read file"
                            
        if len(dat_f.keys())!=len(dat_s.keys()):
            print "specify colomun header!"
        else:
            for i in range(len(dat_f.keys())):
                f = [k for j,k in enumerate(dat_f[i]) if j%step==0]
                exp.setdefault(dat_s[i][0],f)
        exp_n = dat_n/step
        #print env.keys()
        #print exp.keys()
        self.exp = exp
        self.env = env
        self.exp_n = exp_n
        self.env_n = env_n
        #print "row:",row
        #print "dat_n:",dat_n
        #print "exp_n:",exp_n
        #print "exp:",len(exp[u'Zeit'])
        #print dat.keys()
        #self.dat = dat
           
    # env
    def width_env(self):
        if u'Probenbreite' in self.env.keys():
            return self.env[u'Probenbreite']
        else:
            return None
    def thickness_env(self):
        if u'Probendicke' in self.env.keys():
            return self.env[u'Probendicke']
        else:
            return None
    def length_env(self):
        if u'Anfangsmesslänge' in self.env.keys():
            return self.env[u'Anfangsmesslänge']
        else:
            return None
    def speed_env(self):
        if u'Geschwindigkeit im Fließbereich' in self.env.keys():
            return self.env[u'Geschwindigkeit im Fließbereich']
        else:
            return None
    def Emodulus_env(self):
        if u'Elastizitätsmodul' in self.env.keys():
            return self.env[u'Elastizitätsmodul']
        else:
            return None
    def rx1_env(self):
        if u'Senkrechte Anisotropie r{lo x1}' in self.env.keys():
            return self.env[u'Senkrechte Anisotropie r{lo x1}']
        else:
            return None
    def rx2_env(self):
        if u'Senkrechte Anisotropie r{lo x2}' in self.env.keys():
            return self.env[u'Senkrechte Anisotropie r{lo x2}']
        else:
            return None
    def rx3_env(self):
        if u'Senkrechte Anisotropie r{lo x3}' in self.env.keys():
            return self.env[u'Senkrechte Anisotropie r{lo x3}']
        else:
            return None
    def rx4_env(self):
        if u'Senkrechte Anisotropie r{lo x4}' in self.env.keys():
            return self.env[u'Senkrechte Anisotropie r{lo x4}']
        else:
            return None
    def rx5_env(self):
        if u'Senkrechte Anisotropie r{lo x5}' in self.env.keys():
            return self.env[u'Senkrechte Anisotropie r{lo x5}']
        else:
            return None
    def rx6_env(self):
        if u'Senkrechte Anisotropie r{lo x6}' in self.env.keys():
            return self.env[u'Senkrechte Anisotropie r{lo x6}']
        else:
            return None
    def rx12_env(self):
        if u'Senkrechte Anisotropie r{lo x1-x2/dynami' in self.env.keys():
            return self.env[u'Senkrechte Anisotropie r{lo x1-x2/dynami']
        else:
            return None
        
    def init_area(self):
        # width,thickness Unit: mm
        if self.width_env() and self.thickness_env() != None:
            A = float(self.width_env())*float(self.thickness_env())*1e-4
            return A
        else:
            return None
    # exp
    def time_exp(self):
        if u'Zeit' in self.exp.keys():
            return self.exp[u'Zeit']
        else:
            return None
        
    def displacement_exp(self):
        if u'Standardweg' in self.exp.keys():
            return [float(d) for d in self.exp[u'Standardweg']]
        else:
            return None
    def displacement_exp_min(self):
        if self.displacement_exp !=None:
            return min(self.displacement_exp())
    def displacement_exp_max(self):
        if self.displacement_exp !=None:
            return max(self.displacement_exp())      
        
    def force_exp(self):
        if u'Standardkraft' in self.exp.keys():
            return [float(f) for f in self.exp[u'Standardkraft']]
        else:
            return None
        
    def strainE_exp(self):
        if u'Dehnung' in self.exp.keys():
            return [float(e) for e in self.exp[u'Dehnung']]
        else:
            return None
    def strainE_exp_min(self):
        if self.strainE_exp != None:
            return min(self.strainE_exp())
    def strainE_exp_max(self):
        if self.strainE_exp != None:
            return max(self.strainE_exp())
        
    def stressE_exp(self):
        if self.force_exp() and self.init_area() != None:
            A0 = self.init_area()
            f = [float(i)/A0 for i in self.force_exp()]
            return f
    def stressE_exp_min(self):
        if self.stressE_exp != None:
            return min(self.stressE_exp())
    def stressE_exp_max(self):
        if self.stressE_exp != None:
            return max(self.stressE_exp())
        
    def strainT_exp(self):
        if u'wahre Dehnung' in self.exp.keys():
            return [float(s) for s in self.exp[u'wahre Dehnung']]
        elif u'Phi' in self.exp.keys():
            return [float(s) for s in self.exp[u'Phi']]
        else:
            return None
    def strainT_exp_min(self):
        if self.strainT_exp !=None:
            return min(self.strainT_exp())
    def strainT_exp_max(self):
        if self.strainT_exp !=None:
            return max(self.strainT_exp())
        
    def stressT_exp(self):
        if u'wahre Spannung' in self.exp.keys():
            return [float(a) for a in self.exp[u'wahre Spannung']]
        elif u'kf' in self.exp.keys():
            return [float(a) for a in self.exp[u'kf']]
        else:
            return None
    def width_exp(self):
        if u'Probenbreite' in self.exp.keys():
            return self.exp[u'Probenbreite']
        else:
            return None
    def widthA_exp(self):
        if u'Breitenänderun' in self.exp.keys():
            return self.exp[u'Breitenänderun']
        else:
            return None
    def log_exp(self):
        if u'logarithmische' in self.exp.keys():
            return self.exp[u'logarithmische']
        else:
            return None
    def Ani_exp(self):
        if u'senkrechte Ani'in self.exp.keys():
            return self.exp[u'senkrechte Ani']
        else:
            return None
    def temp_exp(self):
        if u'Tem.1' in self.exp.keys():
            return [float(t) for t in self.exp[u'Tem.1']]
        elif u'Tem' in self.exp.keys():
            return [float(t) for t in self.exp[u'Tem']]
        else:
            return None
    def rate_exp(self):
        if u'Phi*' in self.exp.keys():
            return [float(r) for r in self.exp[u'Phi*']]
        else:
            return None
    # Output
        
    def Ec_find(self,Ec=2e-1):
        # Ec 0.2% elastic
        Ec=2e-1
        if self.strainE_exp() != None:
            for i,j in enumerate(self.strainE_exp()):
                j = float(j)
                if j > Ec:
                    return i
                
        elif self.strainT_exp != None:
            # math.log(1+Ec)
            for i,j in enumerate(self.strainT_exp()):
                j = float(j)
                if j > Ec:
                    return i
        else:
            return None

    def Fc_find(self):
        # max stress
        s = []
        if self.force_exp() != None:
            for i in self.force_exp():
                s.append(float(i))
            return s.index(max(s))
        elif self.strainT_exp() and self.stressT_exp() != None:
            for i in range(len(self.strainT_exp())):
                stressT = float(self.stressT_exp()[i])
                strainT = float(self.strainT_exp()[i])
                # stressT = stressE*exp(strainT)
                s.append(stressT*1.0/math.exp(strainT))
            return s.index(max(s))
        else:
            return None

    def strainE(self):
        E = [float(i) for i in self.strainE_exp()[self.Ec_find():self.Fc_find()+1]]
        return E
    def stressE(self):
        A0 = self.init_area()
        F = [float(i)/A0 for i in self.force_exp()[self.Ec_find():self.Fc_find()+1]]
        return F
    def strainT(self):
        E = [float(i) for i in self.strainT_exp()[self.Ec_find():self.Fc_find()+1]]
        return E
    def stressT(self):
        F = [float(i) for i in self.stressT_exp()[self.Ec_find():self.Fc_find()+1]]
        return F
    def displacement(self):
        u = [float(i) for i in self.displacement_exp()[self.Ec_find():self.Fc_find()+1]]
        return u
    def force(self):
        F = [float(i) for i in self.force_exp()[self.Ec_find():self.Fc_find()+1]]
        return F
    def temp(self):
        T = [float(i) for i in self.temp_exp()[self.Ec_find():self.Fc_find()+1]]
        return T
    def rate(self):
        R = [float(i) for i in self.rate_exp()[self.Ec_find():self.Fc_find()+1]]
        return R

##def main():
##    f = File('F:\pyform\Hui\St_WR_100_0_5_005.asc')
##    f = File('F:\pyform\Hui\St_WR_20_0_5_001.asc')
##    f = File('F:\pyform\pyform1\dc05_lank_d1.TRA')
##    f = File('F:\pyform\Hui\data.asc')
##    print np.shape(f.strainT_exp())
##    print np.shape(f.stressT_exp())
##    print f.strainT()[0]
##    print f.force_exp()[0]
##    print f.displacement_exp()[0]
##    print type(f.displacement()[0])
##    print f.force()[0]
##    print f.displacement_exp_max()
##    print f.strainE_exp_max()
##    print f.strainT_exp_max()
##    print len(f.stressT())
##    print len(f.strainT())
##    print f.Ec_find()
##    print f.Fc_find()
##    print np.shape(f.temp_exp())
##    print np.shape(f.rate_exp())
##    print type(f.stressE_exp()[0])
##if __name__=='__main__':
##    main()
        
                    
                
            
