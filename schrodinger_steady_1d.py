#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import argparse
import hermite

# class for schrodinger equation
class scheq :
    # member value
    def __init__(self):
        self.delta_x=0.05
        self.xa =1 # 井戸の境界。 x=±xaに井戸の端がある。
        self.eps_E=0.005 # 収束条件
        self.nn=5 #両端から積分を行う際の両端の座標をxa(と-xa)の何倍にするかを示すパラメータ
        self.xL0, self.xR0  = -self.nn*self.xa, self.nn*self.xa
        self.Nx =  int((self.xR0-self.xL0)/self.delta_x)
        self.delta_x = (self.xR0-self.xL0)/float(self.Nx)
        #uLとuRの一致具合をチェックする位置のインデックス。井戸の境界に選んでいる。
        self.i_match = int((self.xa-self.xL0)/self.delta_x)
        self.nL = self.i_match
        self.nR = self.Nx-self.nL
        self.uL = np.zeros([self.nL],float)
        self.uR = np.zeros([self.nR],float)
        self.E=np.pi**2/4
        self.xL=np.zeros([self.Nx])
        self.xR=np.zeros([self.Nx])
        for i in range(self.Nx):
            self.xL[i] = self.xL0 + i*self.delta_x
            self.xR[i] = self.xR0 - i*self.delta_x
        self.k2L=np.zeros([self.Nx+1])
        self.k2R=np.zeros([self.Nx+1])
        self.Elis = list()
        self.check_Elist = list()
        self.Solved_Eigenvalue = list()
        self.EEmin = 0.0
        self.EEmax = 0.0
        self.delta_EE=0.0
        self.mass = 0.5

    # member func
    def well_potential(self,x): # 井戸型ポテンシャル
        if np.abs(x) >  self.xa :
            v = 1000.0
        else :
            v = 0
        return v

    def harmonic_potential(self,x):
        return (-self.mass**2*x**2)

    # set a cofficient of 1th-order psi
    def setk2(self,E): # for E<0
        for i in range(self.Nx+1):
            xxL = self.xL0 + i*self.delta_x
            xxR = self.xR0 - i*self.delta_x
            self.k2L[i] = E-self.harmonic_potential(xxL)
            self.k2R[i] = E-self.harmonic_potential(xxR)

    # set boundary condtions
    def set_condition(self):
        self.uL[0] = 0.0
        self.uL[1] = 1e-12
        self.uR[0] = 0.0
        self.uR[1] = 1e-12

    # set boundary condtions
    def set_condition_even(self):
        self.uL[0] =  0.0
        self.uR[0] =  0.0
        self.uL[1] =  1e-12
        self.uR[1] =  1e-12

    # set boundary condtions
    def set_condition_odd(self):
        self.uL[0] =  0.0
        self.uR[0] =  0.0
        self.uL[1] = -1e-12
        self.uR[1] =  1e-12

    # numerov method
    def Numerov(self,N,delta_x,k2,u):
        b = (delta_x**2)/12.0
        for i in range(1,N-1):
            u[i+1] = (2*u[i]*(1-5*b*k2[i])-(1+b*k2[i-1])*u[i-1])/(1+b*k2[i+1])

    def plot_potential(self):
        # ポテンシャル関数のプロット
        XXX= np.linspace(xL0,xR0, Nx)
        POT=np.zeros([Nx])
        for i in range(Nx):
            POT[i] = V(xL0 + i*delta_x)
        plt.xlabel('X (Bohr)') # ｘ軸のラベル
        plt.ylabel('V (X) (Ry)') # y軸のラベル
        plt.hlines([E], xL0,xR0, linestyles="dashed")  #Energy
        plt.plot(XXX,POT,'-',color='blue')
        plt.show()

    def plot_k2(self):
        # k^2(x)のプロット
        XXX= np.linspace(xL0,xR0, Nx+1)
        plt.plot(XXX, k2L,'-')
        plt.show()

    def E_eval(self):
        uLdash = (self.uL[-1]-self.uL[-2])/self.delta_x
        uRdash = (self.uR[-2]-self.uR[-1])/self.delta_x
        logderi_L=  uLdash/self.uL[-1]
        logderi_R=  uRdash/self.uR[-1]
        return (logderi_L- logderi_R)/(logderi_L+logderi_R)

    def normarize_func(self,u):
        factor = ((self.xR0-self.xL0)/self.Nx)*(np.sum(u[1:-2]**2))
        return factor

    def plot_eigenfunc(self,color_name):
        uuu=np.concatenate([self.uL[0:self.nL-2],self.uR[::-1]],axis=0)
        XX=np.linspace(self.xL0,self.xR0, len(uuu))

        factor=np.sqrt(self.normarize_func(uuu))

        plt.plot(XX,uuu/factor,'-',color=color_name,label='Psi')
        plt.plot(XX,(uuu/factor)**2,'-',color='red',label='| Psi |^2')
        plt.xlabel('X (Bohr)') # ｘ軸のラベル
        plt.ylabel('') # y軸のラベル
        plt.legend(loc='upper right')
        plt.show()

    def clean(self):
        self.uL = np.zeros([self.nL],float)
        self.uR = np.zeros([self.nR],float)

    # boundary condtion (odd or even function)
    def solve(self,eo):
        self.EEmin = 0.1
        self.EEmax = 20
        self.delta_EE=0.01

        NE = int((self.EEmax-self.EEmin)/self.delta_EE)
        self.Elis = list()
        self.Solved_Eigenvalue = list()
        self.check_Elis = list()
        for i in range(NE+1):
            EE=self.EEmin+i*(self.EEmax-self.EEmin)/NE

            self.nL = self.i_match
            self.nR = self.Nx-self.nL

            self.clean()

            if eo == "odd" :
                self.set_condition_odd()
            elif  eo == "even" :
                self.set_condition_even()
            self.setk2(EE)

            self.Numerov(self.nL,self.delta_x,self.k2L,self.uL)
            self.Numerov(self.nR,self.delta_x,self.k2R,self.uR)

            a1= self.E_eval()
            if a1 :  # a1がTrueのとき
                self.Elis.append(EE)
                self.check_Elis.append(a1)
                if np.abs(a1) <= self.eps_E :
                    print("Eigen_value = ", EE)
                    self.Solved_Eigenvalue.append(EE)
                    self.plot_eigenfunc("blue")

    def plot_delta_func(self,c):
        plt.plot(self.Elis, self.check_Elis, 'o',markersize=3, color=c,linewidth=1)
        plt.grid(True) #グラフの枠を作成
        plt.xlim(self.EEmin, self.EEmax) # 描くxの範囲を[xmin,xmax]にする
        plt.ylim(-10, 10) # 描くyの範囲を[ymin,ymax]にする
        plt.hlines([0], self.EEmin,self.EEmax, linestyles="dashed")  # y=y1とy2に破線を描く
        plt.xlabel('Energy (Ry)') # ｘ軸のラベル
        plt.ylabel('Delta_E_function') # y軸のラベル
        plt.show()


if __name__ == "__main__":
    qm = scheq()
    print("E= ",qm.E)
    print("xL0,xR0, i_match, delta_x=", qm.xL0, qm.xR0, qm.i_match, qm.delta_x)
    print("Nx, nL,nR=",qm.Nx, qm.nL,qm.nR)

    #qm.set_condition()
    qm.setk2(qm.E)

# search solution
    qm.solve("even")
    #qm.plot_delta_func("blue")
    qm.solve("odd")
    #qm.plot_delta_func("red")
