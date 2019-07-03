from scipy.optimize import minimize
from scipy import interpolate
import matplotlib.pyplot as plt
import numpy as np
e = 0.003
theta=0.36
eta=2
alpha=.3
delta=0.1
beta=0.98
#The Discrete Space
n_dis=100
K_t=np.linspace(0.001,6,n_dis)
Value_fun=list(np.ones(n_dis))
Value_fun=interpolate.interp1d(K_t,Value_fun)
Policy=[]
#plt.plot(K_t,Value_fun(K_t))
New_Value_best=[]
for k in range(301):
    for i in range(n_dis):
        #fun = lambda K_t1 : (
        #                -np.log(K_t[i]**theta-K_t1+(1-delta)*K_t[i])
        #                -beta*Value_fun(K_t1)
        #                )
        fun=lambda K_t1:-(
                (K_t[i]**alpha-K_t1+(1-delta)*K_t[i])**(1-eta)/(1-eta)
                +beta*Value_fun(K_t[i])
                )
        cons = (
                {'type': 'ineq', 'fun': lambda K_t1: K_t1-e}
               )
        x0 = np.array((0.1)) # 设置初始值
        res = minimize(fun, x0, method='SLSQP', constraints=cons)
        Policy.append(res.x[0])
        New_Value_best.append(-res.fun)
    if New_Value_best[0]==np.nan:
        New_Value_best[0]=New_Value_best[1]
    Value_fun_new=interpolate.interp1d(K_t,New_Value_best)
    if k%2==0:
        plt.subplot(211)
        plt.plot(K_t,Value_fun(K_t))
    if abs(Value_fun_new(3)-Value_fun(3))<0.05:
        plt.figure()
        plt.plot(K_t,Value_fun(K_t))
        plt.figure()
        plt.plot(K_t,Policy)
        print(k)
        break
    plt.subplot(212)
    plt.plot(K_t,Policy)
    Policy.clear()
    Value_fun=Value_fun_new
    New_Value_best.clear()
