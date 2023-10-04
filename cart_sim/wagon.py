import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Define statespace function
def f(t, u, m_s, m_v, l, a, f0,g):
    i_rv = 0.01
    i_b = 1/12*m_s*(l**2)+m_s*((l/2-a)**2)

    r2 = u[0]*u[3]**2-g*np.sin(u[1])

    q1 = f0*a-(l/2-a)*m_s*g*np.cos(u[1])
    q2 = (u[0]*m_v*r2*np.cos(u[1]))/np.sin(u[1])
    q3 = 2*u[0]*m_v*u[2]*u[3]
    q4 = (u[0]**2)*m_v*(u[3]**2)*np.cos(u[1])
    q5 = np.sin(u[1])
    q6 = i_rv+i_b+u[0]**2*m_v
    phi2 = (q1+q2-q3-(q4)/q5)/(q6)
    dudt = [u[2], u[3], r2,phi2] # [r1,phi1,r2,phi2]
    return dudt


# Define time spans, initial values, and constants
g = 9.82
a = 0.25
l = 1
m_s = 5
m_v = 1
f0 = 200
tspan = np.linspace(0, 1.4, 1000)
#u0 = [0, 0, 0.3, 3, 3, 0]
u0 = [0.3*l,-np.pi/32,0,0] #[r,phi,r1,phi1]


def off_beam(t,u): return l-a-u[0]
off_beam.terminal =True
off_beam.direction = -1

# Solve differential equation



sol = solve_ivp(lambda t, u: f(t, u, m_s, m_v, l, a, 2*f0, g),
                [tspan[0], tspan[-1]], u0, t_eval=tspan, rtol=1e-5, events=off_beam )
# Plot states
plt.figure(1)
for i in range(sol.y.shape[0]):
    plt.plot(sol.t, sol.y[i])
plt.legend(['r', 'phi', 'r\'', 'phi\''])
plt.xlabel('t')
plt.ylabel('tillstånd')

print(sol.y[0][-1]) # r last
print(sol.y[1][-1]) # phi last

#plt.show(block=True)

def f2 (t, u):
  y2 = -g
  dudt = [u[2], u[3], 0,y2] # [x',y',x'',y'']
  return dudt

def hit_ground(t,u): return u[1]
hit_ground.terminal =True
hit_ground.direction = -1
tspan2 = np.linspace(0, 4.4, 1000)


#start values trajectory
last_r = sol.y[0][-1]
last_phi = sol.y[1][-1]

last_r_v = sol.y[2][-1]

print(last_r,last_phi,last_r_v)
u1 = [last_r*np.cos(last_phi),last_r*np.sin(last_phi),last_r_v*np.cos(last_phi),last_r_v*np.sin(last_phi)] #[x,y,x',y']
print(u1)


sol2 = solve_ivp(lambda t, u: f2(t, u),
                [tspan2[0], tspan2[-1]], u1, t_eval=tspan2, rtol=1e-5, events=hit_ground )


# plot trajectory
plt.figure(2)
plt.plot(sol2.y[0], sol2.y[1])
xLabel = 'x '+str(sol2.y[0][-1])
print(xLabel)
plt.xlabel(xLabel)
plt.ylabel('y')
plt.show(block=True)
