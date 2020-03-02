
class Model:
	def __init__(self,\
		growth_rate = [0.18,0.1045],death_rate = 2*[0.0412],\
		capacity = [5e6,3e6],predation_rate = [2.2683e-7,3.422e-9],\
		conversion_rate = 5.32e-8):

		self.r1, self.r2 = growth_rate
		self.d1, self.d2 = death_rate
		self.k1, self.k2 = capacity
		self.a1, self.a2 = predation_rate
		self.b = conversion_rate

	def equilibrium_pts(self):
		
		#Finds equlibrium points as [M,N,Z]
		
		from numpy import array as np_array
		from numpy.linalg import solve as np_solve

		g = self.k2*(self.r2-self.d2)/self.r2

		A = np_array([[self.r1/self.k1,self.a1,0], [-self.a2,0,self.b],\
		[0,self.b,self.r2/self.k2]]);
		c = np_array([self.r1,self.d1,self.r2-self.d2]);

		self.E0,self.E1,self.E2,self.E3  = [3*[0],[self.k1,0,0],[0,0,g],\
		[self.k1,0,g]]
	
		self.E4 = [0,(self.b*self.r2*g-self.r2*self.d1)/(self.k2*self.b**2),\
					self.d1/self.b]

		self.E5 = np_solve(A,c).tolist()

	def bifurcation_pt(self):
		
		#calculates the bifurcation point of the system 
		
		from numpy import exp
		from scipy import optimize

		original_beta = self.b
		def D(beta):
			K1 = 1/self.k1
			K2 = 1/self.k2

			self.b = beta
			self.equilibrium_pts()
			C1 = exp(-self.r1*K1*self.E5[0])
			C2 = exp(-self.r2*K2*self.E5[2])
			
			p2 = -1-C1-C2
			p1 = C1 + C2 + C1*C2 + beta**2*self.E5[1]*(1- C2)/(self.r2*K2)\
			-self.a1*self.a2*self.E5[1]*(1-C1)/(self.r1*K1)
			p0 = -beta**2*self.E5[1]*C1*(1-C2)/(K2*self.r2) - C1*C2\
			 + self.a1*self.a2*self.E5[1]*C2*(1- C1)/(K1*self.r1)
			
			return 1 - p1 - p0**2 + p0*p2

		sol = optimize.root_scalar(D,method='brentq', bracket=[4e-8,1],x0=3e-7)

		self.bif = sol.root

		self.b = original_beta


	def calculate_alpha(self):
		
		#determines the alpha that satisfies the requirement that the tumor and resting cells grow at the same rate

		from scipy import optimize
		from numpy import asscalar
		def a1_bar(beta):
			K1 = 1/self.k1
			K2 = 1/self.k2
			return beta*self.r1*(K1*(-beta*self.r1+ self.d1*K2*self.r2)+K2*self.r2*self.a2)/\
				(self.d2*(beta*K1*self.r1-K2*self.r2*self.a2)+self.r2*\
				((K1*(-beta+self.d1*K2)*self.r1)+K2*self.r2*self.a2))

		b_opt = optimize.fmin(a1_bar,x0=3.8e-8,xtol=1e-10,ftol=1e-10,disp=False)
		self.a1 = asscalar(a1_bar(b_opt))

	def run_simulation(self,M0=2e5,N0=2e5,Z0=1.5e5,max_iter=2000):

		# runs the simulation starting with initial point (M0, N0, Z0)

		import numpy as np
		
		self.data = np.zeros([max_iter+1,3])

		self.data[0,:] = [M0,N0,Z0]

		for n in range(max_iter):
			c1 = self.r1-(self.a1*self.data[n,1])
			c2 = self.r1/self.k1
			c3 = self.r2-self.b*self.data[n,1]-self.d2
			c4 = self.r2/self.k2
			self.data[n+1,0] = self.data[n,0]*c1/((c1-c2*self.data[n,0])*np.exp(-c1)+c2*self.data[n,0])
			self.data[n+1,1] = self.data[n,1] * np.exp(self.b*self.data[n,2]-self.d1-self.a2*self.data[n,0])
			self.data[n+1,2] = self.data[n,2]*c3/((c3-c4*self.data[n,2])*np.exp(-c3)+c4*self.data[n,2])

	def draw_vals(self, ax0):
		
		# updates plot with new simulation values when the beta value is updated
		
		ax0.cla()
		self.chart, = ax0.plot(self.data[:,0].T,self.data[:,1])
		ax0.set_xlim([0,self.k1])
		ax0.set_ylim([0,1e6])
		ax0.set_xlabel('Tumor Cell Population (M)')
		ax0.set_ylabel('Hunting Cell Population (N)')
		ax0.grid('on')
		ax0.ticklabel_format(style='sci',scilimits=(0,0))

	def update_beta(self,new_beta):
		
		# update beta when the slider is moved in the interactive plot 
		
		self.b = new_beta*1e-8
		self.run_simulation()
		self.chart.set_xdata(self.data[:,0])
		self.chart.set_ydata(self.data[:,1])
	
	def plot_data(self):

		import matplotlib.pyplot as plt
		from matplotlib.widgets import Slider

		fig = plt.figure()
		ax0 = plt.axes([0.085,0.2,0.9,0.75])
		ax1 = plt.axes([0.21, 0.02, 0.7, 0.03])

		self.draw_vals(ax0)

		B = Slider(ax1,label='Beta value (1e-8)', valmin=4,valmax=50,valstep=0.1)
		Slider.set_val(B,10)
		
		Slider.on_changed(B,self.update_beta)
		plt.show()


	def plot_data_3D(self):

		# plots the values in time in a 3D plot

		import matplotlib.pyplot as plt
		from mpl_toolkits.mplot3d import Axes3D

		fig = plt.figure()
		ax = Axes3D(fig)
		ax.plot(p.data[:,0],p.data[:,1],p.data[:,2])
		ax.set_xlabel('M')
		ax.set_ylabel('N')
		ax.set_zlabel('Z')
		ax.ticklabel_format(style='sci',scilimits=(0,0))
		plt.show()

	def animate_plots(self):

		# generates a gif of the plot as it changes with the value of beta 

		from matplotlib import pyplot as plt
		from matplotlib.animation import FuncAnimation, PillowWriter
		import numpy as np

		frame_count = 100
		fig, ax = plt.subplots()

		self.draw_vals(ax)
		
		b_range = np.linspace(4,60,frame_count)
		
		def animate(i):
			self.update_beta(b_range[i])
			return [self.chart]

		anim = FuncAnimation(fig, animate, frame_count, blit=True)
		anim.save("bifurcation.gif", writer=PillowWriter(fps=24))
		plt.show()

p = Model()
p.equilibrium_pts()
p.bifurcation_pt()
p.calculate_alpha()
p.run_simulation()
p.animate_plots()
