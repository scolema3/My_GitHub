def symop(var_dict):

	if var_dict['symmetry'] == True:
		if var_dict['pickle_jar'] == NULL:
			reflections = [np.array([[1,0,0],[0,1,0],[0,0,1]]),\
np.array([[-1,0,0],[0,1,0],[0,0,1]]),\
np.array([[1,0,0],[0,-1,0],[0,0,1]]),\
np.array([[1,0,0],[0,1,0],[0,0,-1]]),\
np.array([[1,0,0],[0,-1,0],[0,0,-1]]),\
np.array([[-1,0,0],[0,1,0],[0,0,-1]]),\
np.array([[-1,0,0],[0,-1,0],[0,0,1]]),\
np.array([[-1,0,0],[0,-1,0],[0,0,-1]])]
			def rot_x(Theta):
				r_x = np.array([[1,0,0],[0,np.cos(theta),-1*np.sin(theta)],[0,np.sin(theta),np.cos(theta)]])
				return r_x 
			def rot_y(Theta):
				r_y = np.array([[np.cos(theta),0,np.sin(theta)],[0,1,0],[-1*np.sin(theta),0,np.cos(theta)]])
				return r_y
			def rot_z(Theta):
				r_z = np.array([[np.cos(theta),-1*np.sin(theta),0],[np.sin(theta),np.cos(theta),0],[0,0,1]])
				return r_z
			thet = 2*np.pi*np.array([-1,-2,-3,-4,-6,1,2,3,4,6])
			for i in range(0,len(reflections)):
				for j in range(0,len(thet)):
					for k in range(0,len(thet)):
						for m in range(0,len(thet)):
							symtest=np.dot(reflection[i],rot_x(thet[j]))
							symtest = np.dot(symtest,rot_y(thet[k]))
							symtest = np.dot(symtest,rot_z(thet[m]))
							compare = np.zeros((3,3))
							compare[1,:] = symtest[1,:]
							compare[2,:] = symtest[2,:]
							compare[3,:] = symtest[3,:]
							if np.sum((np.eye(3,3)-compare)**2) == 0:
								sym = np.concatenate((np.eye(3,3),symtest),axis = 0)
			var_dict['symops'] = sym
		else:
			#load pickpled file
		max_n = 3
		n = 2 * max_n + 1
		R = np.flilr(np.indices((n,n,n)).T.reshape(-1,3) - max_n) 
		R = np.concatenate((np.zeros(len(R),1),R),axis = 1)
		expand = np.empty((0,4))
		for i in range(0,len(var_dict['basis'])):
			expand = np.concatenate((expand,R+var_dict['basis'][i,:]))





	else:
		
	
