"""
Lawrence Stewart - ls3914@ic.ac.uk CID=00948972
Assumes cost.f90 has been compiled with f2py to generate the module,
hw2mod.so (filename may also be of form hw2mod.xxx.so where xxx is system-dependent text) using
"""
import numpy as np
import matplotlib.pyplot as plt
from hw2mod import cost
from hw2mod import hw2
import time

#my imports
from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
from scipy.optimize import minimize

#2nd attemptfrom mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D

def visualize(Nx,Ny):
    """Display cost function with and without noise on an Ny x Nx grid
    to see lots of noise set c.noise_amp to 100
    """
    #Set noise amplitude for noisey plot
    cost.c_noise_amp=1

    #Generate arrays for -Nx to Nx and -Ny to Ny
    xvec=np.arange(-Nx,Nx+1,1) #Nx+1 to include Nx in the arange 
    yvec=np.arange(-Ny,Ny+1,1)

    #Create lattice of points
    X,Y=np.meshgrid(xvec,yvec)

   #Stack the points as arrays (i.e generate all pairs)
    Z=np.dstack([X,Y])

    #Create the cost matrix from applying cost.costj to Z and also with noise
    costs=np.zeros(np.shape(X))
    noisey_costs=np.zeros(np.shape(X))
    for i in range(len(Z)):
        for j in range(len(Z[i])):
            costs[i][j]=cost.costj(Z[i][j])

            #add noise
            cost.c_noise = True
            #Calculate noisey cost fucntion
            noisey_costs[i][j]=cost.costj(Z[i][j])
            #Turn noise of 
            cost.c_noise = False


    #Set the noise to false after:
    cost.c_noise=False
    

    #Plot a contour for no noise
    plt.figure(figsize=(14,6))
    plt.suptitle('Lawrence Stewart - Created Using visualize().')
    plt.subplot(121)
    plt.xlabel("x")
    plt.ylabel("y")
    CS = plt.contour(X, Y, costs)
    plt.clabel(CS, inline=1, fontsize=10)
    plt.title('Contour Plot of Noiseless Cost Function') 


    #Noise
    plt.subplot(122)
    CS = plt.contour(X, Y, noisey_costs)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.clabel(CS, inline=1, fontsize=10)
    plt.title('Contour Plot for Noisy Cost Function , Amplitude = %s '%cost.c_noise_amp) 
    plt.show() 
  
    #Plot figure without noise

    # fig = plt.figure(figsize=(9,6))
    # ax = fig.gca(projection='3d')
    # surf = ax.plot_surface(X, Y, costs, cmap='magma',linewidth=0, antialiased=False) #cm.coolwarm
    # fig.colorbar(surf, shrink=0.5, aspect=5)
    # plt.title("Noiseless Cost Function Surface Plot")
    # plt.xlabel("y")
    # plt.ylabel("x")
    # fig.autofmt_xdate()
    # ax.set_zlabel("Cost(x,y)")
    # plt.grid('on')
    
    # #Plot figure with noise
    # fig = plt.figure(figsize=(9,6))
    # ax = fig.gca(projection='3d')
    # surf = ax.plot_surface(X, Y, costs, cmap='magma',linewidth=0, antialiased=False) #cm.coolwarm
    # fig.colorbar(surf, shrink=0.5, aspect=5)
    # plt.title("Cost Function Surface Plot with amplitude = %i" %(cost.c_noise_amp))
    # plt.xlabel("y")
    # plt.ylabel("x")
    # fig.autofmt_xdate()
    # ax.set_zlabel("Cost(x,y)")
    # plt.grid('on')
    # plt.show()

    
    return None



def newton_test(xg,display=False):
    """ Use Newton's method to minimize cost function defined in cost module
    Input variable xg is initial guess for location of minimum. When display
    is true, a figure illustrating the convergence of the method should be
    generated

    Output variables: xf -- computed location of minimum, jf -- computed minimum
    Further output can be added to the tuple, output, as needed. It may also
    be left empty.
    
    NOTE- FOR PLOTTING PURPOSES PLEASE INPUT xg as a list, so the inbuilt latex will work in plot titles--
    """
    output = ()
    dat=hw2.newton(xg)
    xf=dat[0]
    jf=dat[1]

    if display==True:

        #Generate current distance from the minimum
        xpath=hw2.xpath
        distances=[]
        for i in range(len(xpath)):
            temp=[1,1]-xpath[i]
            distances.append(np.sqrt(temp[0]**2+temp[1]**2))

        plt.figure(figsize=(14, 7))  
        plt.suptitle('Lawrence Stewart - Created Using newton_test().')

        #plot the cost function at each point
        plt.subplot(121)
        plt.plot(np.arange(1,len(hw2.jpath)+1,1),hw2.jpath,alpha=0.8,color='r')
        plt.xlabel("Iteration")
        plt.ylabel("Cost")
        ax = plt.gca()
        ax.set_facecolor('#D9E6E8')
        plt.title("Cost at each Iteration of Netwons Method, xg =%s"%xg)
        plt.grid('on')
  

        plt.subplot(122)
        plt.plot(np.arange(1,len(hw2.jpath)+1,1),distances,alpha=0.8,color='r')
        plt.title("Distance from Minimum at Each Iteration, xg =%s"%xg)
        plt.xlabel("Iteration")
        plt.ylabel("Distance")
        ax = plt.gca()
        ax.set_facecolor('#D9E6E8')
        plt.grid('on')
        plt.show()




    return xf,jf,output


def bracket_descent_test(xg,display=False):
    """ Use bracket-descent to minimize cost function defined in cost module
    Input variable xg is initial guess for location of minimum. When display
    is true, 1-2 figures comparing the B-D and Newton steps should be generated

    Output variables: xf -- computed location of minimum, jf -- computed minimum


    Discussion and Explanation:

    Bracket descent and newtons test operate in very different fashions as seen in the figures attatched. hw231.png shows the step size, h,  
    that is taken at each iteration of both the algorithms. Distance is defined with the usual euclidean norm, and one can see that Newtons 
    method takes initally very large steps, approximately 10000 (for the starting point [100,10]). The step size drastically decreases
    as the algorithm converges upon the minimum. The reason for such a high initial step size is the gradient based framework that
    Newtons method operates from, allowing it to initially move in large steps towards the minimum. Bracket Descent remains approximately 
    constant in the step size, which is to be expected due to the triangular method of descent that the algorithm utilizes, (the algorithm
    is bounded in the step size it can do).

    h232.png shows the wallclock  and CPU time for both of the methods. Due to a faster convergance, Newtons method terminates after
    a shorter duration, for both CPU and Wallclock time; approximately 0.00001 wallclock and roughly the same for CPU. Bracket descent
    takes longer to converge, with approximately 0.00005 for wallclock and CPU time. If we were to take points that were further away
    from the global minimum of [1,1], we would see this result extrapolated, due to the constant nature of the stepsize of Bracket 
    Descent. The size of newton's tests intial movements would increase with a further away starting point, and the time would remain
    small.

    """


    if display==False:
        xf,jf=hw2.bracket_descent(xg)



    if display==True:
         
        average_bd_wallclock=0
        average_newton_wallclock=0
        average_newton_cpu_time=0
        average_bd_cpu_time=0

        #Time over an average:
        for i in range(20):

            t1=time.time() #start timer 1
            tp1=time.process_time() #start timer 2

            #Run bracket descent
            xf,jf=hw2.bracket_descent(xg)

            t2 = time.time() #t2-t1 gives wallclock time
            tp2 = time.process_time() #tp2-tp1 gives cpu time -- depends on number of cores!

            bd_wallclock=t2-t1
            bd_cpu_time=tp2-tp1
        
            average_bd_wallclock+=bd_wallclock
            average_bd_cpu_time+=bd_cpu_time




        xpath=hw2.xpath


        bracket_steps=[xpath[i+1]-xpath[i] for i in range(len(xpath)-1)]
        newton_steps=[hw2.newtonstep(xpath[i]) for i in range(len(xpath)-1)]
        bracket_steps_dist=[np.sqrt(bracket_steps[i][0]**2+bracket_steps[i][1]**2) for i in range(len(bracket_steps))]
        newton_steps_dist=[np.sqrt(newton_steps[i][0]**2+newton_steps[i][1]**2) for i in range(len(newton_steps))]
        steps=np.arange(1,len(bracket_steps)+1,1)
        ratio_steps=[] 
        for i in range(len(bracket_steps_dist)):
            ratio_steps.append(newton_steps_dist[i]/bracket_steps_dist[i])
        


        #Run newton for timing as well
        
        #Time 
        for i in range(20):
            t1=time.time() #start timer 1
            tp1=time.process_time() #start timer 2

            #newtons 
            hw2.newton(xg)
            
            t2 = time.time() #t2-t1 gives wallclock time
            tp2 = time.process_time() #tp2-tp1 gives cpu time -- depends on number of cores!

            newton_wallclock=t2-t1
            newton_cpu_time=tp2-tp1

            average_newton_wallclock+=newton_wallclock
            average_newton_cpu_time+=newton_cpu_time


        #divide by 20 to create the averages
        average_newton_cpu_time=average_newton_cpu_time/20
        average_newton_wallclock=average_newton_wallclock/20
        average_bd_wallclock=average_bd_wallclock/20
        average_bd_cpu_time=average_bd_cpu_time/20



        plt.figure()

        # plt.subplot(121)
        plt.title("Step (h) Comparison for Newtons and Bracket Descent with xg=%s"%xg)
        plt.suptitle('Lawrence Stewart - Created Using bracket_descent_test().')
        plt.xlabel('Iteration number')
        plt.ylabel('Size of step h')
        plt.plot(steps,bracket_steps_dist,label="Bracket Descent",alpha=0.8,color='r')
        plt.plot(steps,newton_steps_dist,label="Newtons Method",alpha=0.7)
        ax = plt.gca()
        plt.grid('on')
        ax.set_facecolor('#D9E6E8')
        plt.legend()
       
        
        # plt.subplot(122)
        # plt.plot(steps,ratio_steps,alpha=0.7,color='r')
        # plt.title("Ratio of Newton Step taken over Bracket Step taken at Each Iteration")
        # plt.xlabel("Iteration")
        # plt.ylabel("Newton Step/ Bracket Step")
        # plt.grid('on')
        # ax = plt.gca()
        # ax.set_facecolor('#D9E6E8')
        # plt.show()


        #Plot timings:
        fig, ax = plt.subplots()
        plt.suptitle('Lawrence Stewart - Created Using bracket_descent_test().')
        index = np.arange(2)
        bar_width = 0.35
        opacity = 0.7
        rects1 = plt.bar(index, (average_newton_wallclock,average_newton_cpu_time), bar_width,alpha=opacity,color='c',label='Newton')
        rects2 = plt.bar(index + bar_width, (average_bd_wallclock,average_bd_cpu_time), bar_width,alpha=opacity,color='m',label='BD')
        plt.xticks(index + bar_width, ('Wallclock Time', 'CPU Time'))
        plt.legend()
        ax = plt.gca()
        ax.set_facecolor('#D9E6E8')
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.title("Average Timings for Bracket Descent and Newtons for xg=%s"%xg)
        plt.ylabel("Time (s)")
        plt.show()

  
    return xf,jf


def performance():
    """ Assess performance of B-D and L-BFGS-B methods. Add input/output as
    needed

    Discussion and assessment:

    hw241.png shows the error of the final result for both bracket descent and L-BFGS-B for a varied selection
    of points. For Bracket descent we see a range of errors from 0.0030 to 0.00025, however, for L-BFGS-B we
    see much lower errors (discussed below). hw242.png shows the number of iterations required before convergance
    for both of the methods for a variety of points. We once again see that bracket descent requires a higher number
    of iterations for every point evaluated (roughly double the number of iterations as L-BFGS-B requires.)

    This is because "L-BFGS uses an estimation to the inverse Hessian matrix to steer its search through variable space" (Scipy user manual)
    in a similar way to the newton method used earlier. This allows for it to jump larger steps in the direction 
    of the minimum.

    Interestingly, hw243.png shows that for all the points L-BFGS-B had a larger clock time, approximately
    0.0012 to 0.0035, (compared to a much lower time of bracket descent.) This is due to the estimation of  
    the inverse hessian that must be computed at every step of he L-BFGS-B algorithm, whilst bracket descent 
    relies on a series of simple if statments.

    However, it is worth noting that Bracket descent has a default module tolerance of 1e-6. L-BFGS-B by default
    has a tolerance of 2.2e-9. This will explain part of the reason for the slower nature of L-BFGS-B, and also 
    the much higher precision, along with the reasons above.

    hw244/5.png shows the effect that amplitude has on both algorithms. We see that for a starting value of
    [-1,-10], amplitude seems to have a much larger effect on Bracket descent. This will be due to the fact that
    L-BFGS-B will implement costj once per loop, whilst bracket descent will implement costj 3 time, one for each 
    vertex of the triangle. 

    Describe 3 distinct features that could be (or have been) added to your Python code to improve its suitability for use as scientific software:

    1) An easy fix would be to introduce the option for the user to have access to other cost functions (Mean square error e.c.t)
    2) More optimization algorithms contained in the fortran module would give the python user a
    greater degree of freedom. Also greating a GUI would make the program more accessable to scientists.
    3) Paralellisation, - using the OpenMp techniques described in the lectures we could parallelize some
    of the fortran code.

    (4 extra - Make the tool an online resource. This way users do not need to download f2py, gfortan
    the exact same compiler, python ect ect. They can just run it online (to create this make a CGI.)

    Have a nice day Dr. Ray.





    """
    x1_vals=[-100, -50, -10,-1]
    x2=-3

    BFGS_no_iterations=[] #Number of iterations before convergance for L-BFGS-B
    BFGS_errors=[]  #Error for L-BFGS-B
    Bracket_errors=[] #Error for bracket descent
    Bracket_no_iterations=[] #Number of iterations before conv for bracket descent
    Bracket_clocktimes=[]
    BFGS_clocktimes=[]

    #Generate data:
    for i in range(len(x1_vals)):
        #set the point that we will evaluate
        xg=(x1_vals[i],x2)

        #perform both optimization methods
        t1=time.time() #start timer 1
        tp1=time.process_time() #start timer 2

        #RUN THE ALGORITHM
        data_bracket=hw2.bracket_descent(xg)
        
        t2 = time.time() #t2-t1 gives wallclock time
        tp2 = time.process_time() #tp2-tp1 gives cpu time -- depends on number of cores!
        Bracket_clocktimes.append(t2-t1)
        

        
        t1=time.time() #start timer 1
        tp1=time.process_time() #start timer 2

        #RUN THE ALGORITHM
        data_LBFGSB = minimize(cost.costj, xg, method = 'L-BFGS-B')
        
        t2 = time.time() #t2-t1 gives wallclock time
        tp2 = time.process_time() #tp2-tp1 gives cpu time -- depends on number of cores!
        BFGS_clocktimes.append(t2-t1)


        

        #Add the errors to the lists:
        temp_bracket=[1,1]-data_bracket[0]
        Bracket_errors.append(np.sqrt(temp_bracket[0]**2+temp_bracket[1]**2))
        temp_bfgs=[1,1]-data_LBFGSB.x
        BFGS_errors.append(np.sqrt(temp_bfgs[0]**2+temp_bfgs[1]**2))

        #Add the number of iterations required till convergance for both:
        Bracket_no_iterations.append(len(hw2.jpath))
        BFGS_no_iterations.append(data_LBFGSB.nit)

    #Plot (using pyplot plt.bar online template)

    #BOX PLOT FOR BRACKET ERRORS
    fig, ax = plt.subplots()
    plt.suptitle('Lawrence Stewart - Created Using performance().')
    index = np.arange(4)
    bar_width = 0.35
    opacity = 0.7
    rects1 = plt.bar(index, tuple(Bracket_errors), bar_width,alpha=opacity,color='c',label='Bracket Descent')
    rects2 = plt.bar(index + bar_width, tuple(BFGS_errors), bar_width,alpha=opacity,color='m',label='L-BFGS-B')
    plt.xticks(index + bar_width, ('(-100, -3)', '(-50,-3)', '(-10, -3)', '(-1,-3)'))
    plt.legend()
    ax = plt.gca()
    ax.set_facecolor('#D9E6E8')
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.title("Errors Plots for varied starting points ")
    plt.ylabel("Error")
    plt.show()
   

    #BOX PLOT FOR NO OF ITERATIONS
    fig, ax = plt.subplots()
    plt.suptitle('Lawrence Stewart - Created Using performance().')
    index = np.arange(4)
    bar_width = 0.35
    opacity = 0.7
    rects1 = plt.bar(index, tuple(Bracket_no_iterations), bar_width,alpha=opacity,color='c',label='Bracket Descent')
    rects2 = plt.bar(index + bar_width, tuple(BFGS_no_iterations), bar_width,alpha=opacity,color='m',label='L-BFGS-B')
    plt.xticks(index + bar_width, ('(-100, -3)', '(-50,-3)', '(-10, -3)', '(-1,-3)'))
    plt.legend()
    ax = plt.gca()
    ax.set_facecolor('#D9E6E8')
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.title("Number of iterations before convergance for varied points ")
    plt.ylabel("Number of iterations")
    plt.show()

    #BOX PLOT FOR CLOCKTIMES
    fig, ax = plt.subplots()
    plt.suptitle('Lawrence Stewart - Created Using performance().')
    index = np.arange(4)
    bar_width = 0.35
    opacity = 0.7
    rects1 = plt.bar(index, tuple(Bracket_clocktimes), bar_width,alpha=opacity,color='c',label='Bracket Descent')
    rects2 = plt.bar(index + bar_width, tuple(BFGS_clocktimes), bar_width,alpha=opacity,color='m',label='L-BFGS-B')
    plt.xticks(index + bar_width, ('(-100, -3)', '(-50,-3)', '(-10, -3)', '(-1,-3)'))
    plt.legend()
    ax = plt.gca()
    ax.set_facecolor('#D9E6E8')
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.title("Clocktimes for a selection of starting points")
    plt.ylabel("Time (s)")
    plt.show()
   


    #---------------ANALYSE NOISE---------------

    #Analyse Effects of Noise at point xg:
    xg=[-1,-10]
 

    #Generate the different noise values - [0,100]
    noise_vals=np.arange(50+1)

    #Allocate the lists for storing the errors and number of iterations before convergance
    # noise_Bracket_errors=[]
    # noise_BFGS_errors=[]
    noise_Bracket_iterations=[]
    noise_BFGS_iterations=[]


    bd_wallclocks=[]
    LBFGSB_wallclocks=[]
    #Run the algorithms for specified noise value
    for i in range(len(noise_vals)):
        #Set noise
        cost.c_noise=True
        cost.c_noise_amp=noise_vals[i]
        #print(cost.c_noise_amp)
        
        #Run algorithms

        #Time for bracket descent:
        t1=time.time() #start timer 1
        tp1=time.process_time() #start timer 2

        #RUN THE ALGORITHM
        data_bracket=hw2.bracket_descent(xg)

        t2 = time.time() #t2-t1 gives wallclock time
        tp2 = time.process_time() #tp2-tp1 gives cpu time -- depends on number of cores!
        bd_wallclock=t2-t1
        bd_cpu_time=tp2-tp1

        bd_wallclocks.append(bd_wallclock)


        
        #Time for bracket descent:
        t1=time.time() #start timer 1
        tp1=time.process_time() #start timer 2

        #RUN THE ALGORITHM
        data_LBFGSB = minimize(cost.costj, xg, method = 'L-BFGS-B')
       
        t2 = time.time() #t2-t1 gives wallclock time
        tp2 = time.process_time() #tp2-tp1 gives cpu time -- depends on number of cores!
        LBFGSB_wallclock=t2-t1
        LBFGSB_cpu_time=tp2-tp1
        LBFGSB_wallclocks.append(LBFGSB_wallclock)
        # #Append the error for each method for the selected noise value
        # noise_temp_bracket=[1,1]-data_bracket[0]
        # noise_Bracket_errors.append(np.sqrt(noise_temp_bracket[0]**2+noise_temp_bracket[1]**2))
        # noise_temp_bfgs=[1,1]-data_LBFGSB.x
        # noise_BFGS_errors.append(np.sqrt(noise_temp_bfgs[0]**2+noise_temp_bfgs[1]**2))

        #Add the number of iterations required to converge to the lists
        noise_Bracket_iterations.append(len(hw2.jpath))
        noise_BFGS_iterations.append(data_LBFGSB.nit)

    #Reset the noise back to default
    cost.c_noise_amp=0.0




    plt.figure(figsize=(14,6))
    plt.suptitle('Lawrence Stewart - Created Using performance().')
    plt.subplot(121)
    plt.plot(noise_vals,LBFGSB_wallclocks,label="L-BFGS-B",alpha=0.7,color='r')
    plt.plot(noise_vals,bd_wallclocks,label="Bracket Descent",alpha=0.7,color='m')
    ax = plt.gca()
    ax.set_facecolor('#D9E6E8')
    plt.xlabel("Noise Amplitude")
    plt.title("Wallclock times at varied noise amplitudes with xg= %s"%(xg))
    plt.ylabel("Time (s)")
    plt.grid("on")
    plt.legend()
    # plt.show()

    plt.subplot(122)
    plt.plot(noise_vals,noise_Bracket_iterations,label = 'L-BFGS-B',alpha=0.7,color='r')
    plt.plot(noise_vals,noise_BFGS_iterations,label = 'Bracket Descent',alpha=0.7,color='m')
    plt.title("Number of Iterations As Amplitude Changes, xg=%s"%xg)
    plt.ylabel("Number of Iterations")
    plt.xlabel("Amplitude")
    ax = plt.gca()
    ax.set_facecolor('#D9E6E8')
    plt.grid("on")
    plt.legend()
    plt.show()


    #OLD ONE
    # # plt.figure(figsize=(14,6))
    # plt.suptitle('Lawrence Stewart - Created Using nejhgughst().')
    # plt.subplot(121)
    # plt.plot(noise_vals,noise_BFGS_errors,alpha = 0.7, label = 'L-BFGS-B')
    # plt.plot(noise_vals,noise_Bracket_errors, alpha = 0.7, label = 'Bracket Descent')
    # plt.title("Error of Minimum for varied Amplitude of noise , xGuess = (-1,-10)")
    # plt.ylabel("Error")
    # plt.xlabel("Amplitude")
    # ax = plt.gca()
    # ax.set_facecolor('#D9E6E8')
    # plt.legend()

   
    
 

    return None





if __name__ == '__main__':
    #Add code here to call newton_test, bracket_descent_test, performance

   # 1) Generate the contour plot of the cost
    print("------ Q1) Running Visualise to generate Contour Plots -----",end='\r')
    visualize(10,10)


    # # # 2) ----------------------------------------Newtons ----------------------------------------------------
    print("------ Q2) Newtons Tests for a few initial guesses -----",end='\r')
    xguess1=[-10,10]
    xguess2=[300,100]
    xguess3=[1000,-1250]



    #Time each test:
   
    #TEST1------
    t1=time.time() #start timer 1
    tp1=time.process_time() #start timer 2

    #RUN THE ALGORITHM
    newton_test(xguess1,display=False)
    
    t2 = time.time() #t2-t1 gives wallclock time
    tp2 = time.process_time() #tp2-tp1 gives cpu time -- depends on number of cores!
    xg1_wc=t2-t1
    xg1_cpu_time=tp2-tp1

    #Generate current distance from the minimum
    xpath1=hw2.xpath
    distances1=[]
    for i in range(len(xpath1)):
        temp=[1,1]-xpath1[i]
        distances1.append(np.sqrt(temp[0]**2+temp[1]**2))
    jpath1=np.copy(hw2.jpath)



    #TEST2------
    t1=time.time() #start timer 1
    tp1=time.process_time() #start timer 2

    #RUN THE ALGORITHM
    newton_test(xguess2,display=False)
    
    t2 = time.time() #t2-t1 gives wallclock time
    tp2 = time.process_time() #tp2-tp1 gives cpu time -- depends on number of cores!
    xg2_wc=t2-t1
    xg2_cpu_time=tp2-tp1

    #Generate current distance from the minimum
    xpath2=hw2.xpath
    distances2=[]
    for i in range(len(xpath2)):
        temp=[1,1]-xpath2[i]
        distances2.append(np.sqrt(temp[0]**2+temp[1]**2))
    jpath2=np.copy(hw2.jpath)


    #TEST3------
    t1=time.time() #start timer 1
    tp1=time.process_time() #start timer 2

    #RUN THE ALGORITHM
    newton_test(xguess3,display=False)
    
    t2 = time.time() #t2-t1 gives wallclock time
    tp2 = time.process_time() #tp2-tp1 gives cpu time -- depends on number of cores!
    xg3_wc=t2-t1
    xg3_cpu_time=tp2-tp1

    #Generate current distance from the minimum
    xpath3=hw2.xpath
    distances3=[]
    for i in range(len(xpath3)):
        temp=[1,1]-xpath3[i]
        distances3.append(np.sqrt(temp[0]**2+temp[1]**2))
    jpath3=np.copy(hw2.jpath)

    


    #Plotting for Newtons
    plt.figure(figsize=(14, 7))  
    plt.suptitle('Lawrence Stewart - Created Using newton_test().')

    #plot the cost function at each point
    plt.subplot(121)
    plt.plot(jpath1,alpha=0.8,color='r',label="%s"%xguess1)
    plt.plot(jpath2,alpha=0.8,color='g',label="%s"%xguess2)
    plt.plot(jpath3,alpha=0.8,color='m',label="%s"%xguess3)
    plt.xlabel("Iteration")
    plt.ylabel("Cost")
    ax = plt.gca()
    ax.set_facecolor('#D9E6E8')
    plt.legend()
    plt.title("Cost at each Iteration of Newtons Method with varied starting points")
    plt.grid('on')


    plt.subplot(122)
    plt.plot(distances1,alpha=0.8,color='r',label="%s"%xguess1)
    plt.plot(distances2,alpha=0.8,color='g',label="%s"%xguess2)
    plt.plot(distances3,alpha=0.8,color='m',label="%s"%xguess3)
    plt.title("Distance from Minimum at Each Iteration with varied starting points")
    plt.xlabel("Iteration")
    plt.ylabel("Distance")
    ax = plt.gca()
    ax.set_facecolor('#D9E6E8')
    plt.legend()
    plt.grid('on')
    plt.show()



    #-------------------------------------------------------------------------------------------------------------
    # 3) ----------------------------------------Bracket Descent ----------------------------------------------------
    print("------ Q3) Bracket Descent -----",end='\r')
    print(bracket_descent_test([100,10],True))




    #-------------------------------------------------------------------------------------------------------------
    # 4) ----------------------------------------Performance ----------------------------------------------------
    print("------ Q4) Performance -----",end='\r')
    print(performance())

    
    










