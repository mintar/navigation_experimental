#% /*
#%  * Copyright (c) 2016, David Conner (Christopher Newport University)
#%  * Based on genmprim_unicycleplussideways.m 
#%  * Copyright (c) 2008, Maxim Likhachev
#%  * All rights reserved.
#%  * converted by libermat utility (https://github.com/awesomebytes/libermate)
#%  * 
#%  * 
#%  * Redistribution and use in source and binary forms, with or without
#%  * modification, are permitted provided that the following conditions are met:
#%  * 
#%  *     * Redistributions of source code must retain the above copyright
#%  *       notice, this list of conditions and the following disclaimer.
#%  *     * Redistributions in binary form must reproduce the above copyright
#%  *       notice, this list of conditions and the following disclaimer in the
#%  *       documentation and/or other materials provided with the distribution.
#%  *     * Neither the name of the Carnegie Mellon University nor the names of its
#%  *       contributors may be used to endorse or promote products derived from
#%  *       this software without specific prior written permission.
#%  * 
#%  * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#%  * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#%  * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
#%  * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
#%  * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
#%  * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
#%  * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
#%  * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
#%  * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
#%  * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
#%  * POSSIBILITY OF SUCH DAMAGE.
#%  */

import numpy as np
import matcompat

# if available import pylab (from matlibplot)
# if available import pylab (from matlibplot)
visualize_plt = False
try:
    import matplotlib.pylab as plt
    visualize_plt = True
except ImportError:
    pass

def genmprim_unicycleplussideways(outfilename,visualize_ = False):
    
    visualize = visualize_plt and visualize_  # Plot the primitives

    # Local Variables: basemprimendpts22p5_c, endtheta_c, endx_c, baseendpose_c, additionalactioncostmult, fout, numofsamples, basemprimendpts45_c, primind, basemprimendpts0_c, rv, angle, outfilename, numberofangles, startpt, UNICYCLE_MPRIM_16DEGS, sidestepcostmult, rotation_angle, basemprimendpts_c, forwardandturncostmult, forwardcostmult, turninplacecostmult, endpose_c, backwardcostmult, interpfactor, S, R, tvoverrv, dtheta, intermcells_m, tv, dt, currentangle, numberofprimsperangle, resolution, currentangle_36000int, l, iind, errorxy, interind, endy_c, angleind, endpt
    # Function calls: plot, cos, pi, grid, figure, text, int2str, basemprimendpts33p75_c, pause, axis, sin, pinv, basemprimendpts11p25_c, fprintf, fclose, rem, zeros, genmprim_unicycleplussideways, fopen, round, size
    #%
    #%generates motion primitives and saves them into file
    #%
    #%written by Maxim Likhachev
    #%---------------------------------------------------
    #%
    #%defines
    UNICYCLE_MPRIM_16DEGS = 1
    if UNICYCLE_MPRIM_16DEGS == 1:
        resolution = 0.025
        numberofangles = 16
        #%preferably a power of 2, definitely multiple of 8
        numberofprimsperangle = 9
        #%multipliers (multiplier is used as costmult*cost)
        forwardcostmult = 1.
        backwardcostmult = 5.
        forwardandturncostmult = 3.
        sidestepcostmult = 2.
        turninplacecostmult = 1.
        #%note, what is shown x,y,theta *changes* (that is, dx,dy,dtheta and not absolute numbers)
        #%0 degreees
        basemprimendpts0_c = np.zeros((numberofprimsperangle, 4))
        #%x,y,theta,costmult 
        #%angles are positive counterclockwise
        #%0 theta change
        basemprimendpts0_c[0,:] = np.array(np.hstack((1., 0., 0., forwardcostmult)))
        basemprimendpts0_c[1,:] = np.array(np.hstack((8., 0., 0., forwardcostmult)))
        basemprimendpts0_c[2,:] = np.array(np.hstack((-1., 0., 0., backwardcostmult)))
        #%1/16 theta change
        basemprimendpts0_c[3,:] = np.array(np.hstack((8., 1., 1., forwardandturncostmult)))
        basemprimendpts0_c[4,:] = np.array(np.hstack((8., -1., -1., forwardandturncostmult)))
        #%turn in place
        basemprimendpts0_c[5,:] = np.array(np.hstack((0., 0., 1., turninplacecostmult)))
        basemprimendpts0_c[6,:] = np.array(np.hstack((0., 0., -1., turninplacecostmult)))
        #%sideways maintaining the same heading
        basemprimendpts0_c[7,:] = np.array(np.hstack((0., 1., 0., sidestepcostmult)))
        basemprimendpts0_c[8,:] = np.array(np.hstack((0., -1., 0., sidestepcostmult)))
        #%45 degrees
        basemprimendpts45_c = np.zeros((numberofprimsperangle, 4))
        #%x,y,theta,costmult (multiplier is used as costmult*cost)
        #%angles are positive counterclockwise
        #%0 theta change 
        basemprimendpts45_c[0,:] = np.array(np.hstack((1., 1., 0., forwardcostmult)))
        basemprimendpts45_c[1,:] = np.array(np.hstack((6., 6., 0., forwardcostmult)))
        basemprimendpts45_c[2,:] = np.array(np.hstack((-1., -1., 0., backwardcostmult)))
        #%1/16 theta change
        basemprimendpts45_c[3,:] = np.array(np.hstack((5., 7., 1., forwardandturncostmult)))
        basemprimendpts45_c[4,:] = np.array(np.hstack((7., 5., -1., forwardandturncostmult)))
        #%turn in place
        basemprimendpts45_c[5,:] = np.array(np.hstack((0., 0., 1., turninplacecostmult)))
        basemprimendpts45_c[6,:] = np.array(np.hstack((0., 0., -1., turninplacecostmult)))
        #%sideways maintaining the same heading
        basemprimendpts45_c[7,:] = np.array(np.hstack((-1., 1., 0., sidestepcostmult)))
        basemprimendpts45_c[8,:] = np.array(np.hstack((1., -1., 0., sidestepcostmult)))
        #%22.5 degrees
        basemprimendpts22p5_c = np.zeros((numberofprimsperangle, 4.))
        #%x,y,theta,costmult (multiplier is used as costmult*cost)
        #%angles are positive counterclockwise
        #%0 theta change     
        basemprimendpts22p5_c[0,:] = np.array(np.hstack((2., 1., 0., forwardcostmult)))
        basemprimendpts22p5_c[1,:] = np.array(np.hstack((6., 3., 0., forwardcostmult)))
        basemprimendpts22p5_c[2,:] = np.array(np.hstack((-2., -1., 0., backwardcostmult)))
        #%1/16 theta change
        basemprimendpts22p5_c[3,:] = np.array(np.hstack((5., 4., 1., forwardandturncostmult)))
        basemprimendpts22p5_c[4,:] = np.array(np.hstack((7., 2., -1., forwardandturncostmult)))
        #%turn in place
        basemprimendpts22p5_c[5,:] = np.array(np.hstack((0., 0., 1., turninplacecostmult)))
        basemprimendpts22p5_c[6,:] = np.array(np.hstack((0., 0., -1., turninplacecostmult)))
        #%sideways maintaining the same heading
        basemprimendpts22p5_c[7,:] = np.array(np.hstack((-1., 2., 0., sidestepcostmult)))
        basemprimendpts22p5_c[8,:] = np.array(np.hstack((1., -2., 0., sidestepcostmult)))
    else:
        fprintf(1., 'ERROR: undefined mprims type\n')
        return []
        
    # Common code for unicycle models
    fout = open(outfilename, 'w')
    #%write the header
    fout.write('resolution_m: %f\n'%(resolution))
    fout.write('numberofangles: %d\n'%(numberofangles))
    fout.write('totalnumberofprimitives: %d\n'%(numberofprimsperangle*numberofangles))
    #%iterate over angles
    for angleind in np.arange(1., (numberofangles)+1):
        currentangle = ((angleind-1)*2.*np.pi)/numberofangles
        currentangle_36000int = np.round((angleind-1)*36000./numberofangles)
        if (visualize):
            plt.figure(1)#angleind)  # use angleind to plot each primitive set on different window
            plt.hold(True)
            plt.grid(True)
            plt.axis('equal')
            plt.axis([-10*resolution,10*resolution,-10*resolution,10*resolution])
            plt.title(str(angleind)+" "+str(currentangle_36000int/100.))
        #%iterate over primitives    
        for primind in np.arange(1., (numberofprimsperangle)+1):
            fout.write('primID: %d\n'%(primind-1))
            fout.write('startangle_c: %d\n'%(angleind-1))
            #%current angle
            #%compute which template to use
            if (currentangle_36000int%9000) == 0:
                basemprimendpts_c = basemprimendpts0_c[int(primind)-1,:]
                angle = currentangle
            elif (currentangle_36000int%4500) == 0:
                basemprimendpts_c = basemprimendpts45_c[int(primind)-1,:]
                angle = currentangle-45.*np.pi/180.
                
            elif ((currentangle_36000int-7875)%9000) == 0:
                basemprimendpts_c = 1*basemprimendpts33p75_c[primind, :]  # 1* to force deep copy to avoid reference update below
                basemprimendpts_c[0] = basemprimendpts33p75_c[primind, 2]
                #%reverse x and y
                basemprimendpts_c[1] = basemprimendpts33p75_c[primind, 1]
                basemprimendpts_c[2] = -basemprimendpts33p75_c[primind, 3]
                #%reverse the angle as well
                angle = currentangle-(78.75*np.pi)/180.
                print('78p75\n')
                
            elif ((currentangle_36000int-6750)%9000) == 0:
                basemprimendpts_c = 1*basemprimendpts22p5_c[int(primind)-1,:] # 1* to force deep copy to avoid reference update below
                basemprimendpts_c[0] =  basemprimendpts22p5_c[int(primind)-1,1]
                #%reverse x and y
                basemprimendpts_c[1] =  basemprimendpts22p5_c[int(primind)-1,0]
                basemprimendpts_c[2] = -basemprimendpts22p5_c[int(primind)-1,2]
                #%reverse the angle as well
                #print('%d : %d %d %d onto %d %d %d\n'%(primind-1,basemprimendpts22p5_c[int(primind)-1,0], basemprimendpts22p5_c[int(primind)-1,1], basemprimendpts22p5_c[int(primind)-1,2], basemprimendpts_c[0], basemprimendpts_c[1], basemprimendpts_c[2]))
                angle = currentangle-(67.5*np.pi)/180.
                print('67p5\n')
                
            elif ((currentangle_36000int-5625)%9000) == 0:
                basemprimendpts_c = 1*basemprimendpts11p25_c[primind, :] # 1* to force deep copy to avoid reference update below
                basemprimendpts_c[0] = basemprimendpts11p25_c[primind, 2]
                #%reverse x and y
                basemprimendpts_c[1] = basemprimendpts11p25_c[primind, 1]
                basemprimendpts_c[2] = -basemprimendpts11p25_c[primind, 3]
                #%reverse the angle as well
                angle = currentangle-(56.25*np.pi)/180.
                print('56p25\n')
                
            elif ((currentangle_36000int-3375)%9000) == 0:
                basemprimendpts_c = basemprimendpts33p75_c[primind, :]
                angle = currentangle-(33.75*np.pi)/180.
                print('33p75\n')
                
            elif ((currentangle_36000int-2250)%9000) == 0:
                basemprimendpts_c = basemprimendpts22p5_c[int(primind)-1,:]
                angle = currentangle-(22.5*np.pi)/180.
                print('22p5\n')
                
            elif ((currentangle_36000int-1125)%9000) == 0:
                basemprimendpts_c = basemprimendpts11p25_c[primind, :]
                angle = currentangle-(11.25*np.pi)/180.
                print('11p25\n')
                
            else:
                print('ERROR: invalid angular resolution. angle = %d\n', currentangle_36000int)
                return []
                
            
            #%now figure out what action will be        
            baseendpose_c = basemprimendpts_c[0:3]
            additionalactioncostmult = basemprimendpts_c[3]
            endx_c = np.round((baseendpose_c[0]*np.cos(angle))-(baseendpose_c[1]*np.sin(angle)))
            endy_c = np.round((baseendpose_c[0]*np.sin(angle))+(baseendpose_c[1]*np.cos(angle)))
            endtheta_c = (angleind-1+baseendpose_c[2])%numberofangles
            endpose_c = np.array(np.hstack((endx_c, endy_c, endtheta_c)))
            print "endpose_c=",endpose_c            
            print( 'rotation angle=%f\n'% (angle*180./np.pi))
            #if np.logical_and(baseendpose_c[1] == 0., baseendpose_c[2] == 0.):
                #%fprintf(1, 'endpose=%d %d %d\n', endpose_c(1), endpose_c(2), endpose_c(3));
            
            
            #%generate intermediate poses (remember they are w.r.t 0,0 (and not
            #%centers of the cells)
            numofsamples = 10.
            intermcells_m = np.zeros((numofsamples, 3))
            if UNICYCLE_MPRIM_16DEGS == 1.:
                startpt = np.array(np.hstack((0., 0., currentangle)))
                endpt = np.array(np.hstack(((endpose_c[0]*resolution), 
                                            (endpose_c[1]*resolution), 
                                            (( ((angleind-1+baseendpose_c[2])%numberofangles)*2.*np.pi)/numberofangles))))
                
                print "startpt =",startpt
                print "endpt   =",endpt
                intermcells_m = np.zeros((numofsamples, 3))
                if np.logical_or(np.logical_and(endx_c == 0., endy_c == 0.), baseendpose_c[2] == 0.):
                    #%turn in place or move forward            
                    for iind in np.arange(1., (numofsamples)+1):
                        fraction = float(iind-1)/(numofsamples)                 
                        intermcells_m[int(iind)-1,:] = np.array((startpt[0]+(endpt[0]-startpt[0])*fraction, startpt[1]+(endpt[1]-startpt[1])*fraction, 0))
                        rotation_angle = baseendpose_c[2]*(2.*np.pi/numberofangles)
                        intermcells_m[int(iind)-1,2] = (startpt[2]+rotation_angle*fraction)%(2.*np.pi)
                        #print " ",iind,"  of ",numofsamples," fraction=",fraction," rotation=",rotation_angle
                    
                else:
                    #%unicycle-based move forward or backward  (http://sbpl.net/node/53) 
                    R = np.array(np.vstack((np.hstack((np.cos(startpt[2]),  np.sin(endpt[2])-np.sin(startpt[2]))), 
                                            np.hstack((np.sin(startpt[2]), -np.cos(endpt[2])+np.cos(startpt[2]))) )))

                    S = np.dot(np.linalg.pinv(R), np.array(np.vstack((endpt[0]-startpt[0], endpt[1]-startpt[1])) ))
                    l = S[0]
                    tvoverrv = S[1]
                    rv = (baseendpose_c[2]*2.*np.pi/numberofangles)+ l/tvoverrv
                    tv = tvoverrv*rv

                    #print "R=\n",R
                    #print "Rpi=\n",np.linalg.pinv(R)
                    #print "S=\n",S
                    #print "l=",l
                    #print "tvoverrv=",tvoverrv
                    #print "rv=",rv
                    #print "tv=",tv

                    if l<0.:
                        print('WARNING: l = %d < 0 -> bad action start/end points\n'%(l))
                        l = 0.
                    
                    
                    #%compute rv
                    #%rv = baseendpose_c(3)*2*pi/numberofangles;
                    #%compute tv
                    #%tvx = (endpt(1) - startpt(1))*rv/(sin(endpt(3)) - sin(startpt(3)))
                    #%tvy = -(endpt(2) - startpt(2))*rv/(cos(endpt(3)) - cos(startpt(3)))
                    #%tv = (tvx + tvy)/2.0;              
                    #%generate samples
                    for iind in np.arange(1, numofsamples+1):
                        dt = (iind-1)/(numofsamples)
                        #%dtheta = rv*dt + startpt(3);
                        #%intermcells_m(iind,:) = [startpt(1) + tv/rv*(sin(dtheta) - sin(startpt(3))) ...
                        #%                        startpt(2) - tv/rv*(cos(dtheta) - cos(startpt(3))) ...
                        #%                        dtheta];
                        if (dt*tv)<l:
                            intermcells_m[int(iind)-1,:] = np.array(np.hstack((startpt[0]+dt*tv*np.cos(startpt[2]), 
                                                                               startpt[1]+dt*tv*np.sin(startpt[2]), 
                                                                               startpt[2])))
                        else:
                            dtheta = rv*(dt-l/tv)+startpt[2]
                            intermcells_m[int(iind)-1,:] = np.array(np.hstack((startpt[0] + l*np.cos(startpt[2]) + tvoverrv*(np.sin(dtheta)-np.sin(startpt[2])), 
                                                                               startpt[1] + l*np.sin(startpt[2]) - tvoverrv*(np.cos(dtheta)-np.cos(startpt[2])), 
                                                                               dtheta)))
                            
                        
                    #%correct
                    errorxy = np.array(np.hstack((endpt[0] -intermcells_m[int(numofsamples)-1,0], 
                                                  endpt[1] -intermcells_m[int(numofsamples)-1,1])))
                    #print('l=%f errx=%f erry=%f\n'%(l, errorxy[0], errorxy[1]))
                    interpfactor = np.array(np.hstack((np.arange(0., 1.+(1./(numofsamples)), 1./(numofsamples-1)))))

                    #print "intermcells_m=",intermcells_m
                    #print "interp'=",interpfactor.conj().T

                    intermcells_m[:,0] = intermcells_m[:,0]+errorxy[0]*interpfactor.conj().T
                    intermcells_m[:,1] = intermcells_m[:,1]+errorxy[1]*interpfactor.conj().T
                    
            #%write out
            fout.write('endpose_c: %d %d %d\n'% ( endpose_c[0], endpose_c[1], endpose_c[2]))
            fout.write('additionalactioncostmult: %d\n'%( additionalactioncostmult))
            fout.write('intermediateposes: %d\n'%( matcompat.size(intermcells_m, 1.)))
            for interind in np.arange(1., (matcompat.size(intermcells_m, 1.))+1):
                fout.write('%.4f %.4f %.4f\n'%( intermcells_m[int(interind)-1,0], intermcells_m[int(interind)-1,1], intermcells_m[int(interind)-1,2]))
            
            if (visualize):                
                plt.plot(intermcells_m[:,0], intermcells_m[:,1],linestyle="-",marker="o")
                plt.text(endpt[0], endpt[1], str(endpose_c[2]))
                plt.hold(True)
        #if (visualize):     
        #    plt.waitforbuttonpress()  # uncomment to plot each primitive set one at a time
        
    fout.close()
    if (visualize):
        #    plt.waitforbuttonpress()  # hold until buttom pressed
        plt.show()  # Keep windows open until the program is terminated
    return []