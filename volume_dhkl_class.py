#-------------------------- This class has a function which calculates the volume of the given unit cell and interplaner distance

from __future__ import division
import numpy as np

class Do_Calculate_Crystal_Volume_and_dhkl:

            #------- Calculate the Crystal volume with generalized formula --------------- 

            def Calculate_Crystal_Volume(self, a, b, c, alpha, beta, gama):
                        #---- a,b,c are already passed as the nm distance ----------------------------------
                        crystal_volume = (a*b*c)* np.sqrt(1 - (np.cos(np.deg2rad(alpha))**2) -  (np.cos(np.deg2rad(beta))**2)  - (np.cos(np.deg2rad(gama))**2)   + 2 * (np.cos(np.deg2rad(alpha))) *  (np.cos(np.deg2rad(beta)))  * (np.cos(np.deg2rad(gama))) )
                        
                        return crystal_volume

            #----- Format : Lattice Type, miller index(hkl), lattice paramters( a,b,c,alpha, beta, gama)
            
            def Calculate_Dhkl(self, Lattice_Type, h, k, l,  a, b, c, alpha, beta, gama):
                        #------- a,b,c should be pass in the nm distance
                      
                        sinalpha = np.sin(np.deg2rad(alpha))
                        sinbeta = np.sin(np.deg2rad(beta))
                        singama = np.sin(np.deg2rad(gama))

                        cosalpha = np.cos(np.deg2rad(alpha))
                        cosbeta = np.cos(np.deg2rad(beta))
                        cosgama = np.cos(np.deg2rad(gama))
                        # print(Lattice_Type) 
                        # print(sinalpha,sinbeta,singama,cosalpha,cosbeta,cosgama) 
                        #------------------------- Cubic ( SC, BCC, FCC ) Lattice ------------------------------------------------------------------------------------------------
                        if ( Lattice_Type in ["Cubic", "cubic", "CUBIC" "c", "C", "B", "b", "BCC", "bcc", "FCC", "fcc", "F", "f", "SC", "sc","P"]):
                                    nume = a
                                    denu = np.sqrt(  (h**2)   + (k**2) +(l**2) )
                                    dhkl = (nume/denu)
                                    #print "this is cubic test"
                                    # print(dhkl)
                                    return (1/dhkl)

                        #------------------- Tetragonal crystal structure ------------------------------------------------------------------------------------------------------------
                        if ( Lattice_Type in ["Tetragonal", "TET", "tet", "T","t" ]):
                                    nume = (a * c)
                                    denu = np.sqrt(  ((c**2) *((h**2) +(k**2))) + ((l*a)**2) )
                                    dhkl = (nume/denu)
                                    return (1/dhkl)


                        #----------------- Hexagonal lattice ------------------------------------------------------------------------------------------------------------------------------
                        if ( Lattice_Type in ["hcp", "HCP", "HEXAGONAL", "hexagonal", "Hexagonal", "h", "H" ]):
                                    nume = (np.sqrt(3) * (a * c))
                                    denu = np.sqrt(  ((4*(c**2)) * ( (h**2) + (k**2) + (h*k) ))   +  (3 * ((a*l)**2 ) ) )
                                    dhkl = (nume/denu)
                                    return (1/dhkl)


                        # ----------------- ORTHORHOMBIC lattice ------------------------------------------------------------------------------------------------------------------------------
                        if (Lattice_Type in ["o", "O", "ORTHORHOMBIC", "orthorhombic", "ortho", "Orthorhombic"]):
                                    t = ( ((h**2)/(a**2)) + ((k**2)/(b**2)) + ((l**2)/(c**2)) )
                                    dhkl = (1 / np.sqrt(t))
                                    return (1/dhkl)



                        # ------------------- Rhombohedral lattice ----------------------------------------------------------------------------------------------------------
                        if (Lattice_Type in [ "R", "r", "Rhombhohedral", "RH", "rh", "RHOMBHOHEDRAL"]):

                                    cosalpha = np.cos(np.deg2rad(alpha))
                                    sinalpha = np.sin(np.deg2rad(alpha))

                                    nume = ( (a **2) * ( 1- (3* (cosalpha**2)) + (2*(cosalpha**3)) ) )
                                    denu = ( (((h**2) + (k**2) + (l**2))*(sinalpha**2)) + ((2*( (h*k) + (k*l) + (h*l)))*(cosalpha**2))  - cosalpha )
                                    dhkl = np.sqrt(nume / denu)
                                    return (1/dhkl)


                        #---------------- Monoclinic lattice -------------------------------------------------------------------------------------------------------------------

                        if (Lattice_Type in [ "M", "m", "Monoclinic", "monoclinic",  "MONOCLINIC"]):

                                    cosbeta = np.cos(np.deg2rad(beta))
                                    sinbeta = np.sin(np.deg2rad(beta))

                                    t1 = ((h**2)/(a**2))
                                    t2 = ( ((k**2) * (sinbeta**2)) / (b**2) )
                                    t3 = ((l**2)/(c**2))
                                    t4 = ((-2*h*l*cosbeta)/(a*c))
                                    t5 = (1/sinbeta**2)
                                    tot = (t5 * (t1 + t2 + t3 + t4))
                                    dhkl = (1 / np.sqrt(tot))
                                    return (1/dhkl)

                        #------------ Triclinic lattice ----------------------------------------------------------------------------------------------------------------------------
                        if (Lattice_Type in [ "Tri", "tri", "Triclinic", "triclinic",  "TRICLINIC"]):

                                    nume = ( 1 - (cosalpha**2) - (cosbeta**2) - (cosgama**2) + (2*cosalpha*cosbeta*cosgama))

                                    t1 =  (((h**2)/(a**2))*(sinalpha**2))
                                    t2 =  (((k**2)/(b**2))*(sinbeta**2))
                                    t3 =  (((l**2)/(c**2))*(singama**2))
                                    t4 =  (((2*k*l)/(b*c))*(cosalpha))
                                    t5 = (((2*h*l)/(a*c))*(cosbeta))
                                    t6 = (((2*h*k)/(a*b))*(cosgama))
                                    tot_denu = (t1 + t2 + t3 + t4 + t5 + t6)

                                    dhkl = np.sqrt(nume / tot_denu)

                                    return (1/dhkl)



#--------------------------------- Done -----------------------------------------

