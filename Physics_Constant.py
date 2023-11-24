#------------------------------------- Module class containing the Constant values used in the physics -------------------------------------------------------------------------------------------------------
import numpy as np
class Constants:
            def __init__(self):
                        self.Electron_Charge = 1.6*1e-19                                        #----- e :  coulomb : Electron charge -------------
                        self.Plank_Constant = 6.64*1e-34                                        #---  h :  joule.second: Planck constant --------------
                        self.PI = np.pi                                                         #---- pi:  pi value --------------------------
                        self.Velocity_of_Light = 3*1e8                                          #----- c : metre/second :Velocity of light -------------
                        self.Bohr_Radius = 0.529*1e-10                                          #---- ao : meter :Bohr Radius -------------------
                        self.Electron_Rest_Mass = 9.1*1e-31                                     #--- me : kg :  Electron mass -----------------
                        self.Boltzman_constant=1.380649*1e-23                                   #--- K : Joule/Kelvin Boltzman constant ----------
                        self.Bohr_Magneton = 9.274009*1e-24                                     #--- UB : joule/Tesla : Bohr magneton mu --------
                        self.Gravitational_Acceleration = 9.87                                  #----- g: m/s2 --
                        self.h_cross = 1.054571800*1e-34                                        #---- j.second
                        self.Vacuum_Permeability = 4*self.PI*1e-7                               #--- Newton/ampere**2
                        self.Vacuum_Permitivity = 8.854187*1e-12                                #--
                        self.Rydberg_Constant = 1.0973731*1e-7                                  #-- meter^-1
                        self.Stefen_Boltzman_Constant = 5.670367*1e-8
                        self.Avagadro_Number = 6.0221 * 1e23                                    #-- Na : Avagadro number /mol ---------
                        self.Coulomb_Constant_k = 9*1e9                                         #-- 1/4 pi epsilon -- in electrostatic
                        self.Faraday_Constant_F = 96485.309                                     #-- coulomb/mol ------------------
                        self.Proton_Rest_Mass = 1.6726231*1e-27                                 #-- in Kilogram
                        self.Wein_Displacement_w = 2.897756*1e-3                                #- m.Kelvin
                        self.One_Radian_to_Degree = 57.295779
                        self.One_Degreee_to_Radian = 0.0174532

#------------------------------------------ Class ends here ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------