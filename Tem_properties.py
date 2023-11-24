#------------------------------------------------------------------ Example code for calculating the TEM properties --------------------------------------------------------------------------------------------------
from __future__ import  division
import numpy as np
import matplotlib.pyplot as plt
import os


class Calculate_TEM_Properties():
            def __init__(self):
                        self.work_folder = os.getcwd()
                        #self.Voltage = 1000*1e3
                        self.h = 6.63*1e-34
                        self.e = 1.6*1e-19
                        self.mo = 9.1*1e-31
                        self.c = 3*1e8
                        self.moc2 = self.mo*(self.c**2)


                        #print "These all values are for 300 KV"
#------------------------------------------------------------ Calculating properties starts here -------------------------------------------------------------------------------------------------------------------------------

            def Calculate_Relativistic_WaveLength(self, Voltage):
                        Voltage = Voltage * 1000
                        self.ev = (Voltage * self.e)
                        nume = self.h * self.c
                        denu = (self.e * Voltage) *  ( ( self.e * Voltage ) + (2*self.mo* (self.c**2)))
                        Relativistic_wave_length = (nume/np.sqrt(denu))
                        #print "Relativistic wave length is :", Relativistic_wave_length
                        return Relativistic_wave_length


            def Calculate_Relativistic_Factor(self):
                        Relativistic_Factor = (1 + (self.ev/self.moc2))
                        #print "Relativistic factor is :", Relativistic_Factor
                        return Relativistic_Factor


            def Calculate_K_Wave_Vector(self):
                        K_Wave_Vector = (2*np.pi)/(self.Calculate_Relativistic_WaveLength())
                        #print "Wave vector is :", K_Wave_Vector
                        return K_Wave_Vector


            def Calculate_Relativistic_Mass(self):
                        Relativistic_Mass = (self.mo*self.Calculate_Relativistic_Factor())
                        #print "Relativistic mass of electron is :", Relativistic_Mass
                        return Relativistic_Mass


            def Calculate_Electron_Velocity(self):
                        nume = np.sqrt( self.ev*(self.ev + 2*self.moc2))
                        denu = (self.ev + self.moc2)

                        Velocity = (self.c * (nume/denu))
                        #print "Relativistic velocity is :", Velocity
                        return  Velocity


            def Calculate_Kinetic_Energy(self):
                        Kinetic_Energy = (0.5*self.mo* (self.Calculate_Electron_Velocity()**2 ))
                        #print "The Kinetic Energy is in Kilo Electron Volt: ", Kinetic_Energy
                        return (Kinetic_Energy/(1.6*1e-16))


            def Calculate_Interaction_Parameter(self):
                        nume = (2*np.pi*self.Calculate_Relativistic_Mass()*self.e*self.Calculate_Relativistic_WaveLength())
                        denu = (self.h**2)
                        Interaction_Parameter = (nume/denu)
                        #print "The interaction parameter is :", Interaction_Parameter
                        return Interaction_Parameter


