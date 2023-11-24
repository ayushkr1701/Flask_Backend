#----------------------------This is a module file which extract the information from the input file ------------------------------------------------------------------------------------
#----- ----------------------This is done from the Regular expression library of python -------------------------------------------------------------------------------------------------
#--------- ------------------WILL THIS WORK FOR THE FCC ATOMS, WHERE YOU HAVE ONLY ONE ATOMS IN THE THE CENTER 

import re
import os

#********************************** Class Begins from here ********************************************************* --------------------------------------------------------------------------------------------------------------------

class Extract_Structure_Info:
            def __init__(self):
                        self.My_Source_Folder = os.getcwd()
                        

            def Extract_Info(self, Structure_File_Path):
                        # print(self.My_Source_Folder)
                        Structure_File_Name = os.path.basename(Structure_File_Path)
                        File_Extension = (Structure_File_Name.split("."))[1]

            #----------- Directing towards the Specific format. i.e, WIEN2K, VASP, or other structure file format -------------------------------------------------------------
            #----------- CASE 1: --- WIEN2K Structure file --------------------------------------------------------------------------------------------------
                        if (File_Extension =="struct"):
                                    # print(Structure_File_Path)
                                    return self.Extract_From_WIEN2K_Format(Structure_File_Path)


            #--------- CASE 2: -- XYZ Format -----------------------------------------------------------------------------------------------------------------------------
            #--------- You can add other method here to extend this class --------------------------------------------------------------------------------


#*****************************************************************************************************************
#********************* Function which extract the information from the different type of structure file ----------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

            def Extract_From_WIEN2K_Format(self, WIEN2K_Structure_File_Path):

                        with open("%s" % (WIEN2K_Structure_File_Path)) as WIEN2K_Struct_File:

                        #----------- Getting the Material Name ----------------------------------------------------------------------------
                                    Material_Name = (WIEN2K_Struct_File.readline()).split(" ")[0]
                                    wien2k_struct_file_info = WIEN2K_Struct_File.read()
                                    # print(type(wien2k_struct_file_info))
                                    # print(wien2k_struct_file_info)
                                    
                                    #Readlines_data = WIEN2K_Struct_File.readlines()

                        #---------- Getting the Lattice type ------------------------------------------------------------------------------------------
                                    Lattice_Type_Info = (re.findall(r'.*NONEQUIV*.+\s', wien2k_struct_file_info))
                                    Lattice_Type =  Lattice_Type_Info[0].split(" ")[0]
                                    Inequivalent_Atoms = int( [x for x in Lattice_Type_Info[0].split(" ") if x.isdigit()][0])
                                    # print("Inequivalent atoms :",Inequivalent_Atoms)


                        #----------- Making the list of the multiplicity from the file -------------------------------------------------------

                                    Mult_List = re.findall(r'MULT=\s+([0-9]*)', wien2k_struct_file_info)
                                    # print("Mult_List :", Mult_List)
                                    #print Mult_List
                        #---------- Find the atom Name and the atomic number ---------------------------------------------------------

                                    Atom_Name = re.findall(r'.*NPT.*', wien2k_struct_file_info)
                                    #print Atom_Name
                                    Atom_Name_List = []
                                    Atom_Z_List = []
                                    for i in range(len(Atom_Name)):
                                                Atom_Name_Z = Atom_Name[i].split(" ")
                                                number = len(Atom_Name_Z)
                                                Atom_Name_List.append(Atom_Name_Z[0])
                                                Atom_Z_List.append(Atom_Name_Z[number-1])

                                    # print(Atom_Name_List)
                                    # print(Atom_Z_List)
                                    
                        #----------Getting all the coordinates ---------------------------------------------------------------------------------------------------------------------
                        # ** - Make sure that the coordinate are place after x,y,z without and space, i.e, x=0.xxx, y-0.yyy, z=0.zzz
                        #

                                    X_Coordinate_List = []
                                    Y_Coordinate_List = []
                                    Z_Coordinate_List = []

                                    #X_Coordinate = re.findall(r'X=.*', wien2k_struct_file_info)
                                    X_Coordinate = re.findall(r'X=\w.[0-9]*', wien2k_struct_file_info)
                                    Y_Coordinate = re.findall(r'Y=\w.[0-9]*', wien2k_struct_file_info)
                                    Z_Coordinate = re.findall(r'Z=\w.[0-9]*', wien2k_struct_file_info)

                                    #print X_Coordinate
                                    for i in range(len(X_Coordinate)):
                                                X_Coordinate_List.append( ((X_Coordinate[i].split("="))[1] ))
                                                Y_Coordinate_List.append( Y_Coordinate[i].split("=")[1])
                                                Z_Coordinate_List.append( Z_Coordinate[i].split("=")[1])
                                    
                                   
                        WIEN2K_Struct_File.close()

                        #-------- finding the lattice paramters and angles ----------------------------------------------
                        #-------- A list contains lattice parameeter in bohr unit and lattice angles --------------

                        with open("%s" % (WIEN2K_Structure_File_Path)) as WIEN2K_Struct_File:
                                    Readlines_data = WIEN2K_Struct_File.readlines()
                                    Lattice_parameter_match = [line for line in Readlines_data if "MODE OF CAL" in line]
                                    index = Readlines_data.index(Lattice_parameter_match[0])
                                    lattice_index = (index+1)       #----- line number for the lattice paramter lines ---------------------
                                    Lattice_Parameter_Angle = list(filter( None, (((Readlines_data)[lattice_index]).split(" ")) ))

                        Lattice_Parameter_Angle_List = []
                        #---- Converting bohr to angstron and then to the nano-meter distance -------------------------------------------
                        Lattice_Parameter_Angle_List.append(((float(Lattice_Parameter_Angle[0]) * 0.529177) / 10))
                        Lattice_Parameter_Angle_List.append(((float(Lattice_Parameter_Angle[1]) * 0.529177) / 10))
                        Lattice_Parameter_Angle_List.append(((float(Lattice_Parameter_Angle[2]) * 0.529177) / 10))
                        Lattice_Parameter_Angle_List.append(float(Lattice_Parameter_Angle[3]))
                        Lattice_Parameter_Angle_List.append(float(Lattice_Parameter_Angle[4]))
                        Lattice_Parameter_Angle_List.append(float(Lattice_Parameter_Angle[5]))

                        #print Lattice_Parameter_Angle_List

                        WIEN2K_Struct_File.close()

                        #---------- Return the output to the G optimization code ----------------------------------------------------------------------------------------
                        # print("Material_Name : ",Material_Name)
                        # print("\nLattice_Type : ",Lattice_Type)
                        # print("\nInequivalent Atoms : ",Inequivalent_Atoms)
                        # print("\nMult_List : ",Mult_List)
                        # print("\nAtom_Name_List : ",Atom_Name_List)
                        # print("Atom_Z_List : ",Atom_Z_List)
                        # print("\nLattice_Parameter : ", Lattice_Parameter_Angle_List)
                        # print("\nX_Coordinate_List : ",X_Coordinate_List)
                        # print("\nY_Coordinate_List : ",Y_Coordinate_List)
                        # print("\nZ_Coordinate_List : ",Z_Coordinate_List)

                        
                        return (Material_Name, Lattice_Type, Inequivalent_Atoms, Lattice_Parameter_Angle_List, Mult_List, Atom_Name_List, Atom_Z_List, X_Coordinate_List, Y_Coordinate_List,Z_Coordinate_List)



                        #------------- Close the opened WIEN2K file -----------------------------------------------------------------------------------------------------------


# class_Extract_Structure_Info = Extract_Structure_Info()
# class_Extract_Structure_Info.Extract_Info('./FeGe.struct')
# class_Extract_Structure_Info.Extract_Info('./BaFe12O19.struct')
