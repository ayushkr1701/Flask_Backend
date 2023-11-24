# import json
import re


class Extract_Structure_Info:
    def __init__(self, string_a):
        self.string_a = string_a

    def Extract_Info(self):
        line_split_list = self.string_a.split('\n')

        # initializing a count to change lines while reading the struct file
        count = 0

        # getting material name
        material_name = line_split_list[0]
        material_name = material_name.rstrip()
        if '\r' in material_name:
            material_name = material_name.replace('\r', '')
        count += 1

        # getting lattice type
        lattice_type = line_split_list[count][0]
        if '\r' in lattice_type:
            lattice_type = lattice_type.replace('\r', '')

        # finding number of non eq atoms
        find_index_noneq = line_split_list[count].rfind(':')
        non_eq_atoms = str(line_split_list[count][find_index_noneq+1:])
        non_eq_atoms = non_eq_atoms.strip()
        non_eq_atoms = float(non_eq_atoms)

        non_eq_atoms = str(non_eq_atoms)
        if '\r' in non_eq_atoms:
            non_eq_atoms = non_eq_atoms.replace('\r', '')
        non_eq_atoms = float(non_eq_atoms)
        count += 1

        # getting calc and unit
        match_calc = re.search("CALC", line_split_list[count])
        match_unit = re.search("unit", line_split_list[count])

        mode_of_calculation = line_split_list[count][match_calc.end(
        )+1:match_unit.start()]
        unit = line_split_list[count][match_unit.end()+1:]
        count += 1

        # finding lattice parameters
        lattice_parameters = line_split_list[count].split(" ")
        lattice_parameters_final = []
        counting_variable = 0
        for ele in lattice_parameters:  # removing empty spaces
            if ele.strip():
                if '\r' in ele:
                    ele = ele.replace('\r', '')
                #
                if counting_variable < 3:
                    lattice_parameters_final.append((float(ele)*0.529177)/10)
                else:
                    lattice_parameters_final.append(float(ele))
                counting_variable += 1

        count += 1

        # defining empty list
        atom_name_list = []
        atom_z_list = []
        x_coordinate_list = []
        y_coordinate_list = []
        z_coordinate_list = []
        isplit_list = []
        mult_list = []
        npt_list = []
        r_list = []
        rmt_list = []

        # running a nested for loop to find the values of the above parameters
        for i in range(int(non_eq_atoms)):
            x_start = line_split_list[count].rfind('X') + 2
            y_start = line_split_list[count].rfind('Y')
            z_start = line_split_list[count].rfind('Z')
            x_cr = line_split_list[count][x_start:y_start]
            y_cr = line_split_list[count][y_start+2:z_start]
            z_cr = line_split_list[count][z_start+2:]

            if '\r' in x_cr:
                x_cr = x_cr.replace('\r', '')
            x_coordinate_list.append(float(x_cr))
            if '\r' in y_cr:
                y_cr = y_cr.replace('\r', '')
            y_coordinate_list.append(float(y_cr))
            if '\r' in z_cr:
                z_cr = z_cr.replace('\r', '')

            z_coordinate_list.append(float(z_cr))

            count += 1

            match_mult = re.search("MULT", line_split_list[count])
            match_isplit = re.search("ISPLIT", line_split_list[count])
            mult = line_split_list[count][match_mult.end(
            )+2:match_isplit.start()]
            mult = int(mult)
            mult_list.append(mult)
            isplit = int(line_split_list[count][match_isplit.end()+2:])
            isplit_list.append(isplit)
            for j in range(mult-1):
                count += 1
                x_start = line_split_list[count].rfind('X') + 2
                y_start = line_split_list[count].rfind('Y')
                z_start = line_split_list[count].rfind('Z')
                x_cr = line_split_list[count][x_start:y_start]
                y_cr = line_split_list[count][y_start+2:z_start]
                z_cr = line_split_list[count][z_start+2:]
                if '\r' in x_cr:
                    x_cr = x_cr.replace('\r', '')
                x_coordinate_list.append(float(x_cr))
                if '\r' in y_cr:
                    y_cr = y_cr.replace('\r', '')
                y_coordinate_list.append(float(y_cr))
                if '\r' in z_cr:
                    z_cr = z_cr.replace('\r', '')

                z_coordinate_list.append(float(z_cr))

            count += 1
            match_npt = re.search("NPT= ", line_split_list[count])
            match_ro = re.search("R", line_split_list[count])
            match_rmt = re.search("RMT= ", line_split_list[count])
            match_z = re.search("Z: ", line_split_list[count])

            atom_name = line_split_list[count][:match_npt.start()]
            npt = line_split_list[count][match_npt.end():match_ro.start()]
            ro = line_split_list[count][match_ro.end()+2:match_rmt.start()]
            rmt = line_split_list[count][match_rmt.end():match_z.start()]
            z = line_split_list[count][match_z.end():]
            atom_name = atom_name.strip()
            if '\r' in atom_name:
                atom_name = atom_name.replace('\r', '')
            atom_name_list.append(atom_name)
            npt_list.append(npt)
            r_list.append(ro)
            rmt_list.append(rmt)
            if '\r' in z:
                z = z.replace('\r', '')
            atom_z_list.append(float(z))
            count += 4

        structFileInfo = {'Material_Name': material_name, 'Lattice_Type': lattice_type, 'Inequivalent_Atoms': non_eq_atoms, 'Mult_list': mult_list, 'Atom_Name_List': atom_name_list,
                          'Atom_Z_List': atom_z_list, 'Lattice_Parameter': lattice_parameters_final, 'X_Coordinate_List': x_coordinate_list, 'Y_Coordinate_List': y_coordinate_list, 'Z_Coordinate_List': z_coordinate_list}

        Structure_info_dict = {'output': structFileInfo}
        # print(Structure_info_dict)
        return Structure_info_dict
