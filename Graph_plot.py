import numpy as np 
import matplotlib.pyplot as plt 

def Plot_GPoints_graph(cnm,thickness_range,number_of_points_to_plot,list_containng_values):
        extinctioon_distance_round = []
        for x in range(len(list_containng_values)):
            k = abs(round(list_containng_values[x][-1],2))
            extinctioon_distance_round.append(k)

        unique_extinction_distance_list  = []
        repeated_extinction_distance_list = []
        Optimized_h_list = [ ]
        Optimized_k_list = [ ]
        Optimized_l_list =  [ ]
        Optimized_PSF_list = [ ]
        Optimized_Extinction_list = unique_extinction_distance_list
        for i in range(len(list_containng_values)):
            #   repeating_count = extinctioon_distance_round.count(extinctioon_distance_round[i])
              if extinctioon_distance_round[i] not in unique_extinction_distance_list:
                    unique_extinction_distance_list.append(extinctioon_distance_round[i])
                    Optimized_h_list.append(list_containng_values[i][0])
                    Optimized_k_list.append(list_containng_values[i][1])
                    Optimized_l_list.append(list_containng_values[i][2])
                    Optimized_PSF_list.append(list_containng_values[i][5])
                    # repeated_extinction_distance_list.append(repeating_count)

        #Optimized_Extinction_list.append(list_containng_values[index_list][-1])

        final_2d_list_x = []
        final_2d_list_y = []
        final_parameter_list = [] # 2d list , h , k , l , PSF , Extinction_distance
        for j in range(number_of_points_to_plot):
            # y axis is optimized function 
            # x axis is thickness range
            # optimized_function = [(((Optimized_Extinction_list[j]*1e-9)/ (np.pi*cnm)) *  ((np.sin(((np.pi*x*1e-9)/(Optimized_Extinction_list[j]*1e-9))))**2)*np.abs(Optimized_PSF_list[j]) ) for x in thickness_range]  
            optimized_function = []
            x_axis =[]
            # print(Optimized_Extinction_list[j])
            for z in range(thickness_range):
                 
                 
                 optimized_function.append((((Optimized_Extinction_list[j]*1e-9)/ (np.pi*cnm)) *  ((np.sin(((np.pi*z*1e-9)/(Optimized_Extinction_list[j]*1e-9))))**2)*np.abs(Optimized_PSF_list[j]))) 
                 x_axis.append(z)
                 
            parameter_list = []     
            parameter_list.append(Optimized_h_list[j])
            parameter_list.append(Optimized_k_list[j])
            parameter_list.append(Optimized_l_list[j])
            parameter_list.append(Optimized_PSF_list[j])
            parameter_list.append(Optimized_Extinction_list[j])
            final_parameter_list.append(parameter_list)    
            final_2d_list_y.append(optimized_function)
            final_2d_list_x.append(x_axis)
        # print(final_parameter_list)
        ChartData = {"final_2d_list_x":final_2d_list_x,"final_2d_list_y":final_2d_list_y,"final_parameter_list":final_parameter_list}  
        return ChartData




