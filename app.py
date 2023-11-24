import json

from flask_cors import CORS
from flask import Flask, request,jsonify

from G_Optimization_class import Do_G_Optimization
from Load_structure_info_new import Extract_Structure_Info
from Graph_plot import Plot_GPoints_graph

app = Flask(__name__)
CORS(app)


@app.route('/')
def App():
    return "Backend server for optimization of G vectors for EMCD experiments"

@app.route("/file_data",methods=['POST'])
def File_Read_Function():
    FileInput = request.get_json() 

    FileInputDictonary = json.loads(FileInput['body'])
    
    global extracted_info_data
    class_Extract_Structure_Info = Extract_Structure_Info(FileInputDictonary['fileData'])
    extracted_info_data = class_Extract_Structure_Info.Extract_Info()

    return jsonify(extracted_info_data)

    
@app.route("/g_optimized_values",methods=['GET','POST'])
def G_Optimization_output_function():
    Magentic_Atoms_response = request.get_json() 
    Magnetic_Atoms = json.loads(Magentic_Atoms_response['magnetic_atom_dict'])
    
    global extracted_data_with_magnetic_atoms
    extracted_data_with_magnetic_atoms = {};
    extracted_data_with_magnetic_atoms['magnetic_atoms'] = Magnetic_Atoms['magneticAtoms']
    extracted_data_with_magnetic_atoms['extracted_data'] = extracted_info_data['output']
    extracted_data_with_magnetic_atoms['other_parameters'] = Magnetic_Atoms['otherPara']


    G_Optimization_class = Do_G_Optimization(extracted_data_with_magnetic_atoms)
    Output_G_points_parameter, Output_Optimized_G_Parameters= G_Optimization_class.ON_DO_CALCULATE_G_OPTIMIZATION()
    
    return {"Output_G_points_parameter": Output_G_points_parameter,"Output_Optimized_G_Parameters":Output_Optimized_G_Parameters} 
 
 
@app.route("/thickness_gpoints_values",methods=['GET','POST'])  
def G_Optimized_Chart_Data():
    Thickness_GPoints  = request.get_json() 
    Thickness_GPoints_Dictionary = json.loads(Thickness_GPoints['thickness_and_gpoints'])
    
    cnm = round(Thickness_GPoints_Dictionary['cnm'],2)
    thickness = int(Thickness_GPoints_Dictionary['thickness'])
    gPoints = int(Thickness_GPoints_Dictionary['gPoints'])
    G_Optimized_data = Thickness_GPoints_Dictionary['G_Optimized_data']
    
    Ouptput_Chart_data = Plot_GPoints_graph(cnm, thickness,gPoints, G_Optimized_data)
    
    return Ouptput_Chart_data
    

if __name__ == '__main__':
    app.run(host='0.0.0.0')