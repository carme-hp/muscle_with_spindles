import numpy as np
import variables   

scenario_name = variables.scenario_name


def callback_function_contraction(raw_data):
  t = raw_data[0]["currentTime"]
  if True:

    geometry_data_x = raw_data[0]["data"][0]["components"][0]["values"]
    geometry_data_z = raw_data[0]["data"][0]["components"][2]["values"]
    active_pk2_data_33 = raw_data[0]["data"][4]["components"][2]["values"]
    T_current_traction_x = raw_data[0]["data"][6]["components"][0]["values"]
    T_current_traction_z = raw_data[0]["data"][6]["components"][2]["values"]
    T_material_traction_x = raw_data[0]["data"][7]["components"][0]["values"]
    T_material_traction_z = raw_data[0]["data"][7]["components"][2]["values"]


    number_of_nodes = variables.bs_x * variables.bs_y

    max_displacement_middle = []
    max_displacement_first_quarter = []
    max_displacement_second_quarter = []
    z_end = 0
    active_pk2_data_33_middle = 0
    active_pk2_data_33_first_quarter = 0
    active_pk2_data_33_second_quarter = 0
    max_current_traction_x_center = []
    current_traction_z_end = 0
    current_traction_z_first_quarter = 0
    current_traction_z_second_quarter = 0
    max_material_traction_x_center = []
    material_traction_z_end = 0
    material_current_traction_z_first_quarter = 0
    material_traction_z_second_quarter = 0




    for i in range(number_of_nodes):
      max_displacement_middle.append(geometry_data_x[int(number_of_nodes*(variables.bs_z -1)/2) + i])
      max_displacement_first_quarter.append(geometry_data_x[int(number_of_nodes*(variables.bs_z -1)/4) + i])
      max_displacement_second_quarter.append(geometry_data_x[int(number_of_nodes*(3*variables.bs_z -1)/4) + i])

      z_end += geometry_data_z[number_of_nodes*(variables.bs_z -1) + i]

      active_pk2_data_33_middle += active_pk2_data_33[int(number_of_nodes*(variables.bs_z -1)/2) + i]
      active_pk2_data_33_first_quarter += active_pk2_data_33[int(number_of_nodes*(variables.bs_z -1)/4) + i]  
      active_pk2_data_33_second_quarter += active_pk2_data_33[int(number_of_nodes*(3*variables.bs_z -1)/4) + i]

      max_current_traction_x_center.append(T_current_traction_x[int(number_of_nodes*(variables.bs_z -1)/2) + i])

      current_traction_z_end += T_current_traction_z[number_of_nodes*(variables.bs_z -1) + i]
      current_traction_z_first_quarter += T_current_traction_z[int(number_of_nodes*(variables.bs_z -1)/4) + i]
      current_traction_z_second_quarter += T_current_traction_z[int(number_of_nodes*(3*variables.bs_z -1)/4) + i]

      max_material_traction_x_center.append(T_material_traction_x[int(number_of_nodes*(variables.bs_z -1)/2) + i])

      material_traction_z_end += T_material_traction_z[number_of_nodes*(variables.bs_z -1) + i]
      material_current_traction_z_first_quarter += T_material_traction_z[int(number_of_nodes*(variables.bs_z -1)/4) + i]  
      material_traction_z_second_quarter += T_material_traction_z[int(number_of_nodes*(3*variables.bs_z -1)/4) + i] 


    max_displacement_middle = np.max(max_displacement_middle)
    max_displacement_first_quarter = np.max(max_displacement_first_quarter)
    max_displacement_second_quarter = np.max(max_displacement_second_quarter)

    z_end = z_end/number_of_nodes
    
    active_pk2_data_33_middle = active_pk2_data_33_middle/number_of_nodes
    active_pk2_data_33_first_quarter = active_pk2_data_33_first_quarter/number_of_nodes
    active_pk2_data_33_second_quarter = active_pk2_data_33_second_quarter/number_of_nodes

    max_current_traction_x_center = np.max(max_current_traction_x_center) 
    current_traction_z_end = current_traction_z_end/number_of_nodes 
    current_traction_z_first_quarter = current_traction_z_first_quarter/number_of_nodes 
    current_traction_z_second_quarter = current_traction_z_second_quarter/number_of_nodes

    max_material_traction_x_center = np.max(max_material_traction_x_center)
    material_traction_z_end = material_traction_z_end/number_of_nodes
    material_current_traction_z_first_quarter = material_current_traction_z_first_quarter/number_of_nodes
    material_traction_z_second_quarter = material_traction_z_second_quarter/number_of_nodes


    f = open("out/" + scenario_name + "_output.csv", "a")   # f = open("out/" + scenario_name + "output_" + str(variables.prestretch_force) + "N.csv", "a")
    f.write(str(t))
    f.write(",")
    f.write(str(max_displacement_middle))
    f.write(",")
    f.write(str(max_displacement_first_quarter))
    f.write(",")
    f.write(str(max_displacement_second_quarter))
    f.write(",")
    f.write(str(z_end))
    f.write(",")
    f.write(str(active_pk2_data_33_middle))
    f.write(",")
    f.write(str(active_pk2_data_33_first_quarter))
    f.write(",")
    f.write(str(active_pk2_data_33_second_quarter))
    f.write(",")
    f.write(str(max_current_traction_x_center))
    f.write(",")
    f.write(str(current_traction_z_end))
    f.write(",")
    f.write(str(current_traction_z_first_quarter))
    f.write(",")
    f.write(str(current_traction_z_second_quarter))
    f.write(",")
    f.write(str(max_material_traction_x_center))
    f.write(",")
    f.write(str(material_traction_z_end))
    f.write(",")
    f.write(str(material_current_traction_z_first_quarter))
    f.write(",")
    f.write(str(material_traction_z_second_quarter))
    f.write("\n")
    f.close()
