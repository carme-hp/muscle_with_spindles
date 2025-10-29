# This settings file can be used for two different equations:
# - Isotropic hyperelastic material
# - Linear elasticity
#
# arguments: <scenario_name> <force>


import numpy as np
import sys, os
import importlib

# parse arguments
rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])

# add folders to python path
script_path = os.path.dirname(os.path.abspath(__file__))
var_path = os.path.join(script_path, "variables")
sys.path.insert(0, var_path)

import variables

# if first argument contains "*.py", it is a custom variable definition file, load these values
if ".py" in sys.argv[0]:
  variables_path_and_filename = sys.argv[0]
  variables_path,variables_filename = os.path.split(variables_path_and_filename)  # get path and filename 
  sys.path.insert(0, os.path.join(script_path,variables_path))                    # add the directory of the variables file to python path
  variables_module,_ = os.path.splitext(variables_filename)                       # remove the ".py" extension to get the name of the module
  
  if rank_no == 0:
    print("Loading variables from \"{}\".".format(variables_path_and_filename))
    
  custom_variables = importlib.import_module(variables_module, package=variables_filename)    # import variables module
  variables.__dict__.update(custom_variables.__dict__)
  sys.argv = sys.argv[1:]     # remove first argument, which now has already been parsed
else:
  if rank_no == 0:
    print("Warning: There is no variables file, e.g:\n ./fibers ../settings_fibers.py fibers.py\n")
  exit(0)

script_path = os.path.dirname(os.path.abspath(__file__))
var_path = os.path.join(script_path, "variables")
sys.path.insert(0, var_path)


# parameters
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

  

config = {
  "scenarioName":                 scenario_name,                # scenario name to identify the simulation runs in the log file
  "logFormat":                    "csv",                        # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "solverStructureDiagramFile":   "solver_structure.txt",       # output file of a diagram that shows data connection between solvers
  "mappingsBetweenMeshesLogFile": "mappings_between_meshes_log.txt",    # log file for mappings 

  "Meshes": variables.meshes,
  "Solvers": {
    "linearElasticitySolver": {           # solver for linear elasticity
      "relativeTolerance":  1e-10,
      "absoluteTolerance":  1e-10,         # 1e-10 absolute tolerance of the residual    ,
      "maxIterations":      1e4,
      "solverType":         "gmres",
      "preconditionerType": "none",
      "dumpFilename":       "",
      "dumpFormat":         "matlab",
    }, 
    "diffusionSolver": {
      "solverType":                     "cg",
      "preconditionerType":             "none",
      "relativeTolerance":              1e-10,
      "absoluteTolerance":              1e-10,
      "maxIterations":                  1e4,
      "dumpFilename":                   "",
      "dumpFormat":                     "matlab"
    },
    "mechanicsSolver": {
      "solverType":                     "preonly",
      "preconditionerType":             "lu",
      "relativeTolerance":              1e-10,
      "absoluteTolerance":              1e-10,
      "maxIterations":                  1e4,
      "snesLineSearchType":             "l2",
      "snesRelativeTolerance":          1e-5,
      "snesAbsoluteTolerance":          1e-5,
      "snesMaxIterations":              10,
      "snesMaxFunctionEvaluations":     1e8,
      "snesRebuildJacobianFrequency":   5,
      "dumpFilename":                   "",
      "dumpFormat":                     "matlab"
    }
  },

    

  "Coupling": {
    "timeStepWidth":            variables.dt_3D,
    "logTimeStepWidthAsKey":    "dt_3D",
    "durationLogKey":           "duration_3D",
    "endTime":                  variables.end_time,
    "connectedSlotsTerm1To2":   {1:2},  # transfer stress to MuscleContractionSolver gamma
    "connectedSlotsTerm2To1":   None,   # transfer nothing back

    "Term1": { # fibers (FastMonodomainSolver)
      "MultipleInstances": { 
        "ranksAllComputedInstances":    list(range(n_ranks)),
        "nInstances":                   1,

        "instances": [{
          "ranks": [0],

          "StrangSplitting": {
            "timeStepWidth":            variables.dt_splitting,
            "logTimeStepWidthAsKey":    "dt_splitting",
            "durationLogKey":           "duration_splitting",
            "timeStepOutputInterval":   100,
            "connectedSlotsTerm1To2":   None, #{0:0,1:1,2:2,3:3,4:4},
            "connectedSlotsTerm2To1":   None, #{0:0,1:1,2:2,3:3,4:4},

            "Term1": { # reaction term
              "MultipleInstances": {
                "nInstances":   variables.n_fibers,

                "instances": [{
                  "ranks": [0],

                  "Heun": {
                    "timeStepWidth":            variables.dt_0D,
                    "logTimeStepWidthAsKey":    "dt_0D",
                    "durationLogKey":           "duration_0D",
                    "timeStepOutputInterval":   100,

                    "initialValues":                [],
                    "dirichletBoundaryConditions":  {},
                    "dirichletOutputFilename":      None,
                    "inputMeshIsGlobal":            True,
                    "checkForNanInf":               False,
                    "nAdditionalFieldVariables":    0,
                    "additionalSlotNames":          [],
                    "OutputWriter":                 [],

                    "CellML": {
                      "modelFilename":          variables.input_dir + "hodgkin_huxley-razumova.cellml",
                      "meshName":               "fiber{}".format(fiber), 
                      "stimulationLogFilename": "out/" + scenario_name + "stimulation.log",

                      "statesInitialValues":                        [],
                      "initializeStatesToEquilibrium":              False,
                      "initializeStatesToEquilibriumTimeStepWidth": 1e-4,
                      "optimizationType":                           "vc",
                      "approximateExponentialFunction":             True,
                      "compilerFlags":                              "-fPIC -march=native -Wno-deprecated_declarations -shared",
                      "maximumNumberOfThreads":                     0,

                      "setSpecificStatesCallEnableBegin":       variables.specific_states_call_enable_begin,
                      "setSpecificStatesCallFrequency":         variables.specific_states_call_frequency,
                      "setSpecificStatesRepeatAfterFirstCall":  0.01,
                      "setSpecificStatesFrequencyJitter":       [0] ,
                      "setSpecificStatesCallInterval":          0,
                      "setSpecificStatesFunction":              None,
                      "additionalArgument":                     None, 

                      "mappings": {
                        ("parameter", 0):               "membrane/i_Stim",
                        ("parameter", 1):               "Razumova/l_hs",
                        ("parameter", 2):               ("constant", "Razumova/rel_velo"),
                        ("connectorSlot", "vm"):        "membrane/V",
                        ("connectorSlot", "stress"):    "Razumova/activestress",
                        ("connectorSlot", "alpha"):     "Razumova/activation",
                        ("connectorSlot", "lambda"):    "Razumova/l_hs",
                        ("connectorSlot", "ldot"):      "Razumova/rel_velo"
                      },
                      "parametersInitialValues": [0.0, 1.0, 0.0],
                    },
                  }
                } for fiber in range(variables.n_fibers)] 
              }
            },

            "Term2": { # diffusion term
              "MultipleInstances": {
                "nInstances": variables.n_fibers, 

                "OutputWriter": [
                  {
                    "format":             "Paraview",
                    "outputInterval":     int(1.0 / variables.dt_3D * variables.output_interval),
                    "filename":           "out/" + scenario_name + "/fibers",
                    "fileNumbering":      "incremental",
                    "binary":             True,
                    "fixedFormat":        False,
                    "onlyNodalValues":    True,
                    "combineFiles":       True
                  }
                ],

                "instances": [{
                  "ranks": [0],

                  "ImplicitEuler": {
                    "timeStepWidth":            variables.dt_1D,
                    "logTimeStepWidthAsKey":    "dt_1D",
                    "durationLogKey":           "duration_1D",
                    "timeStepOutputInterval":   100,

                    "nAdditionalFieldVariables":    4,
                    "additionalSlotNames":          ["stress", "alpha", "lambda", "ldot"],

                    "solverName":                       "diffusionSolver",
                    "timeStepWidthRelativeTolerance":   1e-10,

                    "dirichletBoundaryConditions":      {},
                    "dirichletOutputFilename":          None,
                    "inputMeshIsGlobal":                True,
                    "checkForNanInf":                   False,
                    "OutputWriter":                     [],

                    "FiniteElementMethod": {
                      "meshName":           "fiber{}".format(fiber),
                      "inputMeshIsGlobal":  True,
                      "solverName":         "diffusionSolver",
                      "prefactor":          variables.diffusion_prefactor,
                      "slotName":           "vm"
                    }
                  }
                } for fiber in range(variables.n_fibers)]
              }
            }
          }
        }]
      },

      "fiberDistributionFile":                              variables.fiber_distribution_file,
      "firingTimesFile":                                    variables.firing_times_file,
      "valueForStimulatedPoint":                            20.0,
      "onlyComputeIfHasBeenStimulated":                     True,
      "disableComputationWhenStatesAreCloseToEquilibrium":  True,
      "neuromuscularJunctionRelativeSize":                  0.0,################################change for no randomness
      "generateGPUSource":                                  True,
      "useSinglePrecision":                                 False
    },

    "Term2": { # solid mechanics (MuscleContractionSolver)
      "MuscleContractionSolver": {
        "Pmax":                         variables.pmax,
        "slotNames":                    ["lambdaContraction", "ldotContraction", "gammaContraction", "TContraction"],
        #"slotNames":                    ["lambda", "ldot", "gamma", "T"],
        "dynamic":                      True,

        "numberTimeSteps":              1,
        "timeStepOutputInterval":       100,
        "lambdaDotScalingFactor":       1,
        "enableForceLengthRelation":    True,
        "mapGeometryToMeshes":          [],

        "OutputWriter": [
          {
            "format":             "Paraview",
            "outputInterval":     int(1.0 / variables.dt_3D * variables.output_interval),
            "filename":           "out/" + scenario_name + "/mechanics",
            "fileNumbering":      "incremental",
            "binary":             True,
            "fixedFormat":        False,
            "onlyNodalValues":    True,
            "combineFiles":       True
          }
        ],

        "DynamicHyperelasticitySolver": {
          "durationLogKey":         "duration_3D",
          "logTimeStepWidthAsKey":  "dt_3D",
          "numberTimeSteps":        1,
          "materialParameters":     variables.material_parameters,
          "density":                variables.rho,
          "timeStepOutputInterval": 1,

          "meshName":                   "mesh3D",
          "fiberDirectionInElement":    variables.fiber_direction,
          "inputMeshIsGlobal":          True,
          "fiberMeshNames":             [],
          "fiberDirection":             [0,0,1],

          "solverName":                 "mechanicsSolver",
          "displacementsScalingFactor":  1.0,
          "useAnalyticJacobian":        True,
          "useNumericJacobian":         False,
          "dumpDenseMatlabVariables":   False,
          "loadFactorGiveUpThreshold":  1,
          "loadFactors":                [],
          "scaleInitialGuess":          False,
          "extrapolateInitialGuess":    True,
          "nNonlinearSolveCalls":       1,

          "dirichletBoundaryConditions":                            variables.contraction_dirichlet_bc, 
          "neumannBoundaryConditions":                              {} if variables.isometric else variables.contraction_neumann_bc, 
          "updateDirichletBoundaryConditionsFunction":              None,
          "updateDirichletBoundaryConditionsFunctionCallInterval":  1,
          "divideNeumannBoundaryConditionValuesByTotalArea":        True,

          "initialValuesDisplacements": [[0, 0, 0] for _ in range(variables.bs_x * variables.bs_y * variables.bs_z)],
          "initialValuesVelocities":    [[0, 0, 0] for _ in range(variables.bs_x * variables.bs_y * variables.bs_z)],
          "constantBodyForce":          (0, 0, 0),

          "dirichletOutputFilename":    "out/" + scenario_name + "/dirichlet_output",
          "residualNormLogFilename":    "out/" + scenario_name + "/residual_norm_log.txt",
          "totalForceLogFilename":      "out/" + scenario_name + "/total_force_log.txt",

          "OutputWriter": [
            {
              "format": "PythonCallback",
              "callback": callback_function_contraction,
              "outputInterval": 1,
            }
          ],
          "pressure":       { "OutputWriter": [] },
          "dynamic":        { "OutputWriter": [] },
          "LoadIncrements": { "OutputWriter": [] }
        }
      }
    }
  }
}
  