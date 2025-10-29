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

# import write_to_file to have access to its callback function
sys.path.insert(0, script_path)
import write_to_file 

scenario_name = variables.scenario_name                                                          

variables.meshes.update(
{
  "muscleSpindleMesh": {
    "nElements" :         variables.n_muscle_spindles-1 if n_ranks == 1 else variables.n_muscle_spindles*n_ranks,
    "physicalExtent":     0,
    "physicalOffset":     0,
    "logKey":             "muscle_spindle",
    "inputMeshIsGlobal":  True,
    "nRanks":             n_ranks
  }

})

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
  "description":            "overall solver",    
  "endTime":                variables.end_time,
  "timeStepWidth":          variables.dt_3D,     
  "logTimeStepWidthAsKey":  "dt_neurons",
  "durationLogKey":         "duration_total",
  "timeStepOutputInterval": 100,
  "connectedSlotsTerm1To2": {},         
  "connectedSlotsTerm2To1": {},     # todo: connect input for muscle spindles 
  "Term1": {  # muscle with spindles
    # muscle spindle model solver
    "Heun" : {
      "description":                  "muscle spindle",
      "timeStepWidth":                variables.dt_muscle_spindles,
      "logTimeStepWidthAsKey":        "dt_muscle_spindles",
      "durationLogKey":               "duration_muscle_spindles",
      "initialValues":                [],
      "timeStepOutputInterval":       1e4,
      "inputMeshIsGlobal":            True,
      "dirichletBoundaryConditions":  {},
      "dirichletOutputFilename":      None,                                 # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
      "checkForNanInf":               True,             # check if the solution vector contains nan or +/-inf values, if yes, an error is printed. This is a time-consuming check.
      "nAdditionalFieldVariables":    0,
      "additionalSlotNames":          [],
          
      # cellml model of muscle spindle
      "CellML" : {
        "modelFilename":                          variables.muscle_spindle_cellml_file,           # input C++ source file or cellml XML file
        "initializeStatesToEquilibrium":          False,                                          # if the equilibrium values of the states should be computed before the simulation starts
        "initializeStatesToEquilibriumTimestepWidth": 1e-4,                                       # if initializeStatesToEquilibrium is enable, the timestep width to use to solve the equilibrium equation
        
        # optimization parameters
        "optimizationType":                       "vc",                                           # "vc", "simd", "openmp" type of generated optimizated source file
        "approximateExponentialFunction":         True,                                           # if optimizationType is "vc", whether the exponential function exp(x) should be approximate by (1+x/n)^n with n=1024
        "compilerFlags":                          "-fPIC -O3 -march=native -shared ",             # compiler flags used to compile the optimized model code
        "maximumNumberOfThreads":                 0,                                              # if optimizationType is "openmp", the maximum number of threads to use. Default value 0 means no restriction.
        
        # stimulation callbacks, motor neuron is not stimulated by a callback function, but has a constant stimulation current
        "setSpecificStatesFunction":              None,                                           # callback function that sets states like Vm, activation can be implemented by using this method and directly setting Vm values, or by using setParameters/setSpecificParameters
        #"setSpecificStatesCallInterval":         0,                                              # set_specific_states should be called variables.stimulation_frequency times per ms, the factor 2 is needed because every Heun step includes two calls to rhs
        "setSpecificStatesCallInterval":          0,                                              # 0 means disabled
        "setSpecificStatesCallFrequency":         0,                                              # set_specific_states should be called variables.stimulation_frequency times per ms
        "setSpecificStatesFrequencyJitter":       0,                                              # random value to add or substract to setSpecificStatesCallFrequency every stimulation, this is to add random jitter to the frequency
        "setSpecificStatesRepeatAfterFirstCall":  0.01,                                           # [ms] simulation time span for which the setSpecificStates callback will be called after a call was triggered
        "setSpecificStatesCallEnableBegin":       0,                                              # [ms] first time when to call setSpecificStates
        "additionalArgument":                     None,
                              
        "handleResultFunction":                   None, #handle_result,              # callback function that gets all current values and can do something with them
        "handleResultCallInterval":               1,                          # interval in which handle_result will be called
        "handleResultFunctionAdditionalParameter": None,                      # additional last argument for handle_result
        
        "mappings":                               variables.muscle_spindle_mappings,              # mappings between parameters and algebraics/constants and between outputConnectorSlots and states, algebraics or parameters, they are defined in helper.py
        "parametersInitialValues":                variables.muscle_spindle_parameters_initial_values,  # initial values for the parameters, either once for all instances, or for all instances in array of struct ordering with nParameters_ parameters per instance: [inst0p0, inst0p1, ... inst0pn, inst1p0, inst1p1, ...]
        
        "meshName":                               "muscleSpindleMesh",                                       
        "stimulationLogFilename":                 "out/stimulation.log",

        # output writer for states, algebraics and parameters                
        "OutputWriter" : [
          {"format": "PythonFile", "outputInterval": (int)(2./variables.dt_muscle_spindles*variables.output_timestep_neurons), "filename": "out/muscle_spindles", "binary": True, "fixedFormat": False, "combineFiles": True, "onlyNodalValues": True, "fileNumbering": "incremental"}
        ]
      }
    }
  },
  "Term2": {  # muscle contraction
    "Coupling": {
      "description":            "electromechanics", 
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
                "callback": write_to_file.callback_function_contraction,
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
}
}
  