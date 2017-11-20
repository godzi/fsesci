#include "configuration_loader.h"
#include "library.h"


using json = nlohmann::json;
const double E = std::exp(1.0);
const double kBoltz = 1.38064852e-5;// (*pN um *)



/// Assign Simulation Parameters
//SimulationParameters assign_simulation_parameters_from_json(SimulationParameters simp, json jsonobjsimp) {
//	if (!(jsonobjsimp["expTime"].empty())) {
//		simp.expTime = stod(jsonobjsimp["expTime"].get<std::string>());
//	}
	/*if (!(jsonobjsimp["microsteps"].empty())) {
		simp.microsteps = (int)round(stod(jsonobjsimp["microsteps"].get<std::string>()));
	}
	if (!(jsonobjsimp["nTotal"].empty())) {
		simp.nTotal = (int)round(stod(jsonobjsimp["nTotal"].get<std::string>()));
	}*/
//	return simp;
//}

// Assign configuration from json object
Configuration assign_config_from_json(Configuration conf, json jsonobj) {
	//// Assign Logger Parameters from json
	if (!(jsonobj["Name"].empty())) {
		conf.loggerParameters.name = jsonobj["Name"].get<std::string>();
	}
	if (!(jsonobj["LoggerParameters"]["FilePath"].empty())) {
		conf.loggerParameters.filepath= jsonobj["LoggerParameters"]["FilePath"].get<std::string>();
	}
	
	//// Assign Model Parameters from json
	
	if (!(jsonobj["ModelParameters"]["expTime"].empty())) {
		conf.modelParameters.expTime = stod(jsonobj["ModelParameters"]["expTime"].get<std::string>());
			
	}
	if (!(jsonobj["ModelParameters"]["bindseed"].empty())) {
		conf.modelParameters.bindseed = stoul(jsonobj["ModelParameters"]["bindseed"].get<std::string>());

	}
	
	if (!(jsonobj["ModelParameters"]["T"].empty())) {
		conf.modelParameters.T= stod(jsonobj["ModelParameters"]["T"].get<std::string>());
		conf.modelParameters.kT = kBoltz*conf.modelParameters.T;
	}
	if (!(jsonobj["ModelParameters"]["Gdepth"].empty())) {
		conf.modelParameters.G = (conf.modelParameters.kT)*stod(jsonobj["ModelParameters"]["Gdepth"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["Gtotal"].empty())) {
		conf.modelParameters.Gtot = (conf.modelParameters.kT)*stod(jsonobj["ModelParameters"]["Gtotal"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["L"].empty())) {
		conf.modelParameters.L = stod(jsonobj["ModelParameters"]["L"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["sigma"].empty())) {
		conf.modelParameters.sigma = stod(jsonobj["ModelParameters"]["sigma"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["unbindsigmathreshold"].empty())) {
		conf.modelParameters.unbindsigmathreshold = stod(jsonobj["ModelParameters"]["unbindsigmathreshold"].get<std::string>());
	}
	
	if (!(jsonobj["ModelParameters"]["MAPKon"].empty())) {
		conf.modelParameters.MAPKon = stod(jsonobj["ModelParameters"]["MAPKon"].get<std::string>());
	}
	
	if (!(jsonobj["ModelParameters"]["gammaMol"].empty())) {
		conf.modelParameters.gammaMol = stod(jsonobj["ModelParameters"]["gammaMol"].get<std::string>());
		conf.modelParameters.DMol = conf.modelParameters.kT / stod(jsonobj["ModelParameters"]["gammaMol"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["gammaBeadL"].empty())) {
		conf.modelParameters.gammaBeadL = stod(jsonobj["ModelParameters"]["gammaBeadL"].get<std::string>());
		conf.modelParameters.DBeadL = conf.modelParameters.kT / stod(jsonobj["ModelParameters"]["gammaBeadL"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["gammaBeadR"].empty())) {
		conf.modelParameters.gammaBeadR = stod(jsonobj["ModelParameters"]["gammaBeadR"].get<std::string>());
		conf.modelParameters.DBeadR = conf.modelParameters.kT / stod(jsonobj["ModelParameters"]["gammaBeadR"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["gammaMT"].empty())) {
		conf.modelParameters.gammaMT =stod(jsonobj["ModelParameters"]["gammaMT"].get<std::string>());
		conf.modelParameters.DMT = conf.modelParameters.kT / stod(jsonobj["ModelParameters"]["gammaMT"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["DMol"].empty())) {
		conf.modelParameters.DMol = conf.modelParameters.kT / stod(jsonobj["ModelParameters"]["gammaMol"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["DBeadL"].empty())) {
		conf.modelParameters.DBeadL = conf.modelParameters.kT / stod(jsonobj["ModelParameters"]["gammaBeadL"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["DBeadR"].empty())) {
		conf.modelParameters.DBeadR = conf.modelParameters.kT / stod(jsonobj["ModelParameters"]["gammaBeadR"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["DMT"].empty())) {
		conf.modelParameters.DMT = conf.modelParameters.kT / stod(jsonobj["ModelParameters"]["gammaMT"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["trapstiffL"].empty())) {
		conf.modelParameters.trapstiffL = stod(jsonobj["ModelParameters"]["trapstiffL"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["trapstiffR"].empty())) {
		conf.modelParameters.trapstiffR = stod(jsonobj["ModelParameters"]["trapstiffR"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MTstiffWeakSlopeL"].empty())) {
		conf.modelParameters.MTstiffWeakSlopeL = stod(jsonobj["ModelParameters"]["MTstiffWeakSlopeL"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MTstiffWeakBoundaryL"].empty())) {
		conf.modelParameters.MTstiffWeakBoundaryL = stod(jsonobj["ModelParameters"]["MTstiffWeakBoundaryL"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MTstiffParabolicAL"].empty())) {
		conf.modelParameters.MTstiffParabolicAL = stod(jsonobj["ModelParameters"]["MTstiffParabolicAL"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MTstiffParabolicBL"].empty())) {
		conf.modelParameters.MTstiffParabolicBL = stod(jsonobj["ModelParameters"]["MTstiffParabolicBL"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MTstiffParabolicCL"].empty())) {
		conf.modelParameters.MTstiffParabolicCL = stod(jsonobj["ModelParameters"]["MTstiffParabolicCL"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MTstiffStrongBoundaryL"].empty())) {
		conf.modelParameters.MTstiffStrongBoundaryL = stod(jsonobj["ModelParameters"]["MTstiffStrongBoundaryL"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MTstiffStrongSlopeL"].empty())) {
		conf.modelParameters.MTstiffStrongSlopeL = stod(jsonobj["ModelParameters"]["MTstiffStrongSlopeL"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MTstiffStrongIntersectL"].empty())) {
		conf.modelParameters.MTstiffStrongIntersectL = stod(jsonobj["ModelParameters"]["MTstiffStrongIntersectL"].get<std::string>());
	}

	if (!(jsonobj["ModelParameters"]["MTstiffWeakSlopeR"].empty())) {
		conf.modelParameters.MTstiffWeakSlopeR = stod(jsonobj["ModelParameters"]["MTstiffWeakSlopeR"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MTstiffWeakBoundaryR"].empty())) {
		conf.modelParameters.MTstiffWeakBoundaryR = stod(jsonobj["ModelParameters"]["MTstiffWeakBoundaryR"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MTstiffParabolicAR"].empty())) {
		conf.modelParameters.MTstiffParabolicAR = stod(jsonobj["ModelParameters"]["MTstiffParabolicAR"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MTstiffParabolicBR"].empty())) {
		conf.modelParameters.MTstiffParabolicBR = stod(jsonobj["ModelParameters"]["MTstiffParabolicBR"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MTstiffParabolicCR"].empty())) {
		conf.modelParameters.MTstiffParabolicCR = stod(jsonobj["ModelParameters"]["MTstiffParabolicCR"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MTstiffStrongBoundaryR"].empty())) {
		conf.modelParameters.MTstiffStrongBoundaryR = stod(jsonobj["ModelParameters"]["MTstiffStrongBoundaryR"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MTstiffStrongSlopeR"].empty())) {
		conf.modelParameters.MTstiffStrongSlopeR = stod(jsonobj["ModelParameters"]["MTstiffStrongSlopeR"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MTstiffStrongIntersectR"].empty())) {
		conf.modelParameters.MTstiffStrongIntersectR = stod(jsonobj["ModelParameters"]["MTstiffStrongIntersectR"].get<std::string>());
	}
	
	if (!(jsonobj["ModelParameters"]["MTlength"].empty())) {
		conf.modelParameters.MTlength = stod(jsonobj["ModelParameters"]["MTlength"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["molstiff"].empty())) {
		conf.modelParameters.molstiff = stod(jsonobj["ModelParameters"]["molstiff"].get<std::string>());
	}

	if (!(jsonobj["ModelParameters"]["feedbackFreq"].empty())) {
		conf.modelParameters.feedbackFreq = stod(jsonobj["ModelParameters"]["feedbackFreq"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["DmblMoveAmplitude"].empty())) {
		conf.modelParameters.DmblMoveAmplitude = stod(jsonobj["ModelParameters"]["DmblMoveAmplitude"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["prestretchTotalForce"].empty())) {
		conf.modelParameters.prestretchTotalForce = stod(jsonobj["ModelParameters"]["prestretchTotalForce"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["movementTotalForce"].empty())) {
		conf.modelParameters.movementTotalForce = stod(jsonobj["ModelParameters"]["movementTotalForce"].get<std::string>());
	}

	//// Assign Initial Conditions from json
	if (!(jsonobj["InitialConditions"]["xPed"].empty())) {
		conf.initialConditions.xPed = stod(jsonobj["InitialConditions"]["xPed"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["xMol"].empty())) {
		conf.initialConditions.initialState.xMol = stod(jsonobj["InitialConditions"]["xMol"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["xMT"].empty())) {
		conf.initialConditions.initialState.xMT = stod(jsonobj["InitialConditions"]["xMT"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["xBeadl"].empty())) {
		conf.initialConditions.initialState.xBeadl = stod(jsonobj["InitialConditions"]["xBeadl"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["xBeadr"].empty())) {
		conf.initialConditions.initialState.xBeadr = stod(jsonobj["InitialConditions"]["xBeadr"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["xTrapl"].empty())) {
		conf.initialConditions.initialState.xTrapl = stod(jsonobj["InitialConditions"]["xTrapl"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["xTrapr"].empty())) {
		conf.initialConditions.initialState.xTrapr = stod(jsonobj["InitialConditions"]["xTrapr"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["direction"].empty())) {
		conf.currentState.direction = stod(jsonobj["InitialConditions"]["direction"].get<std::string>());
		conf.initialConditions.initialState.direction = stod(jsonobj["InitialConditions"]["direction"].get<std::string>());
	}
	//// Assign Dynamic Coordinates from json initial conditions
	if (!(jsonobj["InitialConditions"]["xTrapl"].empty())) {
		conf.currentState.xTrapl = stod(jsonobj["InitialConditions"]["xTrapl"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["xTrapr"].empty())) {
		conf.currentState.xTrapr = stod(jsonobj["InitialConditions"]["xTrapr"].get<std::string>());
	}

	if (!(jsonobj["InitialConditions"]["xMol"].empty())) {
		conf.currentState.xMol = stod(jsonobj["InitialConditions"]["xMol"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["xMT"].empty())) {
		conf.currentState.xMT = stod(jsonobj["InitialConditions"]["xMT"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["xBeadl"].empty())) {
		conf.currentState.xBeadl = stod(jsonobj["InitialConditions"]["xBeadl"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["xBeadr"].empty())) {
		conf.currentState.xBeadr = stod(jsonobj["InitialConditions"]["xBeadr"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["inititalBindState"].empty())) {
		conf.currentState.bindingState = stod(jsonobj["ModelParameters"]["inititalBindState"].get<std::string>());
		conf.initialConditions.initialState.bindingState = stod(jsonobj["ModelParameters"]["inititalBindState"].get<std::string>());
	}
	////
	return conf;
}


////Simulation params loader
//SimulationParameters load_simulationparams(std::string paramInputFilename) {
//	json fulljson = parse_json_string(readfile(paramInputFilename));
//	json jsonsimp = fulljson["SimulationParameters"];
//	SimulationParameters simp = {};
//	return assign_simulation_parameters_from_json(simp, jsonsimp);
	 
//}

//// Configuration creator new
std::vector <Configuration> load_configuration(std::string paramInputFilename, unsigned nThreads) {
	std::vector <std::string> configs = split(paramInputFilename, ";");
	if (configs.size()% nThreads !=0) {
	
		throw std::runtime_error{ "Number of configs is not multiple of number of cores" };
	}
	std::vector <Configuration> configurationsVector;
	for (const auto& config : configs) {
		json fulljson = parse_json_string(readfile(config));
		json defaultjson = fulljson["Configuration"];
		Configuration iterate;
		iterate = assign_config_from_json(iterate, defaultjson);
		configurationsVector.push_back(iterate);
	}

	return configurationsVector;

}


//// Configuration creator
/*
std::vector <Configuration> load_configuration(std::string paramInputFilename) {
	json fulljson = parse_json_string(readfile(paramInputFilename));
	json defaultjson = fulljson["Configuration"];
	Configuration default;
	default = assign_config_from_json(default, defaultjson);



	std::vector <Configuration> configurationsVector;
	Configuration iterate;

	// iterate the array of configurations
	for (json::iterator jsonit = defaultjson["Configurations"].begin(); jsonit != defaultjson["Configurations"].end(); ++jsonit) {
		iterate = default;
		iterate=assign_config_from_json(iterate, *jsonit);
		configurationsVector.push_back(iterate);
	}

	return configurationsVector;

}
*/