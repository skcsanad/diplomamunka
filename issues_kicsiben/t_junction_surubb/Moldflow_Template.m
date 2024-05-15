%extern "C" double STDCALL SolverUserHb3dSpecificVolume(int aDerivsNeeded, double aTemperature, double aPressure, double aCure, double *aDerivs)

%extern "C" double STDCALL SolverUserHb3dSpecificVolumeAtNode(int derivsNeeded, double temperature, double pressure, double cure, double *derivs, size_t nodeId)

%extern "C" double STDCALL SolverUserHb3dViscosity(double temperature, double shearRate, double pressure)

%extern "C" double STDCALL SolverUserHb3dViscosityAtNode(double temperature, double shearRate, double pressure, size_t nodeId)

%extern "C" void STDCALL SolverUserHb3dTimeStepComplete(double time)



%extern "C" void STDCALL SolverUserHb3dAnalysisSetup()

%extern "C" void STDCALL SolverUserHb3dCleanup()

%extern "C" int STDCALL SolverUserWp3dInitialize()

%extern "C" int STDCALL SolverUserWp3dAnalysisSetup()

%extern "C" double STDCALL SolverUserWp3dSpecificVolume(int aNodeLabel, double aPressure, double aTemperature, int aNumShot)

%extern "C" double STDCALL SolverUserWp3dSpecificVolumeTemperature(int aNodeLabel, double aPressure, double aSpecificVolume, int aNumShot)

%extern "C" void STDCALL SolverUserWp3dCleanup()

