The content of this folder includes data and codes generated as part of the manuscript “Glacially enhanced silicate weathering revealed by Holocene lake records”:

Data (all other data appear in the manuscript were previously reported and can be accessed elsewhere):
	IcelandRiverLakeDissolved.xlsx: river and lake dissolved chemistry (Ge/Si, d30Si, d18O, dD)
	IcelandBiogenicSi.xlsx: chemical and isotopic composition of biogenic silica samples and the associated sampling depth/age, sedimentation rates
	Iceland_WFlux.xlsx: calculated nutrient utilization and catchment-scale weathering flux 
	Iceland_WFlux_closed.xlsx: calculated nutrient utilization and catchment-scale with a closed-system model
	Cole_CompiledWaterData.xlsx: compiled hydrochemical data of Icelandic rivers, see Cole et al., 2022, JGR-Earth Surface

Codes:
	york_curve_fit: linear regression that accounts for uncertainties in both y and x
 
	Mann_Kendall: non-parametric trend detection
 
	CalculateWeatheringFlux.m: calculating weathering flux assuming open system
 
	closedRayleigh.m: calculating weathering flux assuming closed system behavior for the lakes
	MKTest_HAK.m: Resampling and trend detection for HAK
	CorrTest_HVT.m: Resampling and detection of linear correlation between ice area/glacial discharge and chemical weathering fluxes.
	ODE_MC_HAK.m: assess whether the lake behaves as an open system, generate Supplementary Figure 6A
	ODE_MC_HVT.m: assess whether the lake behaves as an open system, generate Supplementary Figure 6B
	lake_model.m: ode function of Equations 5-6 in Supplementary Text S2.
	
