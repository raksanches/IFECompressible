#pragma once

#include <string>

class Material
{
public:
	Material(const int& index, const double& viscosity, const double& density, const double& thermalCond, const double& specHeatp, const double& specHeatv, const double& TempInf);

	~Material();

	int getIndex();

	double getViscosity();

	double getDensity();

	double getThermalCond();

	double getSpecificHeatp();

	double getSpecificHeatv();

	double getUndisturbedTemperature();


	void setIndex(const int& index);
	
	void setViscosity(const double& viscosity);

	void setDensity(const double& density);

private:
	int index_;
	double viscosity_;
	double density_;
	double thermalCond_; 
	double specHeatp_;
	double specHeatv_; 
	double TempInf_;
};

