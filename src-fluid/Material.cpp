#include "Material.h"

Material::Material(const int& index, const double& viscosity, const double& density, const double& thermalCond, const double& specHeatp, const double& specHeatv, const double& TempInf)
{
	index_ = index;
	viscosity_ = viscosity;
	density_ = density;
	thermalCond_ = thermalCond;
	specHeatp_ = specHeatp;
	specHeatv_ = specHeatv;
	TempInf_ = TempInf;
}

Material::~Material() {}

int Material::getIndex()
{
	return index_;
}

double Material::getViscosity()
{
	return viscosity_;
}

double Material::getDensity()
{
	return density_;
}


double Material::getThermalCond()
{
	return thermalCond_;
}

double Material::getSpecificHeatp()
{
return specHeatp_;
}

double Material::getSpecificHeatv()
{
return specHeatv_;
}

double Material::getUndisturbedTemperature()
{
return TempInf_;
}



void Material::setIndex(const int& index)
{
	index_ = index;
}

void Material::setViscosity(const double& viscosity)
{
	viscosity_ = viscosity;
}

void Material::setDensity(const double& density)
{
	density_ = density;
}