#include "GeometricBoundaryCondition.h"

GeometricBoundaryCondition::GeometricBoundaryCondition(const int& index,
	const std::string& object,
	const bool& restrictedX, 
	const bool& restrictedY,
	const double& componentX,
	const double& componentY,
	const std::string& referenceSystem,
	const std::string& method,
	const double& penaltyParameter)
{
	index_ = index;
	(object[0] == 'p') ? point_ = object : line_ = object;
	componentX_ = componentX;
	componentY_ = componentY;
	restrictedX_ = restrictedX;
	restrictedY_ = restrictedY;
	referenceSystem_ = referenceSystem;
	method_ = method;
	penaltyParameter_ = penaltyParameter;
}

GeometricBoundaryCondition::~GeometricBoundaryCondition() {}

int GeometricBoundaryCondition::getIndex()
{
	return index_;
}

std::string GeometricBoundaryCondition::getPointName()
{
	return point_;
}

std::string GeometricBoundaryCondition::getLineName()
{
	return line_;
}

std::string GeometricBoundaryCondition::getReferenceSystem()
{
	return referenceSystem_;
}

bool GeometricBoundaryCondition::getRestrictedX()
{
	return restrictedX_;
}

bool GeometricBoundaryCondition::getRestrictedY()
{
	return restrictedY_;
}

double GeometricBoundaryCondition::getComponentX()
{
	return componentX_;
}

double GeometricBoundaryCondition::getComponentY()
{
	return componentY_;
}

std::string GeometricBoundaryCondition::getMethod()
{
	return method_;
}

double GeometricBoundaryCondition::getPenaltyParameter()
{
	return penaltyParameter_;
}

void GeometricBoundaryCondition::setIndex(const int& index)
{
	index_ = index;
}

void GeometricBoundaryCondition::setPointName(const std::string& name)
{
	point_ = name;
}

void GeometricBoundaryCondition::setLineName(const std::string& name)
{
	line_ = name;
}

void GeometricBoundaryCondition::setReferenceSystem(const std::string& referenceSystem)
{
	referenceSystem_ = referenceSystem;
}

void GeometricBoundaryCondition::setComponentX(double componentX)
{
	componentX_ = componentX;
}

void GeometricBoundaryCondition::setComponentY(double componentY)
{
	componentY_ = componentY;
}

void GeometricBoundaryCondition::setRestrictedX(bool restrictedX)
{
	restrictedX_ = restrictedX;
}

void GeometricBoundaryCondition::setRestrictedY(bool restrictedY)
{
	restrictedY_ = restrictedY;
}

void GeometricBoundaryCondition::setMethod(const std::string& method)
{
	method_ = method;
}

void GeometricBoundaryCondition::setPenaltyParameter(const double& penaltyParameter)
{
	penaltyParameter_ = penaltyParameter;
}
