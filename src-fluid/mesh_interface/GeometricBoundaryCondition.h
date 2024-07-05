#pragma once

#include <string>
#include <vector>

class GeometricBoundaryCondition
{
public:
	GeometricBoundaryCondition(const int& index, const std::string& object,const bool& restrictedX,const bool& restrictedY, const double& componentX, const double& componentY, const std::string& referenceSystem = "GLOBAL", const std::string& method = "STRONG", const double& penaltyParameter = 1.0e6);

	~GeometricBoundaryCondition();

	int getIndex();

	std::string getPointName();

	std::string getLineName();

	std::string getReferenceSystem();

	bool getRestrictedX();

	bool getRestrictedY();

	double getComponentX();

	double getComponentY();

	std::string getMethod();

	double getPenaltyParameter();

	void setIndex(const int& index);

	void setPointName(const std::string& name);

	void setLineName(const std::string& name);

	void setReferenceSystem(const std::string& referenceSystem);

	void setComponentX(double componentX);

	void setComponentY(double componentY);

	void setRestrictedX(bool restrictedX);

	void setRestrictedY(bool restrictedY);

	void setMethod(const std::string& methhod);

	void setPenaltyParameter(const double& penaltyParameter);

private:
	int index_;
	std::string point_;
	std::string line_;
	std::string referenceSystem_;
	double componentX_;
	double componentY_;
	bool restrictedX_;
	bool restrictedY_;
	std::string method_;
	double penaltyParameter_;

};

