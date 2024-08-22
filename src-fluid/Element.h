#pragma once

#include "Node.h"
#include "Material.h"
#include <tuple>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <vector>

using namespace boost::numeric::ublas;

class Element
{
public:
	Element(const int& index, const std::vector<Node*>& nodes, Material* material, const double& thickness, const std::string& elementType);

	~Element();

	int getIndex();

	Node* getNode(const int& index);

	std::vector<Node*> getNodes();

	std::string getType();

	std::vector<Node*> getSideNodes(const int& side);

	std::pair<vector<double>, vector<double> > gaussLegendre(const int& nga);

	bounded_vector<double,6> getRhsMomentum();

	bounded_vector<double,3> getRhsDensity();

	bounded_vector<double,3> getRhsEnergy();

	matrix<double> hammerQuadrature(const int& nh );

	vector<double> boundaryShapeFunction(const double& xsi);

	vector<double> boundaryDerivativeShapeFunction(const double& xsi);

	vector<double> domainShapeFunction(const double& xsi1, const double& xsi2);

	matrix<double>  domainDerivativeShapeFunction(const double& xsi1, const double& xsi2);

	bounded_vector<double, 2> normalVector(const double& xsi);

	bounded_matrix<double, 2, 2> currentJacobianMatrix(const double& xsi1, const double& xsi2);

	double getArea();

	double getH();

	double getDeltaT();

	bounded_matrix<double, 3, 3> elementContributionsMesh(); 

	std::pair<bounded_vector<double, 6>, bounded_matrix<double, 6, 6> > elementContributionsS1(const double& alphaM, const double& alphaF, const bool& sutherland);

	std::pair<bounded_vector<double, 3>, bounded_matrix<double, 3, 3> > elementContributionsS2(const double& alphaM, const double& alphaF, const bool& sutherland); //alphaM e alphaF

	std::pair<bounded_vector<double, 3>, bounded_matrix<double, 3, 3> > elementContributionsS3(const double& alphaM, const double& alphaF, const bool& sutherland);

	bounded_vector<double,6> elementShockLowSpeedS1();

	bounded_vector<double,3> elementShockLowSpeedS2();

	bounded_vector<double,3> elementShockLowSpeedS3();

	void setAnalysisParameters(const double& deltat, const double& gravity, const double& beta, const double& gamma);

	void setTimeStep(double timeStep);

	double getTimeStep();

//	std::vector<int> setNeummanSides();

private:
	int index_;
	int order_;
	std::vector<Node*> nodes_;
	Material* material_;
	double thickness_;
	std::string elementType_;
	int numberOfIntegrationPoints_;
	matrix<double> domainIntegrationPoints_;
	//std::vector<double> boundaryIntegrationPoints_;
	//std::vector<double> boundaryIntegrationWeights_;
	double deltat_;
	double gravity_;
	double beta_;
	double gamma_;
	double area_;
	double timeStep_;
	double meshStiffness_;
	//std::vector<int> neummanSides_;
	//std::pair<vector<double>, vector<double> > pontosGauss_;
};

