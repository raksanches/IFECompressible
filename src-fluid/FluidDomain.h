#pragma once

#include "mesh_interface/Geometry.h"
#include "mesh_interface/Mesh.h"
#include "BoundaryCondition.h"
//#include "triangle.h"


#include <unordered_map>
#include <petscksp.h>
#include <metis.h>
//#include <fstream.h>

class FluidDomain
{
public:

    //Create fluid domain geometry
	FluidDomain(Geometry* geometry, const int &index= 0);


	~FluidDomain();


    //Functions of Fluid Domain

	Node* getNode(const int& index);

	Element* getElement(const int& index);

	Material* getMaterial(const int& index);

	std::vector<BoundaryCondition*> getBoundaryConditions(const std::string& type);

	void setAnalysisParameters(const int& numberOfTimeSteps, const double& deltat, const double& gravity, const double& rhoInf);

	void setInitialVelocities(const bounded_vector<double,2>& uinf);

	void setInitialMeshVelocities(const bounded_vector<double,2>& uMeshIn);

	void setShockCapture(const double& qdif);

	void setSmallVelocityShockCapture(const double& alpha);

	void useSutherland(const bool& sutherland);

	void addNode(const int& index, const bounded_vector<double,2>& initialCoordinate);

	void addSurfaceMaterial(const std::vector<PlaneSurface*> surfaces, const double& viscosity, const double& density, const double& thermalCond, const double& specHeatp, const double& specHeatv, const double& TempInf);

	void addElement(const int& index, const std::vector<int>& nodes, const int& materialIndex, const double& thickness, const std::string& elementType);

	void generateMesh(const std::string& elementType = "T3", const std::string& algorithm = "AUTO", std::string geofile = std::string(), const std::string& gmshPath = std::string(), const bool& plotMesh = true, const bool& showInfo = false);

	void readInput(const std::string& inputFile, const bool& deleteFiles = true);

	std::vector<int> getConstrains();

	std::vector<int> getConstrains1();

	void getNeummanConstrains();

	std::vector<int> getMeshConstrains();

	void addBoundaryConditions();

	void solveCompressibleFlow();

	void solveCompressibleFlowMoving();

	void exportToParaview(const int& timestep);


	void domainDecompositionMETIS();
  
	//fluid variables
  int index_;
	Geometry * geometry_;
	std::vector<Node*> nodes_;
	std::vector<Element*> elements_;
	std::vector<Material*> materials_;
	std::unordered_map<std::string, std::vector<BoundaryCondition*> > boundaryConditions_;
	int numberOfTimeSteps_;
	double deltat_;
	double gravity_;
	double beta_;
	double gamma_;
	idx_t* elementPartition_;
	idx_t* nodePartition_;
	double area_;
	double rhoInf_;
	bounded_vector<double,2> uinf_;
	bounded_vector<double,2> uMeshIn_;
	double qdif_; //shock capture coeficient
	double theta_; //time integrator parameter
	double alpha_;
	bool sutherland_;

	bool analysisType_;
};
