#pragma once

#include <boost/numeric/ublas/vector.hpp>
#include <vector>

using namespace boost::numeric::ublas;

class Node
{

public:
Node(const int& index, const bounded_vector<double, 2>& initialCoordinate);


 int getIndex(); 

 double getDeltatxTimeStep(); 

 int getDOFNode(); 
 
 void setDOFNode(int& inode);
 
 bounded_vector<double, 2> getCurrentCoordinate(); 


 void setConstrain(const bool& isConstrained); 

 void setNeummanConstrain(const bool& isConstrained); 

void setMeshConstrain(const bool& isConstrained);

 void setCurrentCoordinate(const bounded_vector<double, 2>& currentCoordinate);

 void incrementCurrentCoordinate(const int& direction, const double& value);

 void setCurrentVelocity(const bounded_vector<double, 2>& currentVelocity);

 void setCurrentMeshVelocity(const bounded_vector<double, 2>& currentMeshVelocity);

 void setCurrentMomentum(const bounded_vector<double, 2>& currentMomentum);

 void setCurrentInternalEnergy(const double& currentInternalEnergy);

 void setDeltaMomentum(const bounded_vector<double,2>& deltarhou);

 void setDeltaMesh(const bounded_vector<double,2>& deltaMesh);

 void setDeltaDensity(const double& deltarho);

 void setDeltaEnergy(const double& deltarhoe);
 
 void setCurrentPressure(const double& currentPressure);

 void setCurrentTemperature(const double& currentTemperature);

 void setCurrentDensity(const double& currentDensity);

 void setLumpedMass(const double& mlump);

 void clearLumpedMass();

 double getLumpedMass();

 void incrementCurrentPressure(const double& value);

 bounded_vector<double, 2> getCurrentVelocity();

 bounded_vector<double, 2> getCurrentMomentum();

 bounded_vector<double, 2> getCurrentMeshVelocity();
 bounded_vector<double, 2> getDeltaMesh();

 double getCurrentInternalEnergy();

 bounded_vector<double, 2> getDeltaMomentum();

 double getCurrentPressure();

 double getCurrentTemperature();

 double getDeltaDensity();

 double getDeltaEnergy();

 double getCurrentDensity();

 bool getConstrain();

bool getNeummanConstrain();

bool getMeshConstrain();

void setDeltat(const double& deltat);

 void updateVariables();

 void updateMesh();

 void computeSecondaryVariables(const double& calorv, const double& calorp);


private: 
 int index_; 
 bounded_vector<double, 2> initialCoordinate_; 

 bounded_vector<double, 2> currentCoordinate_; 
 
 bounded_vector<double, 2> currentVelocity_;

 bounded_vector<double, 2> currentMomentum_;

 bounded_vector<double, 2> currentMeshVelocity_;

 bounded_vector<double, 2> deltarhou_;

 bounded_vector<double, 2> deltaMesh_;


 double currentPressure_, currentDensity_, deltarho_, currentInternalEnergy_, currentTemperature_, deltarhoe_, mLump_, deltat_;

 bool constrained_, NeummanConstrained_, MeshConstrained_; 

 int DOFNode_; 

};
