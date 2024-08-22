#include "Node.h"


Node::Node(const int& index, const bounded_vector<double, 2>& initialCoordinate)
{
    
    index_ = index;
    initialCoordinate_ = initialCoordinate;
    currentCoordinate_ = initialCoordinate;
    for (size_t i = 0; i < 2; i++)
    {
        currentVelocity_(i) = 0.0; 

        currentMomentum_(i) = 0.0; 

        currentMeshVelocity_(i) = 0.0; 
    }  
    currentDensity_ = 0.;
    currentPressure_ = 0.0;
    currentTemperature_ = 0.;
    currentInternalEnergy_ = 0.;
    constrained_ = false;
    NeummanConstrained_=false;
    MeshConstrained_=false;
    mLump_=0.;
    deltaMesh_(0)=0.;
    deltaMesh_(1)=0.;
}


 int Node::getDOFNode()
 {
     return DOFNode_;
 }

void Node::setDOFNode(int& inode)
{
    DOFNode_= inode;
}

int Node::getIndex()
{
    return index_;
}


 bounded_vector<double, 2> Node::getCurrentCoordinate()
    {
     return currentCoordinate_;
    }

void Node::incrementCurrentCoordinate(const int& direction, const double& value)
{
    currentCoordinate_(direction) += value;
}
 void Node::setCurrentVelocity(const bounded_vector<double, 2>& currentVelocity)
 {
    currentVelocity_ = currentVelocity;
 }

 void Node::setCurrentMeshVelocity(const bounded_vector<double, 2>& currentMeshVelocity)
 {
    currentMeshVelocity_ = currentMeshVelocity;
 }

 void Node::setCurrentTemperature(const double& currentTemperature)
 {
    currentTemperature_ = currentTemperature;
 }

 void Node::setLumpedMass(const double& mlump)
 {
 mLump_=mlump;
 }

void Node::clearLumpedMass()
{

mLump_=0.;

}

double Node::getLumpedMass()
{

return mLump_;
}



void Node::setCurrentMomentum(const bounded_vector<double, 2>& currentMomentum)
 {
    currentMomentum_ = currentMomentum;
 }

 void Node::setCurrentPressure(const double& currentPressure)
 {
    currentPressure_ = currentPressure;
 }

  void Node::setCurrentDensity(const double& currentDensity)
 {
    currentDensity_ = currentDensity;
 }

 
 void Node::incrementCurrentPressure(const double& value)
 {
    currentPressure_ += value;
 }

 void Node::setConstrain(const bool& isConstrained)
 {
     constrained_ = isConstrained;
 }


 void Node::setNeummanConstrain(const bool& isConstrained)
 {
     NeummanConstrained_ = isConstrained;
 }

 void Node::setMeshConstrain(const bool& isConstrained)
 {
     MeshConstrained_ = isConstrained;
 }


 bounded_vector<double, 2> Node::getCurrentVelocity()
 {
    return currentVelocity_;
 }

 bounded_vector<double, 2> Node::getCurrentMeshVelocity()
 {
    return currentMeshVelocity_;
 }

 bounded_vector<double, 2> Node::getDeltaMesh()
 {
    return deltaMesh_;
 }


 bounded_vector<double, 2> Node::getCurrentMomentum()
 {
    return currentMomentum_;
 }

void Node::setCurrentInternalEnergy(const double& currentInternalEnergy)
 {
    currentInternalEnergy_ = currentInternalEnergy;
    
 }
 

double Node::getCurrentInternalEnergy()
 {
    
    return currentInternalEnergy_;
 }

void Node::setDeltaMomentum(const bounded_vector<double,2>& deltarhou)
 {
    deltarhou_ = deltarhou;
    
 }

void Node::setDeltaMesh(const bounded_vector<double,2>& deltaMesh)
 {
    deltaMesh_ = deltaMesh;
    
 }

 void Node::setDeltaDensity(const double& deltarho)
 {
    deltarho_ = deltarho;
    
 }

void Node::setDeltaEnergy(const double& deltarhoe)
 {
    deltarhoe_ = deltarhoe;
    
 }

bounded_vector<double, 2> Node::getDeltaMomentum()
 {
    return deltarhou_;
 }


 double Node::getDeltaDensity()
{
    return deltarho_;
}

 double Node::getDeltaEnergy()
{
    return deltarhoe_;
}

 double Node::getCurrentPressure()
{
    return currentPressure_;
}

double Node::getCurrentDensity()
{
    return currentDensity_;
}

bool Node::getConstrain()
{
    return constrained_;
}

bool Node::getNeummanConstrain()
{
    return NeummanConstrained_;
}

bool Node::getMeshConstrain()
{
    return MeshConstrained_;
}

void Node::setDeltat(const double& deltat)
{ 
   deltat_ = deltat;

}

void Node::updateVariables()
{ 
    currentDensity_ += deltarho_;

    currentInternalEnergy_ += deltarhoe_;

    currentMomentum_ += deltarhou_;

}

void Node::updateMesh()
{ 
    for (size_t i = 0; i < 2; i++)
    {
        
        currentMeshVelocity_(i) = deltaMesh_(i)/deltat_;
        currentCoordinate_(i) = currentCoordinate_(i)+deltaMesh_(i);  
        // atualizar as coordenadas da malha: x=x+deltamesh
    } 

}

void Node:: computeSecondaryVariables(const double& calorv, const double& calorp)
{ 
   for (size_t i = 0; i < 2; i++)
    {
        currentVelocity_(i) = currentMomentum_(i) / currentDensity_; 
    } 
   
   currentTemperature_ = (currentInternalEnergy_/currentDensity_ - (currentVelocity_(0)*currentVelocity_(0)+currentVelocity_(1)*currentVelocity_(1))/2.)/calorv;

   currentPressure_ = (calorp-calorv) * currentDensity_ * currentTemperature_;
}

 double Node::getCurrentTemperature()
 {
    return currentTemperature_;
 }