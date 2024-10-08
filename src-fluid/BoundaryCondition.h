#pragma once

#include <vector>
#include <string>

class BoundaryCondition
{
    public:
    BoundaryCondition(const int& index, const int& nodeIndex, bool& restrictedX, bool& restrictedY, const double& componentX, const double& componentY,
                      const std::string& method, const std::string& referenceSystem, const double& penaltyParameter = 1.0e6);
    
    BoundaryCondition(const int& index, const int& elementIndex, const int& elementSide, bool& restrictedX, bool& restrictedY,const double& componentX,
                      const double& componentY, const std::string& method, const std::string& referenceSystem, const double& penaltyParameter = 1.0e6);

    ~BoundaryCondition();

    int getNodeIndex();

    double getComponent(const int& direction);

    bool getRestrictedDir(const int& direction);

    private:
    int index_;
    int nodeIndex_;
    int elementIndex_;
    int elementSide_;
    double componentX_;
	double componentY_;
    bool restrictedX_;
    bool restrictedY_;
    std::string method_;
    double penaltyParameter_;
    std::string referenceSystem_;

};