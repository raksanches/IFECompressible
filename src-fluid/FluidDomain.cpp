#include "FluidDomain.h"
std::vector<int> intersection(std::vector<int> vector1, std::vector<int> vector2)
{
	std::sort(vector1.begin(), vector1.end());
	std::sort(vector2.begin(), vector2.end());

	std::vector<int> vector3;
	std::set_intersection(vector1.begin(), vector1.end(), vector2.begin(), vector2.end(), std::back_inserter(vector3));

	return vector3;
}

std::vector<std::string> split(std::string str, std::string delim)
{
	std::istringstream is(str);
	std::vector<std::string> values;
	std::string token;
	while (getline(is, token, ' '))
		values.push_back(token);
	return values;
}

FluidDomain::FluidDomain(Geometry* geometry, const int &index)
{
	geometry_ = geometry;
	index_ = index;
}

FluidDomain::~FluidDomain(){}

//Functions of Fluid Domain
Node* FluidDomain::getNode(const int& index)
{
	return nodes_[index];
}

Element* FluidDomain::getElement(const int& index)
{
	return elements_[index];
}

Material* FluidDomain::getMaterial(const int& index)
{
	return materials_[index];
}

std::vector<BoundaryCondition*> FluidDomain::getBoundaryConditions(const std::string& type)
{
	return boundaryConditions_[type];
}

void FluidDomain::setAnalysisParameters(const int& numberOfTimeSteps, const double& deltat, const double& gravity, const double& theta)
{
	numberOfTimeSteps_ = numberOfTimeSteps;
	deltat_ = deltat;
	gravity_ = gravity;
	theta_ = theta;
}

void FluidDomain::setInitialVelocities(const bounded_vector<double,2>& uinf)
{
	uinf_ = uinf;
}

void FluidDomain::setInitialMeshVelocities(const bounded_vector<double,2>& uMeshIn)
{
	uMeshIn_ = uMeshIn;
}

void FluidDomain::setShockCapture(const double& qdif)
{
	qdif_ = qdif;
}


void FluidDomain::setSmallVelocityShockCapture(const double& alpha)
{
	alpha_ = alpha;
}

void FluidDomain::useSutherland(const bool& sutherland)
{
	sutherland_ = sutherland;
}


void FluidDomain::addNode(const int& index, const bounded_vector<double,2>& initialCoordinate)
{
	Node* n = new Node(index, initialCoordinate);
	nodes_.push_back(n);
}

void FluidDomain::addElement(const int& index, const std::vector<int>& nodesIndex, const int& materialIndex, const double& thickness, const std::string& elementType)
{
	std::vector<Node*> nodes;
	nodes.reserve(nodesIndex.size());
	for (const int& i : nodesIndex)
		nodes.push_back(nodes_[i]);
	Element* e = new Element(index, nodes, materials_[materialIndex], thickness, elementType);
	e->setAnalysisParameters(deltat_, gravity_, beta_, gamma_);
	elements_.push_back(e);
}

void FluidDomain::addSurfaceMaterial(const std::vector<PlaneSurface*> surfaces, const double& viscosity, const double& density, const double& thermalCond, const double& specHeatp, const double& specHeatv, const double& TempInf)
{
	int index = materials_.size();
	Material* m = new Material(index, viscosity, density, thermalCond,  specHeatp, specHeatv, TempInf);
	materials_.push_back(m);
	for (PlaneSurface* surface : surfaces)
		surface->setMaterial(m);
}

void FluidDomain::readInput(const std::string& inputFile, const bool& deleteFiles)
{
	//defyning the maps that are used to store the elements information
	std::unordered_map<int, std::string> gmshElement = { {1, "line"}, {2, "triangle"}, {3, "quadrilateral"}, {8, "line3"}, {9, "triangle6"}, {10, "quadrilateral9"}, {15, "vertex"}, {16, "quadrilateral8"}, {20, "triangle9"}, {21, "triangle10"}, {26, "line4"}, {36, "quadrilateral16"}, {39, "quadrilateral12"} };
	std::unordered_map<std::string, int> numNodes = { {"vertex", 1}, {"line", 2}, {"triangle", 3}, {"quadrilateral", 4}, {"line3", 3}, {"triangle6", 6}, {"quadrilateral8", 8}, {"quadrilateral9", 9}, {"line4", 4}, {"triangle", 9}, {"triangle10", 10}, {"quadrilateral12", 12}, {"quadrilateral16", 16}};
	std::unordered_map<std::string, std::string> supportedElements = { {"triangle", "T3"}, {"triangle6", "T6"}, {"triangle10", "T10"}, {"quadrilateral", "Q4"}, {"quadrilateral8", "Q8"}, {"quadrilateral9", "Q9"}, {"quadrilateral12", "Q12"}, {"quadrilateral16", "Q16"} };
	std::unordered_map<Line*, std::vector< std::vector<int> >> lineElements;

	//opening the .msh file
	std::ifstream file(inputFile);
	std::string line;
	std::getline(file, line); std::getline(file, line); std::getline(file, line); std::getline(file, line);

	//reading physical entities
	int number_physical_entities;
	file >> number_physical_entities;
	std::getline(file, line);
	std::unordered_map<int, std::string> physicalEntities;
	physicalEntities.reserve(number_physical_entities);
	for (int i = 0; i < number_physical_entities; i++)
	{
		std::getline(file, line);
		std::vector<std::string> tokens = split(line, " ");
		int index;
		std::istringstream(tokens[1]) >> index;
		physicalEntities[index] = tokens[2].substr(1, tokens[2].size() - 2);
	}
	std::getline(file, line); std::getline(file, line);

	//reading nodes
	int number_nodes;
	file >> number_nodes;
	nodes_.reserve(number_nodes);
	std::getline(file, line);
	for (int i = 0; i < number_nodes; i++)
	{
		std::getline(file, line);
		std::vector<std::string> tokens = split(line, " ");
		bounded_vector<double,2> coord;
		std::istringstream(tokens[1]) >> coord(0);
		std::istringstream(tokens[2]) >> coord(1);
		addNode(i, coord);
	}
	std::getline(file, line); std::getline(file, line);

	//reading elements
	int number_elements;
	file >> number_elements;
	elements_.reserve(number_elements);
	std::getline(file, line);
	int cont = 0;
	for (int i = 0; i < number_elements; i++)
	{
		std::getline(file, line);
		std::vector<std::string> tokens = split(line, " ");
		std::vector<int> values(tokens.size(), 0);
		for (size_t j = 0; j < tokens.size(); j++)
			std::istringstream(tokens[j]) >> values[j];
		std::string elementType = gmshElement[values[1]];
		int number_nodes_per_element = numNodes[elementType];
		std::vector<int> elementNodes;
		elementNodes.reserve(number_nodes_per_element);
		for (size_t j = 5 ; j < values.size(); j++)
			elementNodes.push_back(values[j]-1);
		std::string name = physicalEntities[values[3]];
		//Adding 2D elements to surfaces
		if (name[0] == 's')
		{
			if (supportedElements.find(elementType) == supportedElements.end())
			{
				std::cout << elementType << " is not supported.\n";
				exit(EXIT_FAILURE);
			}
			PlaneSurface* object = geometry_->getPlaneSurface(name);
			int materialIndex = object->getMaterial()->getIndex();
			double thickness = object->getThickness();
			addElement(cont, elementNodes, materialIndex, thickness, supportedElements[elementType]);
			object->addElementToSurface(elements_[cont]);
			//Verifying if an element side touches a line
			for (auto it = lineElements.begin(); it != lineElements.end();)
			{
				auto& key = it->first;
				auto& value = it->second;
				int it2 = 0;
				for (auto it1 = value.begin(); it1 != value.end();)
				{
					auto nodes = lineElements[key][it2];
					std::sort(nodes.begin(), nodes.end());
					if (intersection(nodes, elementNodes) == nodes)
					{
						key->appendPlaneElement(elements_[cont]);
						it1 = value.erase(it1);
					}
					else
					{
						it1++;
						it2++;
					}
				}
				if (value.size() == 0)
					it = lineElements.erase(it);
				else
					it++;

			}
			
			cont += 1;
		}
		//Adding 1D elements to lines
		else if (name[0] == 'l')
		{
			Line* object = geometry_->getLine(name);
			std::vector<Node*> nodes;
			nodes.reserve(elementNodes.size());
			for (int index : elementNodes)
				nodes.push_back(nodes_[index]);
			object->appendNodes(nodes);
			if (lineElements.count(object) != 1)
				lineElements[object] = std::vector< std::vector<int> >();
			lineElements[object].push_back(elementNodes);
		}
		//Adding a node to point
		else
		{
			Point* object = geometry_->getPoint(name);
			object->addNodeToPoint(getNode(elementNodes[0]));
		}
	}
	//Closing the file
	file.close();
	if (deleteFiles)
		system((remove + inputFile).c_str());
}

std::vector<int> FluidDomain::getConstrains()
{
	std::vector<int> dof;
	//setting dof that are constrained
	int icount=0;
	for (BoundaryCondition* bc : boundaryConditions_["POINT_DISPLACEMENT"])
	{icount++;

		int nodeIndex = bc->getNodeIndex();

				nodes_[nodeIndex]->setConstrain(true);

		for (size_t i = 0; i < 2; i++)
		{

			if (bc->getRestrictedDir(i)==true)
			{
				dof.push_back(2*nodeIndex + i);	
			}


/*			if (bc->getComponent(i).size() != 0)
			{
				std::vector<double> valuei = bc ->getComponent(i);
				if(valuei[0]>0.01)dof.push_back(2*nodeIndex + i);
				
			}
*/

		}
	

	}
	return dof;
}





void FluidDomain::getNeummanConstrains()
{
	//setting dof that are constrained
	for (BoundaryCondition* bc : boundaryConditions_["NEUMANN_POINT"])
	{
		int nodeIndex = bc->getNodeIndex();
		nodes_[nodeIndex]->setNeummanConstrain(true);	
	}

}

std::vector<int> FluidDomain::getConstrains1()
{
	std::vector<int> dof;
	//setting dof that are constrained
	for (BoundaryCondition* bc : boundaryConditions_["POINT_DISPLACEMENT"])
	{
		int nodeIndex = bc->getNodeIndex();
		int dofIndex = nodeIndex;// nodes_[nodeIndex] -> getDOFNode();
		nodes_[nodeIndex]->setConstrain(true);

		/*for (size_t i = 0; i < 2; i++)
		{
			if (bc->getComponent(i).size() != 0)
			{
				dof.push_back(2*dofIndex + i);
			}
		}
*/
	}
	return dof;
}

std::vector<int> FluidDomain::getMeshConstrains()
{
//setting nodes that are fixed

//	for (BoundaryCondition* bc : boundaryConditions_["MESH_DISPLACEMENT"])
//	{
//		int nodeIndex = bc->getNodeIndex();
//		nodes_[nodeIndex]->setMeshConstrain(true);	
//	}
	std::vector<int> dof;
	//setting dof that are constrained
	int icount=0;
	for (BoundaryCondition* bc : boundaryConditions_["MESH_DISPLACEMENT"])
	{icount++;

		int nodeIndex = bc->getNodeIndex();

				nodes_[nodeIndex]->setMeshConstrain(true);

		for (size_t i = 0; i < 2; i++)
		{

			if (bc->getRestrictedDir(i)==true)
			{
				dof.push_back(2*nodeIndex + i);	
			}


//			if (bc->getComponent(i).size() != 0)
//			{
//				std::vector<double> valuei = bc ->getComponent(i);
//				if(valuei[0]>0.01)dof.push_back(2*nodeIndex + i);
//				
//			}


		}
	

	}
	return dof;
}

std::vector<double> FluidDomain::getMeshDisplacementValues()
{
//setting nodes that are fixed

//	for (BoundaryCondition* bc : boundaryConditions_["MESH_DISPLACEMENT"])
//	{
//		int nodeIndex = bc->getNodeIndex();
//		nodes_[nodeIndex]->setMeshConstrain(true);	
//	}
	std::vector<double> dofValue;
	//setting dof that are constrained
	int icount=0;
	for (BoundaryCondition* bc : boundaryConditions_["MESH_DISPLACEMENT"])
	{icount++;

		int nodeIndex = bc->getNodeIndex();

		for (size_t i = 0; i < 2; i++)
		{

			if (bc->getRestrictedDir(i)==true)
			{
				double auxvalue = bc->getComponent(i);
				dofValue.push_back(auxvalue);	
			}


//			if (bc->getComponent(i).size() != 0)
//			{
//				std::vector<double> valuei = bc ->getComponent(i);
//				if(valuei[0]>0.01)dof.push_back(2*nodeIndex + i);
//				
//			}


		}
	

	}
	return dofValue;
}

/*std::vector<double> BoundaryCondition::getComponent(const int& direction)
{
    if (direction == 0)
        return componentX_;
    else
        return componentY_;
}
*/


void FluidDomain::addBoundaryConditions()
{
	//transfering neumann conditions
	for (auto gbc : geometry_->getBoundaryConditions("NEUMANN"))
	{
		//transfering point loads applied to points
		if (gbc->getLineName().length() == 0)
		{
			std::string pointName = gbc->getPointName();
			std::string referenceSystem = gbc->getReferenceSystem();
			double componentX = gbc->getComponentX();
			double componentY = gbc->getComponentY();
			bool restrictedX = gbc->getRestrictedX();
			bool restrictedY = gbc->getRestrictedY(); 
			std::string method = gbc->getMethod();
			double penaltyParameter = gbc->getPenaltyParameter();
			std::string type = "NEUMANN_POINT";
			if (boundaryConditions_.count(type) == 0)
				boundaryConditions_[type] = std::vector<BoundaryCondition*>();
			int index = boundaryConditions_[type].size();
			Point* point = geometry_->getPoint(pointName);
			Node* node = point->getPointNode();
			BoundaryCondition* bc = new BoundaryCondition(index, node->getIndex(), restrictedX, restrictedY, componentX, componentY, method, referenceSystem, penaltyParameter);
			boundaryConditions_[type].push_back(bc);
		}

		if (gbc->getPointName().length() == 0)
		{
			std::string lineName = gbc->getLineName();
			std::string referenceSystem = gbc->getReferenceSystem();
			double componentX = gbc->getComponentX();
			double componentY = gbc->getComponentY();
			bool restrictedX = gbc->getRestrictedX();
			bool restrictedY = gbc->getRestrictedY();
			std::string method = gbc->getMethod();
			double penaltyParameter = gbc->getPenaltyParameter();
			std::string type = "NEUMANN_POINT";
			if (boundaryConditions_.count(type) == 0)
				boundaryConditions_[type] = std::vector<BoundaryCondition*>();
			int index = boundaryConditions_[type].size();
			Line* line = geometry_->getLine(lineName);
			for (Node* node : line->getLineNodes())
			{
			  BoundaryCondition* bc = new BoundaryCondition(index, node->getIndex(), restrictedX, restrictedY, componentX, componentY, method, referenceSystem, penaltyParameter);
			  boundaryConditions_[type].push_back(bc);
			  index++;
			}
		}
		//transfering distributed loads
		// else
		// {
		// 	std::string lineName = gbc->getLineName();
		// 	std::string referenceSystem = gbc->getReferenceSystem();
		// 	std::vector<double> componentX = gbc->getComponentX();
		// 	std::vector<double> componentY = gbc->getComponentY();
		// 	std::string method = gbc->getMethod();
		// 	double penaltyParameter = gbc->getPenaltyParameter();
		// 	std::string type = "DIST_LOAD";
		// 	if (boundaryConditions_.count(type) == 0)
		// 		boundaryConditions_[type] = std::vector<BoundaryCondition*>();
		// 	int index = boundaryConditions_[type].size();
		// 	Line* line = geometry_->getLine(lineName);
		// 	Node* node = point->getPointNode();
		// 	BoundaryCondition* bc = new BoundaryCondition(index, node->getIndex(), componentX, componentY, method, referenceSystem, penaltyParameter);
		// 	boundaryConditions_[type].push_back(bc);
		// }
	}

	for (auto gbc : geometry_->getBoundaryConditions("DIRICHLET"))
	{
		//transfering point displacements applied to points
		if (gbc->getLineName().length() == 0)
		{
			std::string pointName = gbc->getPointName();
			std::string referenceSystem = gbc->getReferenceSystem();
			double componentX = gbc->getComponentX();
			double componentY = gbc->getComponentY();
			bool restrictedX = gbc->getRestrictedX();
			bool restrictedY = gbc->getRestrictedY();
			std::string method = gbc->getMethod();
			double penaltyParameter = gbc->getPenaltyParameter();
			std::string type = "POINT_DISPLACEMENT";
			if (boundaryConditions_.count(type) == 0)
				boundaryConditions_[type] = std::vector<BoundaryCondition*>();
			int index = boundaryConditions_[type].size();
			Point* point = geometry_->getPoint(pointName);
			Node* node = point->getPointNode();
			BoundaryCondition* bc = new BoundaryCondition(index, node->getIndex(), restrictedX, restrictedY, componentX, componentY, method, referenceSystem, penaltyParameter);
			boundaryConditions_[type].push_back(bc);
		}
		//transfering point displacements applied to lines
		if (gbc->getPointName().length() == 0)
		{
			std::string lineName = gbc->getLineName();
			std::string referenceSystem = gbc->getReferenceSystem();
			double componentX = gbc->getComponentX();
			double componentY = gbc->getComponentY();
			bool restrictedX = gbc->getRestrictedX();
			bool restrictedY = gbc->getRestrictedY();
			std::string method = gbc->getMethod();
			double penaltyParameter = gbc->getPenaltyParameter();
			std::string type = "POINT_DISPLACEMENT";
			if (boundaryConditions_.count(type) == 0)
				boundaryConditions_[type] = std::vector<BoundaryCondition*>();
			int index = boundaryConditions_[type].size();
			Line* line = geometry_->getLine(lineName);
			for (Node* node : line->getLineNodes())
			{
			  BoundaryCondition* bc = new BoundaryCondition(index, node->getIndex(), restrictedX, restrictedY, componentX, componentY, method, referenceSystem, penaltyParameter);
				boundaryConditions_[type].push_back(bc);
				index++;
			}
		}

	}
		for (auto gbc : geometry_->getBoundaryConditions("MESHDISPLACEMENT"))
	{
		//transfering point displacements applied to points
		if (gbc->getLineName().length() == 0)
		{
			std::string pointName = gbc->getPointName();
			std::string referenceSystem = gbc->getReferenceSystem();
			double componentX = gbc->getComponentX();
			double componentY = gbc->getComponentY();
			bool restrictedX = gbc->getRestrictedX();
			bool restrictedY = gbc->getRestrictedY();
			std::string method = gbc->getMethod();
			double penaltyParameter = gbc->getPenaltyParameter();
			std::string type = "MESH_DISPLACEMENT";
			if (boundaryConditions_.count(type) == 0)
				boundaryConditions_[type] = std::vector<BoundaryCondition*>();
			int index = boundaryConditions_[type].size();
			Point* point = geometry_->getPoint(pointName);
			Node* node = point->getPointNode();
			BoundaryCondition* bc = new BoundaryCondition(index, node->getIndex(), restrictedX, restrictedY, componentX, componentY, method, referenceSystem, penaltyParameter);
			boundaryConditions_[type].push_back(bc);
		}
		//transfering point displacements applied to lines
		if (gbc->getPointName().length() == 0)
		{
			std::string lineName = gbc->getLineName();
			std::string referenceSystem = gbc->getReferenceSystem();
			double componentX = gbc->getComponentX();
			double componentY = gbc->getComponentY();
			bool restrictedX = gbc->getRestrictedX();
			bool restrictedY = gbc->getRestrictedY();
			std::string method = gbc->getMethod();
			double penaltyParameter = gbc->getPenaltyParameter();
			std::string type = "MESH_DISPLACEMENT";
			if (boundaryConditions_.count(type) == 0)
				boundaryConditions_[type] = std::vector<BoundaryCondition*>();
			int index = boundaryConditions_[type].size();
			Line* line = geometry_->getLine(lineName);
			for (Node* node : line->getLineNodes())
			{
			  BoundaryCondition* bc = new BoundaryCondition(index, node->getIndex(), restrictedX, restrictedY, componentX, componentY, method, referenceSystem, penaltyParameter);
				boundaryConditions_[type].push_back(bc);
				index++;
			}
		}

	}



}

void FluidDomain::generateMesh(const std::string& elementType, const std::string& algorithm, std::string geofile, const std::string& gmshPath, const bool& plotMesh, const bool& showInfo)
{	

    	int rank, size;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    
    	std::pair<std::string, bool> pair; pair.second = false;
  	if (rank == 0)
  	{
		pair = createMesh(geometry_, elementType, algorithm, geofile, gmshPath, plotMesh, showInfo);

    
		for (int i = 1; i < size; i++)
		{
			MPI_Send(pair.first.c_str(), pair.first.length()+1, MPI_CHAR, i, 0, PETSC_COMM_WORLD);
			if (i == size-1)
			{
				MPI_Send(&pair.second, 1, MPI_C_BOOL, i, 1, PETSC_COMM_WORLD);
				pair.second = false;
			}
		}
	}
	else
	{
		MPI_Status status;
		MPI_Probe(0, 0, PETSC_COMM_WORLD, &status);
		int count;
		MPI_Get_count(&status, MPI_CHAR, &count);
		char buf[count+1];
		MPI_Recv(&buf, count+1, MPI_CHAR, 0, 0, PETSC_COMM_WORLD, &status);
		pair.first = buf;
		if (rank == size-1)
			MPI_Recv(&pair.second, 1, MPI_C_BOOL, 0, 1, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	readInput(pair.first, pair.second);

	addBoundaryConditions();

 	 domainDecompositionMETIS();
}



//*****************************************************
// Solução do problema de escoamento compressível
//*****************************************************
void FluidDomain::solveCompressibleFlow()
{
	//	inicia loop no tempo
    std::cout<<"Starting time loop"<<std::endl; 
	
	Mat               A;  
    Vec               b, x, All, lumpedMass, xlump, AllLump, rhsSG, AllrhsSG;
    PetscErrorCode    ierr;
    PetscInt          Istart, Iend, Idof, Ione, iterations, *dof;
    KSP               ksp;
    PC                pc;
    VecScatter        ctx;
	VecScatter        ctxlump;
	VecScatter        ctxrhsSG;
    PetscScalar       val, value;
	MatNullSpace      nullsp;

	int rank;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	int matIndex = materials_.size();

	double rhoInf = materials_[0]->getDensity();
	double TempInf = materials_[0]->getUndisturbedTemperature();
	double calorv = materials_[0]->getSpecificHeatv();
	double calorp = materials_[0]->getSpecificHeatp();

	bounded_vector<double,2> uprescribed_;
	for (Node* n : nodes_)
		{
		 //	n->setCurrentMeshVelocity(uMeshIn_); 
			bounded_vector<double, 2> coord = n->getCurrentCoordinate();
//	********************************************************
//	Modifications to manually introduce specific initial conditions (example 1)
//	********************************************************
			
			
	/*		double enerInf = TempInf*calorv+.5*(uinf_(0)*uinf_(0)+uinf_(1)*uinf_(1)); 
			//n->setCurrentVelocity(uinf_);
			
			if(coord(0)<5.){
				rhoInf = 11./(enerInf*(calorp/calorv-1.));
				n->setCurrentDensity(rhoInf);
				
			}
			if(coord(0)>4.999999999){
				rhoInf = 10./(enerInf*(calorp/calorv-1.));
				n->setCurrentDensity(rhoInf);
			}
			
			double pressInf = (calorp-calorv)*rhoInf*TempInf;
			
	*/		
			

//	********************************************************
//	Modifications to manually introduce specific initial conditions (example 2 & example 3)
//	********************************************************
			n->setCurrentVelocity(uinf_);
		
			
			double pressInf = (calorp-calorv)*rhoInf*TempInf;
			double enerInf = TempInf*calorv+.5*(uinf_(0)*uinf_(0)+uinf_(1)*uinf_(1)); 
			
			
			

//	********************************************************
//	end Modifications to manually introduce specific initial conditions
//	********************************************************

			
			n->setCurrentPressure(pressInf);
			n->setCurrentVelocity(uinf_);
			
			//n->setCurrentMomentum(0.);
			n->setCurrentInternalEnergy(rhoInf*enerInf);			
			bounded_vector<double, 2> momentumInf = uinf_*rhoInf;
			n->setCurrentMomentum(momentumInf);
			n->setCurrentDensity(rhoInf);
			n->setCurrentTemperature(TempInf);
			
			
			//	********************************************************
			//	Modifications to manually introduce specific initial conditions (example 2)
			//	********************************************************

			if(coord(1)<.20001){
			if((coord(0)>.5999) && (coord(0)<0.6001)){
				uprescribed_(0)=0.;
				uprescribed_(1)=0.;
				n->setCurrentVelocity(uprescribed_);
				n->setCurrentMomentum(uprescribed_);
			}
			}

			//	********************************************************
			//	end Modifications to manually introduce specific initial conditions
			//	********************************************************

			
			//	********************************************************
			//	Modifications to manually introduce specific initial conditions (example 3)
			//	********************************************************

			/*if(coord(0)>0.29999){

				if (coord(1)<0.00001){
					n->setCurrentTemperature(0.555555556);
					uprescribed_(0)=0.;
					uprescribed_(1)=0.;
					n->setCurrentVelocity(uprescribed_);
					n->setCurrentMomentum(uprescribed_);
				}
			}*/

			//	********************************************************
			//	end Modifications to manually introduce specific initial conditions
			//	********************************************************
			
			//	********************************************************
			//	Modifications to manually introduce specific initial conditions (example 4)
			//	********************************************************

			
			/*if (coord(1)>0.99999){
				n->setCurrentTemperature(100.);
			}*/
			

			//	********************************************************
			//	end Modifications to manually introduce specific initial conditions
			//	********************************************************
		}

	exportToParaview(0);

	
			//if(coord(1)<.00000001){
//				bounded_vector<double, 2> rouprescrito;
//				rouprescrito(0)=0.;
//				rouprescrito(1)=0.;
			//	n->setCurrentInternalEnergy(rhoInf*calorv*(TempInf+50.));


	//			n->setCurrentMomentum(rouprescrito);
	//			n->setCurrentVelocity(rouprescrito);

			
/*			if(coord(1)>.999999){
				bounded_vector<double, 2> rouprescrito;
				rouprescrito(0)=0.;
				rouprescrito(1)=0.;
				n->setCurrentMomentum(rouprescrito);
				n->setCurrentVelocity(rouprescrito);

			}
*/
				
	//Array of constrained degrees of freedom
	getNeummanConstrains();
	std::vector<int> temp = getConstrains(); 
	PetscMalloc1(temp.size(),&dof);
	//std::cout<<" restriced dofs "<<" "<<temp.size()<<std::endl;



	for (size_t i = 0; i < temp.size(); i++)
	{
		dof[i] = temp[i];

	}
	  




	for (int itimeStep = 0; itimeStep < numberOfTimeSteps_; itimeStep++)
	{
	    if(rank == 0) std::cout<<"Time step "<<itimeStep<<std::endl; 

		//	************************************
		//	Solves step1 (explicit)
		//	************************************
		//Create PETSc sparse parallel matrix
        ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                            2 * nodes_.size() , 2 * nodes_.size(),
                            100,NULL,300,NULL,&A); 
		//	CHKERRQ(ierr);
    
            ierr = MatGetOwnershipRange(A, &Istart, &Iend);//CHKERRQ(ierr);
			
            //Create PETSc vectors
          	ierr = VecCreate(PETSC_COMM_WORLD,&b); //CHKERRQ(ierr);
            ierr = VecSetSizes(b,PETSC_DECIDE,2*nodes_.size()); //CHKERRQ(ierr);
            ierr = VecSetFromOptions(b); //CHKERRQ(ierr);
            ierr = VecDuplicate(b,&x); //CHKERRQ(ierr);
            ierr = VecDuplicate(b,&All); //CHKERRQ(ierr);
			ierr = VecCreate(PETSC_COMM_WORLD,&lumpedMass);
            ierr = VecSetSizes(lumpedMass,PETSC_DECIDE,nodes_.size()); //CHKERRQ(ierr);
			ierr = VecSetFromOptions(lumpedMass); //CHKERRQ(ierr);
			ierr = VecDuplicate(lumpedMass,&xlump); //CHKERRQ(ierr);
            ierr = VecDuplicate(lumpedMass,&AllLump); //CHKERRQ(ierr);
			ierr = VecCreate(PETSC_COMM_WORLD,&rhsSG);
            ierr = VecSetSizes(rhsSG,PETSC_DECIDE,2*nodes_.size()); //CHKERRQ(ierr);
			ierr = VecSetFromOptions(rhsSG); //CHKERRQ(ierr);
			ierr = VecDuplicate(rhsSG,&AllrhsSG); //CHKERRQ(ierr);
	
		//std::cout<< "Passo 00" << std::endl;
		for (Element *el : elements_) // loop over elements
		{
			if (elementPartition_[el->getIndex()] == rank) // if element in rank
			{
				double  Area = el->getArea();
				double mlump=Area/3.;

				for (size_t i = 0; i < el->getNodes().size(); i++)
				{
					int inode = el->getNode(i)->getIndex();
								//std::cout<<"mlump "<<mlump<<std::endl;
					ierr = VecSetValues(lumpedMass, 1, &inode, &mlump, ADD_VALUES);
				}

				std::pair<bounded_vector<double, 6>, bounded_matrix<double, 6, 6> > elementMatrices;
				//std::cout << "Element " << el->getIndex() << "\n";
				elementMatrices = el->elementContributionsS1(qdif_, theta_, sutherland_);

				for (size_t i = 0; i < el->getNodes().size(); i++)
				{ 
					int idof = 2 * el->getNode(i)->getIndex();
					ierr = VecSetValues(b, 1, &idof, &elementMatrices.first(2*i), ADD_VALUES);

					idof = 2 * el->getNode(i)->getIndex() + 1;
					ierr = VecSetValues(b, 1, &idof, &elementMatrices.first(2*i+1), ADD_VALUES);
					
					for (size_t j = 0; j < el->getNodes().size(); j++)
					{

						int dof1 = 2 * el->getNode(i)->getIndex();
						int dof2 = 2 * el->getNode(j)->getIndex();
						ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2*i, 2*j), ADD_VALUES);

						dof1 = 2 * el->getNode(i)->getIndex();
						dof2 = 2 * el->getNode(j)->getIndex()+1;
						ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2*i, 2*j+1), ADD_VALUES);

						dof1 = 2 * el->getNode(i)->getIndex()+1;
						dof2 = 2 * el->getNode(j)->getIndex();
						ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2*i+1, 2*j), ADD_VALUES);

						dof1 = 2 * el->getNode(i)->getIndex()+1;
						dof2 = 2 * el->getNode(j)->getIndex()+1;
						ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2*i+1, 2*j+1), ADD_VALUES); 

					} // j node
				} // i node
			} // if element in rank
		
		} // loop over elements
	
	 	//Assemble matrices and vectors
    	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); //CHKERRQ(ierr);
    	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); //CHKERRQ(ierr);

    	ierr = VecAssemblyBegin(b); //CHKERRQ(ierr);
    	ierr = VecAssemblyEnd(b); //CHKERRQ(ierr);
	 	ierr = VecAssemblyBegin(lumpedMass); //CHKERRQ(ierr);
    	ierr = VecAssemblyEnd(lumpedMass); //CHKERRQ(ierr);

		ierr = VecScatterCreateToAll(lumpedMass, &ctxlump, &AllLump);
		//CHKERRQ(ierr);
		ierr = VecScatterBegin(ctxlump, lumpedMass, AllLump, INSERT_VALUES, SCATTER_FORWARD);
		//CHKERRQ(ierr);
		ierr = VecScatterEnd(ctxlump, lumpedMass, AllLump, INSERT_VALUES, SCATTER_FORWARD);
		//CHKERRQ(ierr);
		ierr = VecScatterDestroy(&ctxlump);

			for (size_t i = 0; i < nodes_.size(); i++)
			{
				double mlump;
				Idof = i;

				ierr = VecGetValues(AllLump, Ione, &Idof, &val);
				//CHKERRQ(ierr);
				mlump = val;
	
				nodes_[i]->setLumpedMass(mlump); 

	
			}
			
		//MatView(A,PETSC_VIEWER_STDOUT_WORLD);//CHKERRQ(ierr);
		//VecView(b,PETSC_VIEWER_STDOUT_WORLD);//CHKERRQ(ierr);
		//VecView(x,PETSC_VIEWER_STDOUT_WORLD);//CHKERRQ(ierr);	
		//-------------------------------------------				
		//Applying boundary conditions on Momentum

		MatZeroRowsColumns(A, temp.size(), dof, 1.0, x, b);
		//End Applying boundary conditions on Momentum
		//-------------------------------------------
		
		//-------------------------------------------
		//Solving system.

		//Create KSP context to solve the linear system
		ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
		//CHKERRQ(ierr);
		ierr = KSPSetOperators(ksp, A, A);
		//CHKERRQ(ierr);
		//Solve using MUMPS
		#if defined(PETSC_HAVE_MUMPS)
			ierr = KSPSetType(ksp, KSPPREONLY);
			ierr = KSPGetPC(ksp, &pc);
			ierr = PCSetType(pc, PCLU);
		#endif
		ierr = KSPSetFromOptions(ksp);
		//CHKERRQ(ierr);
		ierr = KSPSetUp(ksp);

		//Solve linear system
		ierr = KSPSolve(ksp, b, x);
		//CHKERRQ(ierr);
		ierr = KSPGetTotalIterations(ksp, &iterations);

		//MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
		//VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
		//VecView(x,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

		//Gathers the solution vector to the master process
		ierr = VecScatterCreateToAll(x, &ctx, &All);
		//CHKERRQ(ierr);
		ierr = VecScatterBegin(ctx, x, All, INSERT_VALUES, SCATTER_FORWARD);
		//CHKERRQ(ierr);
		ierr = VecScatterEnd(ctx, x, All, INSERT_VALUES, SCATTER_FORWARD);
		//CHKERRQ(ierr);
		ierr = VecScatterDestroy(&ctx);
		//CHKERRQ(ierr);
		//End solving system
		//---------------------------------------------------------------
			
		//---------------------------------------------------------------
		//Updating momentum
		Ione = 1;

		for (size_t i = 0; i < nodes_.size(); i++)
		{
			bounded_vector<double,2> deltarhou;
			Idof = 2 * i;
			ierr = VecGetValues(All, Ione, &Idof, &val);
			//CHKERRQ(ierr);
			deltarhou(0) = val;

			Idof = 2 * i + 1;
			ierr = VecGetValues(All, Ione, &Idof, &val);
			//CHKERRQ(ierr);
			deltarhou(1) = val;

//	********************************************************
//	Modifications to manually introduce specific initial conditions (example 3)
//	********************************************************
			/*bounded_vector<double, 2> coord = nodes_[i]->getCurrentCoordinate(); 
			if(coord(0)>0.29999){

				if (coord(1)<0.00001){
					deltarhou(0) = 0.;
				}
			}*/
		
		

//	********************************************************
//	end Modifications to manually introduce specific initial conditions
//	********************************************************
						
		nodes_[i]->setDeltaMomentum(deltarhou); 
		}
		//End updating momentum
		//---------------------------------------------------------------
		ierr = KSPDestroy(&ksp); //CHKERRQ(ierr);
		ierr = VecDestroy(&b); //CHKERRQ(ierr);
		ierr = VecDestroy(&x); //CHKERRQ(ierr);
		ierr = VecDestroy(&All); //CHKERRQ(ierr);
		ierr = MatDestroy(&A); //CHKERRQ(ierr);

		MPI_Barrier(PETSC_COMM_WORLD);


	//****************************************
	// Small Speed Shock Capture for momentum
	//****************************************

		for (Element *el : elements_) // loop over elements
		{
			if (elementPartition_[el->getIndex()] == rank) // if element in rank
			{
				bounded_vector<double, 6> rhsS;
				//std::cout << "Element " << el->getIndex() << "\n";
				rhsS = el->elementShockLowSpeedS1();

				for (size_t i = 0; i < el->getNodes().size(); i++)
				{ 
					int idof = 2 * el->getNode(i)->getIndex();
						ierr = VecSetValues(rhsSG, 1, &idof, &rhsS(2*i), ADD_VALUES);
						idof = 2 * el->getNode(i)->getIndex() + 1;
						ierr = VecSetValues(rhsSG, 1, &idof, &rhsS(2*i+1), ADD_VALUES);
				} // i node
			} // if element in rank
		
		} // loop over elements


		ierr = VecAssemblyBegin(rhsSG); //CHKERRQ(ierr);
    	ierr = VecAssemblyEnd(rhsSG); //CHKERRQ(ierr);

		ierr = VecScatterCreateToAll(rhsSG, &ctxrhsSG, &AllrhsSG);
		//CHKERRQ(ierr);
		ierr = VecScatterBegin(ctxrhsSG, rhsSG, AllrhsSG, INSERT_VALUES, SCATTER_FORWARD);
		//CHKERRQ(ierr);
		ierr = VecScatterEnd(ctxrhsSG, rhsSG, AllrhsSG, INSERT_VALUES, SCATTER_FORWARD);
		//CHKERRQ(ierr);
		ierr = VecScatterDestroy(&ctxrhsSG);

		
		Ione = 1;
		for (size_t i = 0; i < nodes_.size(); i++)
			{
				double mlump=nodes_[i]->getLumpedMass();
				bounded_vector<double,2> drou1=nodes_[i]->getDeltaMomentum();
				bounded_vector<double,2> drouS;

				Idof = 2*i;

				ierr = VecGetValues(AllrhsSG, Ione, &Idof, &val);
				//CHKERRQ(ierr);
				double drou1s = (1./(1.+0.5*alpha_))*drou1(0)+(alpha_/(1.+0.5*alpha_))*(1./mlump)*val;

				Idof = 2*i+1;

				ierr = VecGetValues(AllrhsSG, Ione, &Idof, &val);
				//CHKERRQ(ierr);
				double drou2s = (1./(1.+0.5*alpha_))*drou1(1)+(alpha_/(1.+0.5*alpha_))*(1./mlump)*val;
				drouS(0)=drou1s;
				drouS(1)=drou2s;

				//	********************************************************
				//	Modifications to manually introduce specific initial conditions (example 4)
				//	********************************************************

				/*drouS(0)=0.;
				drouS(1)=0.;*/
					

				//	********************************************************
				//	end Modifications to manually introduce specific initial conditions
				//	********************************************************
				nodes_[i]->setDeltaMomentum(drouS);
	
			}


		ierr = VecDestroy(&rhsSG); //CHKERRQ(ierr);
		ierr = VecDestroy(&AllrhsSG); //CHKERRQ(ierr);


		MPI_Barrier(PETSC_COMM_WORLD);


	//*********************************
	//	Solves step2 (implicit or explicit)
	//*********************************
	
		//Create PETSc sparse parallel matrix
        ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                            nodes_.size() , nodes_.size(),
                            100,NULL,300,NULL,&A); 

		//	CHKERRQ(ierr);
            
         ierr = MatGetOwnershipRange(A, &Istart, &Iend);//CHKERRQ(ierr);
		             
        //Create PETSc vectors
        ierr = VecCreate(PETSC_COMM_WORLD,&b); //CHKERRQ(ierr);
        ierr = VecSetSizes(b,PETSC_DECIDE,nodes_.size()); //CHKERRQ(ierr);
        ierr = VecSetFromOptions(b); //CHKERRQ(ierr);
        ierr = VecDuplicate(b,&x); //CHKERRQ(ierr);
        ierr = VecDuplicate(b,&All); //CHKERRQ(ierr); 
		ierr = VecCreate(PETSC_COMM_WORLD,&rhsSG);
        ierr = VecSetSizes(rhsSG,PETSC_DECIDE,nodes_.size()); //CHKERRQ(ierr);
		ierr = VecSetFromOptions(rhsSG); //CHKERRQ(ierr);
		ierr = VecDuplicate(rhsSG,&AllrhsSG); //CHKERRQ(ierr);

		//std::cout<< "Passo 00" << std::endl;
		for (Element *el : elements_) // loop over elements
		{
			
			// Integração numérica das matrizes locais e montagem na matriz global
			// Integhação numérica do vetor incógnita e montagem no vetor global

			if (elementPartition_[el->getIndex()] == rank) // if element in rank
			{

				std::pair<bounded_vector<double, 3>, bounded_matrix<double, 3, 3> > elementMatrices;
				//std::cout << "Element " << el->getIndex() << "\n";
				elementMatrices = el->elementContributionsS2(qdif_, theta_, sutherland_);

				for (size_t i = 0; i < el->getNodes().size(); i++)
				{ 
					int idof = el->getNode(i)->getIndex();

					ierr = VecSetValues(b, 1, &idof, &elementMatrices.first(i), ADD_VALUES);

					for (size_t j = 0; j < el->getNodes().size(); j++)
						{
							int dof1 = el->getNode(i)->getIndex();
							int dof2 = el->getNode(j)->getIndex();
							ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(i, j), ADD_VALUES);
						}
				}
			}
				
		}

		// resolver sistema linear

	 		//Assemble matrices and vectors
    		ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); //CHKERRQ(ierr);
    		ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); //CHKERRQ(ierr);

    		ierr = VecAssemblyBegin(b); //CHKERRQ(ierr);
    		ierr = VecAssemblyEnd(b); //CHKERRQ(ierr);

//			MatView(A,PETSC_VIEWER_STDOUT_WORLD);//CHKERRQ(ierr);
//			VecView(b,PETSC_VIEWER_STDOUT_WORLD);//CHKERRQ(ierr);
//			VecView(x,PETSC_VIEWER_STDOUT_WORLD);//CHKERRQ(ierr);	
			
			//aplicando condição de contorno

			//MatZeroRowsColumns(A, temp.size(), dof, 1.0, x, b);

			//Create KSP context to solve the linear system
			ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
			//CHKERRQ(ierr);
			ierr = KSPSetOperators(ksp, A, A);
			//CHKERRQ(ierr);
			//Solve using MUMPS
			#if defined(PETSC_HAVE_MUMPS)
			ierr = KSPSetType(ksp, KSPPREONLY);
			ierr = KSPGetPC(ksp, &pc);
			ierr = PCSetType(pc, PCLU);
			#endif
			ierr = KSPSetFromOptions(ksp);
			//CHKERRQ(ierr);
			ierr = KSPSetUp(ksp);

			//Solve linear system
			ierr = KSPSolve(ksp, b, x);
			//CHKERRQ(ierr);
			ierr = KSPGetTotalIterations(ksp, &iterations);

			//MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
			//VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
			//VecView(x,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

			//Gathers the solution vector to the master process
			ierr = VecScatterCreateToAll(x, &ctx, &All);
			//CHKERRQ(ierr);
			ierr = VecScatterBegin(ctx, x, All, INSERT_VALUES, SCATTER_FORWARD);
			//CHKERRQ(ierr);
			ierr = VecScatterEnd(ctx, x, All, INSERT_VALUES, SCATTER_FORWARD);
			//CHKERRQ(ierr);
			ierr = VecScatterDestroy(&ctx);
			//CHKERRQ(ierr);

			Ione = 1;

			for (size_t i = 0; i < nodes_.size(); i++)
			{
				double deltarho;
				Idof = i;

				ierr = VecGetValues(All, Ione, &Idof, &val);
				//CHKERRQ(ierr);
				deltarho = val;

				//	********************************************************
				//	Modifications to manually introduce specific initial conditions (example 4)
				//	********************************************************

				//deltarho = 0.;
					

				//	********************************************************
				//	end Modifications to manually introduce specific initial conditions
				//	********************************************************

				nodes_[i]->setDeltaDensity(deltarho); 

				//	********************************************************
				//	Modifications to manually introduce specific initial conditions (example 3)
				//	********************************************************

				/*bounded_vector<double, 2> coord = nodes_[i]->getCurrentCoordinate();
				if(coord(0)<.000001){

					nodes_[i]->setDeltaDensity(0.);

				}*/

				//	********************************************************
				//	end Modifications to manually introduce specific initial conditions
				//	********************************************************
			}

			// calcular o erro (norma do vetor correção)
			ierr = KSPDestroy(&ksp); //CHKERRQ(ierr);
			ierr = VecDestroy(&b); //CHKERRQ(ierr);
			ierr = VecDestroy(&x); //CHKERRQ(ierr);
			ierr = VecDestroy(&All); //CHKERRQ(ierr);
			ierr = MatDestroy(&A); //CHKERRQ(ierr);
			// 	std::cout<<"Iteration "<<iterationStep2<<std::endl;

			MPI_Barrier(PETSC_COMM_WORLD);

	//****************************************
	// Small Speed Shock Capture for density
	//****************************************

		for (Element *el : elements_) // loop over elements
		{
			if (elementPartition_[el->getIndex()] == rank) // if element in rank
			{
			
				bounded_vector<double, 3> rhsS;
				//std::cout << "Element " << el->getIndex() << "\n";
				rhsS = el->elementShockLowSpeedS2();

				for (size_t i = 0; i < el->getNodes().size(); i++)
				{ 
					int idof = el->getNode(i)->getIndex();
						ierr = VecSetValues(rhsSG, 1, &idof, &rhsS(i), ADD_VALUES);
						
				} // i node
			} // if element in rank
		
		} // loop over elements


		ierr = VecAssemblyBegin(rhsSG); //CHKERRQ(ierr);
    	ierr = VecAssemblyEnd(rhsSG); //CHKERRQ(ierr);

		ierr = VecScatterCreateToAll(rhsSG, &ctxrhsSG, &AllrhsSG);
		//CHKERRQ(ierr);
		ierr = VecScatterBegin(ctxrhsSG, rhsSG, AllrhsSG, INSERT_VALUES, SCATTER_FORWARD);
		//CHKERRQ(ierr);
		ierr = VecScatterEnd(ctxrhsSG, rhsSG, AllrhsSG, INSERT_VALUES, SCATTER_FORWARD);
		//CHKERRQ(ierr);
		ierr = VecScatterDestroy(&ctxrhsSG);
		//Assembly. Scatter... rhSG

		
		Ione = 1;
		for (size_t i = 0; i < nodes_.size(); i++)
			{
				double mlump=nodes_[i]->getLumpedMass();
				double drho1=nodes_[i]->getDeltaDensity();
				double drhoS;

				Idof = i;

				ierr = VecGetValues(AllrhsSG, Ione, &Idof, &val);
				//CHKERRQ(ierr);
				drhoS = (1./(1.+0.5*alpha_))*drho1+(alpha_/(1.+0.5*alpha_))*(1./mlump)*val;

				

				nodes_[i]->setDeltaDensity(drhoS);
	
			}

		ierr = VecDestroy(&rhsSG); //CHKERRQ(ierr);
		ierr = VecDestroy(&AllrhsSG); //CHKERRQ(ierr);

		MPI_Barrier(PETSC_COMM_WORLD);


	//***********************************************************************
	//	Solves step3 (explicit)
	//***********************************************************************
		//Create PETSc sparse parallel matrix
        ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                            nodes_.size() , nodes_.size(),
                            100,NULL,300,NULL,&A); 
		//	CHKERRQ(ierr);
            
        ierr = MatGetOwnershipRange(A, &Istart, &Iend);//CHKERRQ(ierr);
            
        //Create PETSc vectors
        ierr = VecCreate(PETSC_COMM_WORLD,&b); //CHKERRQ(ierr);
        ierr = VecSetSizes(b,PETSC_DECIDE,nodes_.size()); //CHKERRQ(ierr);
        ierr = VecSetFromOptions(b); //CHKERRQ(ierr);
        ierr = VecDuplicate(b,&x); //CHKERRQ(ierr);
        ierr = VecDuplicate(b,&All); //CHKERRQ(ierr); 
		ierr = VecCreate(PETSC_COMM_WORLD,&rhsSG);
        ierr = VecSetSizes(rhsSG,PETSC_DECIDE,nodes_.size()); //CHKERRQ(ierr);
		ierr = VecSetFromOptions(rhsSG); //CHKERRQ(ierr);
		ierr = VecDuplicate(rhsSG,&AllrhsSG); //CHKERRQ(ierr);
		
		//std::cout<< "Passo 00" << std::endl;
		for (Element *el : elements_) // loop over elements
		{
			
			// Integração numérica das matrizes locais e montagem na matriz global
			// Integhação numérica do vetor incógnita e montagem no vetor global

			if (elementPartition_[el->getIndex()] == rank) // if element in rank
					{

						std::pair<bounded_vector<double, 3>, bounded_matrix<double, 3, 3> > elementMatrices;
						//std::cout << "Element " << el->getIndex() << "\n";
						elementMatrices = el->elementContributionsS3(qdif_, theta_, sutherland_);

						for (size_t i = 0; i < el->getNodes().size(); i++)
						{ 
							int idof = el->getNode(i)->getIndex();
							ierr = VecSetValues(b, 1, &idof, &elementMatrices.first(i), ADD_VALUES);

			
							for (size_t j = 0; j < el->getNodes().size(); j++)
							{
								int dof1 = el->getNode(i)->getIndex();
								int dof2 = el->getNode(j)->getIndex();
								ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(i, j), ADD_VALUES);



							}
							

						}
					}
				
			
		}

			// resolver sistema linear

	 		//Assemble matrices and vectors
    		ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); //CHKERRQ(ierr);
    		ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); //CHKERRQ(ierr);
	
    		ierr = VecAssemblyBegin(b); //CHKERRQ(ierr);
    		ierr = VecAssemblyEnd(b); //CHKERRQ(ierr);

			//MatView(A,PETSC_VIEWER_STDOUT_WORLD);//CHKERRQ(ierr);
			//VecView(b,PETSC_VIEWER_STDOUT_WORLD);//CHKERRQ(ierr);
			//VecView(x,PETSC_VIEWER_STDOUT_WORLD);//CHKERRQ(ierr);

			//aplicando condição de contorno

			//MatZeroRowsColumns(A, temp.size(), dof, 1.0, x, b);

			//	resolver o sistema aqui.

			//Create KSP context to solve the linear system
			ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
			//CHKERRQ(ierr);
			ierr = KSPSetOperators(ksp, A, A);
			//CHKERRQ(ierr);
			//Solve using MUMPS
			#if defined(PETSC_HAVE_MUMPS)
			ierr = KSPSetType(ksp, KSPPREONLY);
			ierr = KSPGetPC(ksp, &pc);
			ierr = PCSetType(pc, PCLU);
			#endif
			ierr = KSPSetFromOptions(ksp);
			//CHKERRQ(ierr);
			ierr = KSPSetUp(ksp);

			//Solve linear system
			ierr = KSPSolve(ksp, b, x);

			//CHKERRQ(ierr);
			ierr = KSPGetTotalIterations(ksp, &iterations);

			//MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
			//VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
			//VecView(x,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

			//Gathers the solution vector to the master process
			ierr = VecScatterCreateToAll(x, &ctx, &All);
			//CHKERRQ(ierr);
			ierr = VecScatterBegin(ctx, x, All, INSERT_VALUES, SCATTER_FORWARD);
			//CHKERRQ(ierr);
			ierr = VecScatterEnd(ctx, x, All, INSERT_VALUES, SCATTER_FORWARD);
			//CHKERRQ(ierr);
			ierr = VecScatterDestroy(&ctx);
			//CHKERRQ(ierr);



			Ione = 1;

			for (size_t i = 0; i < nodes_.size(); i++)
			{
				double deltarhoe;
				Idof = i;
				ierr = VecGetValues(All, Ione, &Idof, &val);
				//CHKERRQ(ierr);
				deltarhoe = val;

				nodes_[i]->setDeltaEnergy(deltarhoe);
				
				//	********************************************************
				//	Modifications to manually introduce specific initial conditions (example 3)
				//	********************************************************

				/*bounded_vector<double, 2> coord = nodes_[i]->getCurrentCoordinate(); 
				if(coord(0)<0.0001){

				nodes_[i]->setDeltaEnergy(0.);

				}*/

				//	********************************************************
				//	end Modifications to manually introduce specific initial conditions
				//	********************************************************

			}

			


			// calcular o erro (norma do vetor correção)
			ierr = KSPDestroy(&ksp); //CHKERRQ(ierr);
			ierr = VecDestroy(&b); //CHKERRQ(ierr);
			ierr = VecDestroy(&x); //CHKERRQ(ierr);
			ierr = VecDestroy(&All); //CHKERRQ(ierr);
			ierr = MatDestroy(&A); //CHKERRQ(ierr);
			// 	std::cout<<"Iteration "<<iterationStep2<<std::endl;

		MPI_Barrier(PETSC_COMM_WORLD);


	//****************************************
	// Small Speed Shock Capture for energy
	//****************************************

		for (Element *el : elements_) // loop over elements
		{
			if (elementPartition_[el->getIndex()] == rank) // if element in rank
			{
			
				bounded_vector<double, 3> rhsS;
				//std::cout << "Element " << el->getIndex() << "\n";
				//rhsS = el->elementShockLowSpeedS3();

				for (size_t i = 0; i < el->getNodes().size(); i++)
				{ 
					double valuesS=0.;
					int idof = el->getNode(i)->getIndex();
						ierr = VecSetValues(rhsSG, 1, &idof, &valuesS, ADD_VALUES);
						
				} // i node
			} // if element in rank
		
		} // loop over elements

		ierr = VecAssemblyBegin(rhsSG); //CHKERRQ(ierr);
    	ierr = VecAssemblyEnd(rhsSG); //CHKERRQ(ierr);

		ierr = VecScatterCreateToAll(rhsSG, &ctxrhsSG, &AllrhsSG);
		//CHKERRQ(ierr);
		ierr = VecScatterBegin(ctxrhsSG, rhsSG, AllrhsSG, INSERT_VALUES, SCATTER_FORWARD);
		//CHKERRQ(ierr);
		ierr = VecScatterEnd(ctxrhsSG, rhsSG, AllrhsSG, INSERT_VALUES, SCATTER_FORWARD);
		//CHKERRQ(ierr);
		ierr = VecScatterDestroy(&ctxrhsSG);
		//Assembly. Scatter... rhSG
		
		Ione = 1;
		for (size_t i = 0; i < nodes_.size(); i++)
			{
				double mlump=nodes_[i]->getLumpedMass();
				double drhoe1=nodes_[i]->getDeltaEnergy();
				double drhoeS;

				Idof = i;

				ierr = VecGetValues(AllrhsSG, Ione, &Idof, &val);
				//CHKERRQ(ierr);
				drhoeS = (1./(1.+0.5*alpha_))*drhoe1+(alpha_/(1.+0.5*alpha_))*(1./mlump)*val;

				

				nodes_[i]->setDeltaEnergy(drhoeS);
	
			}


		ierr = VecDestroy(&rhsSG); //CHKERRQ(ierr);
		ierr = VecDestroy(&AllrhsSG); //CHKERRQ(ierr);



		//update variables
		int matIndex = 0; //materials_.size();
		for (Node* n : nodes_)
		{
			n->updateVariables();
			
 
		}


	

		for (Node* n : nodes_)
		{			
			n->computeSecondaryVariables(materials_[matIndex] ->getSpecificHeatv(), materials_[matIndex]->getSpecificHeatp());

			//	********************************************************
			//	Modifications to manually introduce specific initial conditions (example 3)
			//	********************************************************

			/*bounded_vector<double, 2> coord = n->getCurrentCoordinate(); 
			if(coord(0)>0.29999){

				if (coord(1)<0.00001){
					n->setCurrentTemperature(0.555555556);
				}
			}*/

			//	********************************************************
			//	end Modifications to manually introduce specific initial conditions
			//	********************************************************

			//	********************************************************
			//	Modifications to manually introduce specific initial conditions (example 4)
			//	********************************************************

			/*bounded_vector<double, 2> coord = n->getCurrentCoordinate(); 
			if (coord(1)>0.99999){
				n->setCurrentTemperature(100.);
			}*/
			

			//	********************************************************
			//	end Modifications to manually introduce specific initial conditions
			//	********************************************************
		}
	
			ierr = VecDestroy(&lumpedMass); //CHKERRQ(ierr);
			ierr = VecDestroy(&xlump); //CHKERRQ(ierr);
			ierr = VecDestroy(&AllLump); //CHKERRQ(ierr);
			exportToParaview(itimeStep+1);
	
}

}


//	************************************
//	SOLVE COMPRESSIBLE FLOW WITH MOVING BOUNDARY
//	************************************

void FluidDomain::solveCompressibleFlowMoving()
{
	//	inicia loop no tempo
    std::cout<<"Starting time loop"<<std::endl; 
	
	Mat               A;  
    Vec               b, x, All, lumpedMass, xlump, AllLump, rhsSG, AllrhsSG, dispMeshValues;
    PetscErrorCode    ierr;
    PetscInt          Istart, Iend, Idof, Ione, iterations, *dof, *dofMesh;
    KSP               ksp;
    PC                pc;
    VecScatter        ctx;
	VecScatter        ctxlump;
	VecScatter        ctxrhsSG;
    PetscScalar       val, value;
	MatNullSpace      nullsp;

	int rank;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	int matIndex = materials_.size();

	double rhoInf = materials_[0]->getDensity();
	double TempInf = materials_[0]->getUndisturbedTemperature();
	double calorv = materials_[0]->getSpecificHeatv();
	double calorp = materials_[0]->getSpecificHeatp();

	/*// Set element mesh moving parameters
    double vMax = 0., vMin = 1.e10;
    for (Element* e : elements_){
        double v = e -> getJacobian();
        if (v > vMax) vMax = v;
        if (v < vMin) vMin = v;
    };
    for (int i = 0; i < numElem; i++){
        double v = elements_[i] -> getJacobian();
        double eta = 1 + (1. - vMin / vMax) / (v / vMax);
        elements_[i] -> setMeshMovingParameter(eta);
    };*/

	bounded_vector<double,2> uprescribed_;
	for (Node* n : nodes_)
		{
			n->setCurrentMeshVelocity(uMeshIn_);
			bounded_vector<double, 2> coord = n->getCurrentCoordinate();
//	********************************************************
//	Modifications to manually introduce specific initial conditions (example 1)
//	********************************************************
			
			
	/*		double enerInf = TempInf*calorv+.5*(uinf_(0)*uinf_(0)+uinf_(1)*uinf_(1)); 
			//n->setCurrentVelocity(uinf_);
			
			if(coord(0)<5.){
				rhoInf = 11./(enerInf*(calorp/calorv-1.));
				n->setCurrentDensity(rhoInf);
				
			}
			if(coord(0)>4.999999999){
				rhoInf = 10./(enerInf*(calorp/calorv-1.));
				n->setCurrentDensity(rhoInf);
			}
			
			double pressInf = (calorp-calorv)*rhoInf*TempInf;
			
	*/		
			

//	********************************************************
//	Modifications to manually introduce specific initial conditions (example 2 & example 3)
//	********************************************************
			n->setCurrentVelocity(uinf_);
		
			
			double pressInf = (calorp-calorv)*rhoInf*TempInf;
			double enerInf = TempInf*calorv+.5*(uinf_(0)*uinf_(0)+uinf_(1)*uinf_(1)); 
			
			
			

//	********************************************************
//	end Modifications to manually introduce specific initial conditions
//	********************************************************

			
			n->setCurrentPressure(pressInf);
			n->setCurrentVelocity(uinf_);
			
			//n->setCurrentMomentum(0.);
			n->setCurrentInternalEnergy(rhoInf*enerInf);			
			bounded_vector<double, 2> momentumInf = uinf_*rhoInf;
			n->setCurrentMomentum(momentumInf);
			n->setCurrentDensity(rhoInf);
			n->setCurrentTemperature(TempInf);
			
			
			//	********************************************************
			//	Modifications to manually introduce specific initial conditions (example 2)
			//	********************************************************

			if(coord(1)<.20001){
			if((coord(0)>.5999) && (coord(0)<0.6001)){
				uprescribed_(0)=0.;
				uprescribed_(1)=0.;
				n->setCurrentVelocity(uprescribed_);
				n->setCurrentMomentum(uprescribed_);
			}
			}

			//	********************************************************
			//	end Modifications to manually introduce specific initial conditions
			//	********************************************************

			
			//	********************************************************
			//	Modifications to manually introduce specific initial conditions (example 3)
			//	********************************************************

			/*if(coord(0)>0.29999){

				if (coord(1)<0.00001){
					n->setCurrentTemperature(0.555555556);
					uprescribed_(0)=0.;
					uprescribed_(1)=0.;
					n->setCurrentVelocity(uprescribed_);
					n->setCurrentMomentum(uprescribed_);
				}
			}*/

			//	********************************************************
			//	end Modifications to manually introduce specific initial conditions
			//	********************************************************
			
			//	********************************************************
			//	Modifications to manually introduce specific initial conditions (example 4)
			//	********************************************************

			
			/*if (coord(1)>0.99999){
				n->setCurrentTemperature(100.);
			}*/
			

			//	********************************************************
			//	end Modifications to manually introduce specific initial conditions
			//	********************************************************
		}

	exportToParaview(0);

	
			//if(coord(1)<.00000001){
//				bounded_vector<double, 2> rouprescrito;
//				rouprescrito(0)=0.;
//				rouprescrito(1)=0.;
			//	n->setCurrentInternalEnergy(rhoInf*calorv*(TempInf+50.));


	//			n->setCurrentMomentum(rouprescrito);
	//			n->setCurrentVelocity(rouprescrito);

			
/*			if(coord(1)>.999999){
				bounded_vector<double, 2> rouprescrito;
				rouprescrito(0)=0.;
				rouprescrito(1)=0.;
				n->setCurrentMomentum(rouprescrito);
				n->setCurrentVelocity(rouprescrito);

			}
*/
				
	//Array of constrained degrees of freedom
	getNeummanConstrains();
	std::vector<int> temp = getConstrains(); 
	PetscMalloc1(temp.size(),&dof);
	//std::cout<<" restriced dofs "<<" "<<temp.size()<<std::endl;

	for (size_t i = 0; i < temp.size(); i++)
	{
		dof[i] = temp[i];

	}
	  
	std::vector<int> tempMesh = getMeshConstrains(); 
	PetscMalloc1(tempMesh.size(),&dofMesh);

	std::vector<double> tempMeshValues = getMeshDisplacementValues(); 

	//for (int i =0; i< tempMeshValues.size(); i++)
	//{
	//	std::cout<<i<<" "<<tempMesh[i]<<" consvalue "<<tempMeshValues[i]<<std::endl;		
	//}
	
	//PetscMalloc1(tempMeshValues.size(),&dispMeshValues);*/

	//std::cout<<" restriced dofs "<<" "<<temp.size()<<std::endl;

	for (size_t i = 0; i < tempMesh.size(); i++)
	{
		dofMesh[i] = tempMesh[i];

	}


	for (int itimeStep = 0; itimeStep < numberOfTimeSteps_; itimeStep++)
	{
	    if(rank == 0) std::cout<<"Time step "<<itimeStep<<std::endl; 

		//	************************************
		//	Solves step1 (explicit)
		//	************************************
		//Create PETSc sparse parallel matrix
        ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                            2 * nodes_.size() , 2 * nodes_.size(),
                            100,NULL,300,NULL,&A); 
		//	CHKERRQ(ierr);
    
            ierr = MatGetOwnershipRange(A, &Istart, &Iend);//CHKERRQ(ierr);
			
            //Create PETSc vectors
          	ierr = VecCreate(PETSC_COMM_WORLD,&b); //CHKERRQ(ierr);
            ierr = VecSetSizes(b,PETSC_DECIDE,2*nodes_.size()); //CHKERRQ(ierr);
            ierr = VecSetFromOptions(b); //CHKERRQ(ierr);
            ierr = VecDuplicate(b,&x); //CHKERRQ(ierr);
            ierr = VecDuplicate(b,&All); //CHKERRQ(ierr);
			ierr = VecCreate(PETSC_COMM_WORLD,&lumpedMass);
            ierr = VecSetSizes(lumpedMass,PETSC_DECIDE,nodes_.size()); //CHKERRQ(ierr);
			ierr = VecSetFromOptions(lumpedMass); //CHKERRQ(ierr);
			ierr = VecDuplicate(lumpedMass,&xlump); //CHKERRQ(ierr);
            ierr = VecDuplicate(lumpedMass,&AllLump); //CHKERRQ(ierr);
			ierr = VecCreate(PETSC_COMM_WORLD,&rhsSG);
            ierr = VecSetSizes(rhsSG,PETSC_DECIDE,2*nodes_.size()); //CHKERRQ(ierr);
			ierr = VecSetFromOptions(rhsSG); //CHKERRQ(ierr);
			ierr = VecDuplicate(rhsSG,&AllrhsSG); //CHKERRQ(ierr);
	
		//std::cout<< "Passo 00" << std::endl;
		for (Element *el : elements_) // loop over elements
		{
			if (elementPartition_[el->getIndex()] == rank) // if element in rank
			{
				double  Area = el->getArea();
				double mlump=Area/3.;

				for (size_t i = 0; i < el->getNodes().size(); i++)
				{
					int inode = el->getNode(i)->getIndex();
								//std::cout<<"mlump "<<mlump<<std::endl;
					ierr = VecSetValues(lumpedMass, 1, &inode, &mlump, ADD_VALUES);
				}

				std::pair<bounded_vector<double, 6>, bounded_matrix<double, 6, 6> > elementMatrices;
				//std::cout << "Element " << el->getIndex() << "\n";
				elementMatrices = el->elementContributionsS1(qdif_, theta_, sutherland_);

				for (size_t i = 0; i < el->getNodes().size(); i++)
				{ 
					int idof = 2 * el->getNode(i)->getIndex();
					ierr = VecSetValues(b, 1, &idof, &elementMatrices.first(2*i), ADD_VALUES);

					idof = 2 * el->getNode(i)->getIndex() + 1;
					ierr = VecSetValues(b, 1, &idof, &elementMatrices.first(2*i+1), ADD_VALUES);
					
					for (size_t j = 0; j < el->getNodes().size(); j++)
					{

						int dof1 = 2 * el->getNode(i)->getIndex();
						int dof2 = 2 * el->getNode(j)->getIndex();
						ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2*i, 2*j), ADD_VALUES);

						dof1 = 2 * el->getNode(i)->getIndex();
						dof2 = 2 * el->getNode(j)->getIndex()+1;
						ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2*i, 2*j+1), ADD_VALUES);

						dof1 = 2 * el->getNode(i)->getIndex()+1;
						dof2 = 2 * el->getNode(j)->getIndex();
						ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2*i+1, 2*j), ADD_VALUES);

						dof1 = 2 * el->getNode(i)->getIndex()+1;
						dof2 = 2 * el->getNode(j)->getIndex()+1;
						ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2*i+1, 2*j+1), ADD_VALUES); 

					} // j node
				} // i node
			} // if element in rank
		
		} // loop over elements
	
	 	//Assemble matrices and vectors
    	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); //CHKERRQ(ierr);
    	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); //CHKERRQ(ierr);

    	ierr = VecAssemblyBegin(b); //CHKERRQ(ierr);
    	ierr = VecAssemblyEnd(b); //CHKERRQ(ierr);
	 	ierr = VecAssemblyBegin(lumpedMass); //CHKERRQ(ierr);
    	ierr = VecAssemblyEnd(lumpedMass); //CHKERRQ(ierr);

		ierr = VecScatterCreateToAll(lumpedMass, &ctxlump, &AllLump);
		//CHKERRQ(ierr);
		ierr = VecScatterBegin(ctxlump, lumpedMass, AllLump, INSERT_VALUES, SCATTER_FORWARD);
		//CHKERRQ(ierr);
		ierr = VecScatterEnd(ctxlump, lumpedMass, AllLump, INSERT_VALUES, SCATTER_FORWARD);
		//CHKERRQ(ierr);
		ierr = VecScatterDestroy(&ctxlump);

			for (size_t i = 0; i < nodes_.size(); i++)
			{
				double mlump;
				Idof = i;

				ierr = VecGetValues(AllLump, Ione, &Idof, &val);
				//CHKERRQ(ierr);
				mlump = val;
	
				nodes_[i]->setLumpedMass(mlump); 

	
			}
			
		//MatView(A,PETSC_VIEWER_STDOUT_WORLD);//CHKERRQ(ierr);
		//VecView(b,PETSC_VIEWER_STDOUT_WORLD);//CHKERRQ(ierr);
		//VecView(x,PETSC_VIEWER_STDOUT_WORLD);//CHKERRQ(ierr);	
		//-------------------------------------------				
		//Applying boundary conditions on Momentum

		MatZeroRowsColumns(A, temp.size(), dof, 1.0, x, b);
		//End Applying boundary conditions on Momentum
		//-------------------------------------------
		
		//-------------------------------------------
		//Solving system.

		//Create KSP context to solve the linear system
		ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
		//CHKERRQ(ierr);
		ierr = KSPSetOperators(ksp, A, A);
		//CHKERRQ(ierr);
		//Solve using MUMPS
		#if defined(PETSC_HAVE_MUMPS)
			ierr = KSPSetType(ksp, KSPPREONLY);
			ierr = KSPGetPC(ksp, &pc);
			ierr = PCSetType(pc, PCLU);
		#endif
		ierr = KSPSetFromOptions(ksp);
		//CHKERRQ(ierr);
		ierr = KSPSetUp(ksp);

		//Solve linear system
		ierr = KSPSolve(ksp, b, x);
		//CHKERRQ(ierr);
		ierr = KSPGetTotalIterations(ksp, &iterations);

		//MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
		//VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
		//VecView(x,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

		//Gathers the solution vector to the master process
		ierr = VecScatterCreateToAll(x, &ctx, &All);
		//CHKERRQ(ierr);
		ierr = VecScatterBegin(ctx, x, All, INSERT_VALUES, SCATTER_FORWARD);
		//CHKERRQ(ierr);
		ierr = VecScatterEnd(ctx, x, All, INSERT_VALUES, SCATTER_FORWARD);
		//CHKERRQ(ierr);
		ierr = VecScatterDestroy(&ctx);
		//CHKERRQ(ierr);
		//End solving system
		//---------------------------------------------------------------
			
		//---------------------------------------------------------------
		//Updating momentum
		Ione = 1;

		for (size_t i = 0; i < nodes_.size(); i++)
		{
			bounded_vector<double,2> deltarhou;
			Idof = 2 * i;
			ierr = VecGetValues(All, Ione, &Idof, &val);
			//CHKERRQ(ierr);
			deltarhou(0) = val;

			Idof = 2 * i + 1;
			ierr = VecGetValues(All, Ione, &Idof, &val);
			//CHKERRQ(ierr);
			deltarhou(1) = val;

//	********************************************************
//	Modifications to manually introduce specific initial conditions (example 3)
//	********************************************************
			/*bounded_vector<double, 2> coord = nodes_[i]->getCurrentCoordinate(); 
			if(coord(0)>0.29999){

				if (coord(1)<0.00001){
					deltarhou(0) = 0.;
				}
			}*/
		
		

//	********************************************************
//	end Modifications to manually introduce specific initial conditions
//	********************************************************
						
		nodes_[i]->setDeltaMomentum(deltarhou); 
		}
		//End updating momentum
		//---------------------------------------------------------------
		ierr = KSPDestroy(&ksp); //CHKERRQ(ierr);
		ierr = VecDestroy(&b); //CHKERRQ(ierr);
		ierr = VecDestroy(&x); //CHKERRQ(ierr);
		ierr = VecDestroy(&All); //CHKERRQ(ierr);
		ierr = MatDestroy(&A); //CHKERRQ(ierr);
		//ierr = VecDestroy(&dispMeshValues); //CHKERRQ(ierr);

		MPI_Barrier(PETSC_COMM_WORLD);


	//****************************************
	// Small Speed Shock Capture for momentum
	//****************************************

		for (Element *el : elements_) // loop over elements
		{
			if (elementPartition_[el->getIndex()] == rank) // if element in rank
			{
				bounded_vector<double, 6> rhsS;
				//std::cout << "Element " << el->getIndex() << "\n";
				rhsS = el->elementShockLowSpeedS1();

				for (size_t i = 0; i < el->getNodes().size(); i++)
				{ 
					int idof = 2 * el->getNode(i)->getIndex();
						ierr = VecSetValues(rhsSG, 1, &idof, &rhsS(2*i), ADD_VALUES);
						idof = 2 * el->getNode(i)->getIndex() + 1;
						ierr = VecSetValues(rhsSG, 1, &idof, &rhsS(2*i+1), ADD_VALUES);
				} // i node
			} // if element in rank
		
		} // loop over elements


		ierr = VecAssemblyBegin(rhsSG); //CHKERRQ(ierr);
    	ierr = VecAssemblyEnd(rhsSG); //CHKERRQ(ierr);

		ierr = VecScatterCreateToAll(rhsSG, &ctxrhsSG, &AllrhsSG);
		//CHKERRQ(ierr);
		ierr = VecScatterBegin(ctxrhsSG, rhsSG, AllrhsSG, INSERT_VALUES, SCATTER_FORWARD);
		//CHKERRQ(ierr);
		ierr = VecScatterEnd(ctxrhsSG, rhsSG, AllrhsSG, INSERT_VALUES, SCATTER_FORWARD);
		//CHKERRQ(ierr);
		ierr = VecScatterDestroy(&ctxrhsSG);

		
		Ione = 1;
		for (size_t i = 0; i < nodes_.size(); i++)
			{
				double mlump=nodes_[i]->getLumpedMass();
				bounded_vector<double,2> drou1=nodes_[i]->getDeltaMomentum();
				bounded_vector<double,2> drouS;

				Idof = 2*i;

				ierr = VecGetValues(AllrhsSG, Ione, &Idof, &val);
				//CHKERRQ(ierr);
				double drou1s = (1./(1.+0.5*alpha_))*drou1(0)+(alpha_/(1.+0.5*alpha_))*(1./mlump)*val;

				Idof = 2*i+1;

				ierr = VecGetValues(AllrhsSG, Ione, &Idof, &val);
				//CHKERRQ(ierr);
				double drou2s = (1./(1.+0.5*alpha_))*drou1(1)+(alpha_/(1.+0.5*alpha_))*(1./mlump)*val;
				drouS(0)=drou1s;
				drouS(1)=drou2s;

				//	********************************************************
				//	Modifications to manually introduce specific initial conditions (example 4)
				//	********************************************************

				/*drouS(0)=0.;
				drouS(1)=0.;*/
					

				//	********************************************************
				//	end Modifications to manually introduce specific initial conditions
				//	********************************************************
				nodes_[i]->setDeltaMomentum(drouS);
	
			}


		ierr = VecDestroy(&rhsSG); //CHKERRQ(ierr);
		ierr = VecDestroy(&AllrhsSG); //CHKERRQ(ierr);


		MPI_Barrier(PETSC_COMM_WORLD);

	//*********************************
	//	Solves step2 (implicit or explicit)
	//*********************************
	
		//Create PETSc sparse parallel matrix
        ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                            nodes_.size() , nodes_.size(),
                            100,NULL,300,NULL,&A); 

		//	CHKERRQ(ierr);
            
         ierr = MatGetOwnershipRange(A, &Istart, &Iend);//CHKERRQ(ierr);
		             
        //Create PETSc vectors
        ierr = VecCreate(PETSC_COMM_WORLD,&b); //CHKERRQ(ierr);
        ierr = VecSetSizes(b,PETSC_DECIDE,nodes_.size()); //CHKERRQ(ierr);
        ierr = VecSetFromOptions(b); //CHKERRQ(ierr);
        ierr = VecDuplicate(b,&x); //CHKERRQ(ierr);
        ierr = VecDuplicate(b,&All); //CHKERRQ(ierr); 
		ierr = VecCreate(PETSC_COMM_WORLD,&rhsSG);
        ierr = VecSetSizes(rhsSG,PETSC_DECIDE,nodes_.size()); //CHKERRQ(ierr);
		ierr = VecSetFromOptions(rhsSG); //CHKERRQ(ierr);
		ierr = VecDuplicate(rhsSG,&AllrhsSG); //CHKERRQ(ierr);

		//std::cout<< "Passo 00" << std::endl;
		for (Element *el : elements_) // loop over elements
		{
			
			// Integração numérica das matrizes locais e montagem na matriz global
			// Integhação numérica do vetor incógnita e montagem no vetor global

			if (elementPartition_[el->getIndex()] == rank) // if element in rank
			{

				std::pair<bounded_vector<double, 3>, bounded_matrix<double, 3, 3> > elementMatrices;
				//std::cout << "Element " << el->getIndex() << "\n";
				elementMatrices = el->elementContributionsS2(qdif_, theta_, sutherland_);

				for (size_t i = 0; i < el->getNodes().size(); i++)
				{ 
					int idof = el->getNode(i)->getIndex();

					ierr = VecSetValues(b, 1, &idof, &elementMatrices.first(i), ADD_VALUES);

					for (size_t j = 0; j < el->getNodes().size(); j++)
						{
							int dof1 = el->getNode(i)->getIndex();
							int dof2 = el->getNode(j)->getIndex();
							ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(i, j), ADD_VALUES);
						}
				}
			}
				
		}
		// resolver sistema linear

	 		//Assemble matrices and vectors
    		ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); //CHKERRQ(ierr);
    		ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); //CHKERRQ(ierr);

    		ierr = VecAssemblyBegin(b); //CHKERRQ(ierr);
    		ierr = VecAssemblyEnd(b); //CHKERRQ(ierr);

//			MatView(A,PETSC_VIEWER_STDOUT_WORLD);//CHKERRQ(ierr);
//			VecView(b,PETSC_VIEWER_STDOUT_WORLD);//CHKERRQ(ierr);
//			VecView(x,PETSC_VIEWER_STDOUT_WORLD);//CHKERRQ(ierr);	
			
			//aplicando condição de contorno

			//MatZeroRowsColumns(A, temp.size(), dof, 1.0, x, b);

			//Create KSP context to solve the linear system
			ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
			//CHKERRQ(ierr);
			ierr = KSPSetOperators(ksp, A, A);
			//CHKERRQ(ierr);
			//Solve using MUMPS
			#if defined(PETSC_HAVE_MUMPS)
			ierr = KSPSetType(ksp, KSPPREONLY);
			ierr = KSPGetPC(ksp, &pc);
			ierr = PCSetType(pc, PCLU);
			#endif
			ierr = KSPSetFromOptions(ksp);
			//CHKERRQ(ierr);
			ierr = KSPSetUp(ksp);

			//Solve linear system
			ierr = KSPSolve(ksp, b, x);
			//CHKERRQ(ierr);
			ierr = KSPGetTotalIterations(ksp, &iterations);

			//MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
			//VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
			//VecView(x,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

			//Gathers the solution vector to the master process
			ierr = VecScatterCreateToAll(x, &ctx, &All);
			//CHKERRQ(ierr);
			ierr = VecScatterBegin(ctx, x, All, INSERT_VALUES, SCATTER_FORWARD);
			//CHKERRQ(ierr);
			ierr = VecScatterEnd(ctx, x, All, INSERT_VALUES, SCATTER_FORWARD);
			//CHKERRQ(ierr);
			ierr = VecScatterDestroy(&ctx);
			//CHKERRQ(ierr);

			Ione = 1;

			for (size_t i = 0; i < nodes_.size(); i++)
			{
				double deltarho;
				Idof = i;

				ierr = VecGetValues(All, Ione, &Idof, &val);
				//CHKERRQ(ierr);
				deltarho = val;

				//	********************************************************
				//	Modifications to manually introduce specific initial conditions (example 4)
				//	********************************************************

				//deltarho = 0.;
					

				//	********************************************************
				//	end Modifications to manually introduce specific initial conditions
				//	********************************************************

				nodes_[i]->setDeltaDensity(deltarho); 

				//	********************************************************
				//	Modifications to manually introduce specific initial conditions (example 3)
				//	********************************************************

				/*bounded_vector<double, 2> coord = nodes_[i]->getCurrentCoordinate();
				if(coord(0)<.000001){

					nodes_[i]->setDeltaDensity(0.);

				}*/

				//	********************************************************
				//	end Modifications to manually introduce specific initial conditions
				//	********************************************************
			}

			// calcular o erro (norma do vetor correção)
			ierr = KSPDestroy(&ksp); //CHKERRQ(ierr);
			ierr = VecDestroy(&b); //CHKERRQ(ierr);
			ierr = VecDestroy(&x); //CHKERRQ(ierr);
			ierr = VecDestroy(&All); //CHKERRQ(ierr);
			ierr = MatDestroy(&A); //CHKERRQ(ierr);
			// 	std::cout<<"Iteration "<<iterationStep2<<std::endl;

			MPI_Barrier(PETSC_COMM_WORLD);

	//****************************************
	// Small Speed Shock Capture for density
	//****************************************

		for (Element *el : elements_) // loop over elements
		{
			if (elementPartition_[el->getIndex()] == rank) // if element in rank
			{
			
				bounded_vector<double, 3> rhsS;
				//std::cout << "Element " << el->getIndex() << "\n";
				rhsS = el->elementShockLowSpeedS2();

				for (size_t i = 0; i < el->getNodes().size(); i++)
				{ 
					int idof = el->getNode(i)->getIndex();
						ierr = VecSetValues(rhsSG, 1, &idof, &rhsS(i), ADD_VALUES);
						
				} // i node
			} // if element in rank
		
		} // loop over elements


		ierr = VecAssemblyBegin(rhsSG); //CHKERRQ(ierr);
    	ierr = VecAssemblyEnd(rhsSG); //CHKERRQ(ierr);

		ierr = VecScatterCreateToAll(rhsSG, &ctxrhsSG, &AllrhsSG);
		//CHKERRQ(ierr);
		ierr = VecScatterBegin(ctxrhsSG, rhsSG, AllrhsSG, INSERT_VALUES, SCATTER_FORWARD);
		//CHKERRQ(ierr);
		ierr = VecScatterEnd(ctxrhsSG, rhsSG, AllrhsSG, INSERT_VALUES, SCATTER_FORWARD);
		//CHKERRQ(ierr);
		ierr = VecScatterDestroy(&ctxrhsSG);
		//Assembly. Scatter... rhSG

		
		Ione = 1;
		for (size_t i = 0; i < nodes_.size(); i++)
			{
				double mlump=nodes_[i]->getLumpedMass();
				double drho1=nodes_[i]->getDeltaDensity();
				double drhoS;

				Idof = i;

				ierr = VecGetValues(AllrhsSG, Ione, &Idof, &val);
				//CHKERRQ(ierr);
				drhoS = (1./(1.+0.5*alpha_))*drho1+(alpha_/(1.+0.5*alpha_))*(1./mlump)*val;

				

				nodes_[i]->setDeltaDensity(drhoS);
	
			}

		ierr = VecDestroy(&rhsSG); //CHKERRQ(ierr);
		ierr = VecDestroy(&AllrhsSG); //CHKERRQ(ierr);

		MPI_Barrier(PETSC_COMM_WORLD);

	
	//***********************************************************************
	//	Solves step3 (explicit)
	//***********************************************************************
		//Create PETSc sparse parallel matrix
        ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                            nodes_.size() , nodes_.size(),
                            100,NULL,300,NULL,&A); 
		//	CHKERRQ(ierr);
            
        ierr = MatGetOwnershipRange(A, &Istart, &Iend);//CHKERRQ(ierr);
            
        //Create PETSc vectors
        ierr = VecCreate(PETSC_COMM_WORLD,&b); //CHKERRQ(ierr);
        ierr = VecSetSizes(b,PETSC_DECIDE,nodes_.size()); //CHKERRQ(ierr);
        ierr = VecSetFromOptions(b); //CHKERRQ(ierr);
        ierr = VecDuplicate(b,&x); //CHKERRQ(ierr);
        ierr = VecDuplicate(b,&All); //CHKERRQ(ierr); 
		ierr = VecCreate(PETSC_COMM_WORLD,&rhsSG);
        ierr = VecSetSizes(rhsSG,PETSC_DECIDE,nodes_.size()); //CHKERRQ(ierr);
		ierr = VecSetFromOptions(rhsSG); //CHKERRQ(ierr);
		ierr = VecDuplicate(rhsSG,&AllrhsSG); //CHKERRQ(ierr);
		
		//std::cout<< "Passo 00" << std::endl;
		for (Element *el : elements_) // loop over elements
		{
			
			// Integração numérica das matrizes locais e montagem na matriz global
			// Integhação numérica do vetor incógnita e montagem no vetor global

			if (elementPartition_[el->getIndex()] == rank) // if element in rank
					{

						std::pair<bounded_vector<double, 3>, bounded_matrix<double, 3, 3> > elementMatrices;
						//std::cout << "Element " << el->getIndex() << "\n";
						elementMatrices = el->elementContributionsS3(qdif_, theta_, sutherland_);

						for (size_t i = 0; i < el->getNodes().size(); i++)
						{ 
							int idof = el->getNode(i)->getIndex();
							ierr = VecSetValues(b, 1, &idof, &elementMatrices.first(i), ADD_VALUES);

			
							for (size_t j = 0; j < el->getNodes().size(); j++)
							{
								int dof1 = el->getNode(i)->getIndex();
								int dof2 = el->getNode(j)->getIndex();
								ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(i, j), ADD_VALUES);



							}
							

						}
					}
				
			
		}

			// resolver sistema linear

	 		//Assemble matrices and vectors
    		ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); //CHKERRQ(ierr);
    		ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); //CHKERRQ(ierr);
	
    		ierr = VecAssemblyBegin(b); //CHKERRQ(ierr);
    		ierr = VecAssemblyEnd(b); //CHKERRQ(ierr);

			//MatView(A,PETSC_VIEWER_STDOUT_WORLD);//CHKERRQ(ierr);
			//VecView(b,PETSC_VIEWER_STDOUT_WORLD);//CHKERRQ(ierr);
			//VecView(x,PETSC_VIEWER_STDOUT_WORLD);//CHKERRQ(ierr);

			//aplicando condição de contorno

			//MatZeroRowsColumns(A, temp.size(), dof, 1.0, x, b);

			//	resolver o sistema aqui.

			//Create KSP context to solve the linear system
			ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
			//CHKERRQ(ierr);
			ierr = KSPSetOperators(ksp, A, A);
			//CHKERRQ(ierr);
			//Solve using MUMPS
			#if defined(PETSC_HAVE_MUMPS)
			ierr = KSPSetType(ksp, KSPPREONLY);
			ierr = KSPGetPC(ksp, &pc);
			ierr = PCSetType(pc, PCLU);
			#endif
			ierr = KSPSetFromOptions(ksp);
			//CHKERRQ(ierr);
			ierr = KSPSetUp(ksp);

			//Solve linear system
			ierr = KSPSolve(ksp, b, x);

			//CHKERRQ(ierr);
			ierr = KSPGetTotalIterations(ksp, &iterations);

			//MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
			//VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
			//VecView(x,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

			//Gathers the solution vector to the master process
			ierr = VecScatterCreateToAll(x, &ctx, &All);
			//CHKERRQ(ierr);
			ierr = VecScatterBegin(ctx, x, All, INSERT_VALUES, SCATTER_FORWARD);
			//CHKERRQ(ierr);
			ierr = VecScatterEnd(ctx, x, All, INSERT_VALUES, SCATTER_FORWARD);
			//CHKERRQ(ierr);
			ierr = VecScatterDestroy(&ctx);
			//CHKERRQ(ierr);



			Ione = 1;

			for (size_t i = 0; i < nodes_.size(); i++)
			{
				double deltarhoe;
				Idof = i;
				ierr = VecGetValues(All, Ione, &Idof, &val);
				//CHKERRQ(ierr);
				deltarhoe = val;

				nodes_[i]->setDeltaEnergy(deltarhoe);
				
				//	********************************************************
				//	Modifications to manually introduce specific initial conditions (example 3)
				//	********************************************************

				/*bounded_vector<double, 2> coord = nodes_[i]->getCurrentCoordinate(); 
				if(coord(0)<0.0001){

				nodes_[i]->setDeltaEnergy(0.);

				}*/

				//	********************************************************
				//	end Modifications to manually introduce specific initial conditions
				//	********************************************************

			}

			


			// calcular o erro (norma do vetor correção)
			ierr = KSPDestroy(&ksp); //CHKERRQ(ierr);
			ierr = VecDestroy(&b); //CHKERRQ(ierr);
			ierr = VecDestroy(&x); //CHKERRQ(ierr);
			ierr = VecDestroy(&All); //CHKERRQ(ierr);
			ierr = MatDestroy(&A); //CHKERRQ(ierr);
			// 	std::cout<<"Iteration "<<iterationStep2<<std::endl;

		MPI_Barrier(PETSC_COMM_WORLD);


	//****************************************
	// Small Speed Shock Capture for energy
	//****************************************

		for (Element *el : elements_) // loop over elements
		{
			if (elementPartition_[el->getIndex()] == rank) // if element in rank
			{
			
				bounded_vector<double, 3> rhsS;
				//std::cout << "Element " << el->getIndex() << "\n";
				//rhsS = el->elementShockLowSpeedS3();

				for (size_t i = 0; i < el->getNodes().size(); i++)
				{ 
					double valuesS=0.;
					int idof = el->getNode(i)->getIndex();
						ierr = VecSetValues(rhsSG, 1, &idof, &valuesS, ADD_VALUES);
						
				} // i node
			} // if element in rank
		
		} // loop over elements

		ierr = VecAssemblyBegin(rhsSG); //CHKERRQ(ierr);
    	ierr = VecAssemblyEnd(rhsSG); //CHKERRQ(ierr);

		ierr = VecScatterCreateToAll(rhsSG, &ctxrhsSG, &AllrhsSG);
		//CHKERRQ(ierr);
		ierr = VecScatterBegin(ctxrhsSG, rhsSG, AllrhsSG, INSERT_VALUES, SCATTER_FORWARD);
		//CHKERRQ(ierr);
		ierr = VecScatterEnd(ctxrhsSG, rhsSG, AllrhsSG, INSERT_VALUES, SCATTER_FORWARD);
		//CHKERRQ(ierr);
		ierr = VecScatterDestroy(&ctxrhsSG);
		//Assembly. Scatter... rhSG
		
		Ione = 1;
		for (size_t i = 0; i < nodes_.size(); i++)
			{
				double mlump=nodes_[i]->getLumpedMass();
				double drhoe1=nodes_[i]->getDeltaEnergy();
				double drhoeS;

				Idof = i;

				ierr = VecGetValues(AllrhsSG, Ione, &Idof, &val);
				//CHKERRQ(ierr);
				drhoeS = (1./(1.+0.5*alpha_))*drhoe1+(alpha_/(1.+0.5*alpha_))*(1./mlump)*val;

				

				nodes_[i]->setDeltaEnergy(drhoeS);
	
			}


		ierr = VecDestroy(&rhsSG); //CHKERRQ(ierr);
		ierr = VecDestroy(&AllrhsSG); //CHKERRQ(ierr);



		//update variables
		int matIndex = 0; //materials_.size();
		for (Node* n : nodes_)
		{
			n->updateVariables();
			
 
		}


	

		for (Node* n : nodes_)
		{			
			n->computeSecondaryVariables(materials_[matIndex] ->getSpecificHeatv(), materials_[matIndex]->getSpecificHeatp());

			//	********************************************************
			//	Modifications to manually introduce specific initial conditions (example 3)
			//	********************************************************

			/*bounded_vector<double, 2> coord = n->getCurrentCoordinate(); 
			if(coord(0)>0.29999){

				if (coord(1)<0.00001){
					n->setCurrentTemperature(0.555555556);
				}
			}*/

			//	********************************************************
			//	end Modifications to manually introduce specific initial conditions
			//	********************************************************

			//	********************************************************
			//	Modifications to manually introduce specific initial conditions (example 4)
			//	********************************************************

			/*bounded_vector<double, 2> coord = n->getCurrentCoordinate(); 
			if (coord(1)>0.99999){
				n->setCurrentTemperature(100.);
			}*/
			

			//	********************************************************
			//	end Modifications to manually introduce specific initial conditions
			//	********************************************************
		}
	
			ierr = VecDestroy(&lumpedMass); //CHKERRQ(ierr);
			ierr = VecDestroy(&xlump); //CHKERRQ(ierr);
			ierr = VecDestroy(&AllLump); //CHKERRQ(ierr);
			

		 //std::cout<<"TO enter mesh"<<std::endl; 
		//	************************************
		//	Solves Mesh Movement
		//	************************************
		//Create PETSc sparse parallel matrix and vectors
			
        ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                            2 * nodes_.size() , 2 * nodes_.size(),
                            100,NULL,300,NULL,&A); 
		//	CHKERRQ(ierr);
    
            ierr = MatGetOwnershipRange(A, &Istart, &Iend);//CHKERRQ(ierr);
			
            //Create PETSc vectors
          	ierr = VecCreate(PETSC_COMM_WORLD,&b); //CHKERRQ(ierr);
            ierr = VecSetSizes(b,PETSC_DECIDE,2*nodes_.size()); //CHKERRQ(ierr);
            ierr = VecSetFromOptions(b); //CHKERRQ(ierr);
            ierr = VecDuplicate(b,&x); //CHKERRQ(ierr);
            ierr = VecDuplicate(b,&All); //CHKERRQ(ierr);
            ierr = VecDuplicate(b,&dispMeshValues); //CHKERRQ(ierr);

		//for (size_t i = 0; i < tempMeshValues.size(); i++)
		//{
		//	dispMeshValues[i] = tempMeshValues[i];

		//}

//std::cout<<"TO enter Elements mesh"<<std::endl; 
	
		for (Element *el : elements_) // loop over elements
		{
			if (elementPartition_[el->getIndex()] == rank) // if element in rank
			{

				bounded_matrix<double, 3, 3> elementMatrix;
						
				elementMatrix = el->elementContributionsMesh();

			 	for (size_t i = 0; i < el->getNodes().size(); i++)
			 	{ 
			 		for (size_t j = 0; j < el->getNodes().size(); j++)
			 		{
			 			int dof1 = 2 * el->getNode(i)->getIndex();
			 			int dof2 = 2 * el->getNode(j)->getIndex();
			 			ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrix(i, j), ADD_VALUES);

			 			dof1 = 2 * el->getNode(i)->getIndex()+1;
			 			dof2 = 2 * el->getNode(j)->getIndex()+1;
			 			ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrix(i, j), ADD_VALUES); 

			 		} // j node
			 	} // i node
			} // if element in rank
		
		} // loop over elements
	
		
	 	//Assemble matrices and vectors
    	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); //CHKERRQ(ierr);
    	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); //CHKERRQ(ierr);


		
			
		//MatView(A,PETSC_VIEWER_STDOUT_WORLD);//CHKERRQ(ierr);
		//VecView(b,PETSC_VIEWER_STDOUT_WORLD);//CHKERRQ(ierr);
		//VecView(x,PETSC_VIEWER_STDOUT_WORLD);//CHKERRQ(ierr);	
		//-------------------------------------------				
		//Applying boundary conditions on Mesh

		ierr = MatZeroRows(A, tempMesh.size(), dofMesh, 1.0,NULL,NULL);

		if(rank == 0 ){
		for (size_t i = 0; i < tempMeshValues.size(); i++)
		{
		int index = tempMesh[i];
								//std::cout<<"mlump "<<mlump<<std::endl;
					ierr = VecSetValues(b, 1, &index, &tempMeshValues[i], INSERT_VALUES);
		}
		}
    	ierr = VecAssemblyBegin(b); //CHKERRQ(ierr);
    	ierr = VecAssemblyEnd(b); //CHKERRQ(ierr);
	 	

		//End Applying boundary conditions on Momentum
		//-------------------------------------------
		
		//-------------------------------------------
		//Solving system.

		//Create KSP context to solve the linear system
		ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
		//CHKERRQ(ierr);
		ierr = KSPSetOperators(ksp, A, A);
		//CHKERRQ(ierr);
		//Solve using MUMPS
		#if defined(PETSC_HAVE_MUMPS)
			ierr = KSPSetType(ksp, KSPPREONLY);
			ierr = KSPGetPC(ksp, &pc);
			ierr = PCSetType(pc, PCLU);
		#endif
		ierr = KSPSetFromOptions(ksp);
		//CHKERRQ(ierr);
		ierr = KSPSetUp(ksp);

		//Solve linear system
		ierr = KSPSolve(ksp, b, x);
		//CHKERRQ(ierr);
		ierr = KSPGetTotalIterations(ksp, &iterations);

		//MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
		//VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
		//VecView(x,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

		//Gathers the solution vector to the master process
		ierr = VecScatterCreateToAll(x, &ctx, &All);
		//CHKERRQ(ierr);
		ierr = VecScatterBegin(ctx, x, All, INSERT_VALUES, SCATTER_FORWARD);
		//CHKERRQ(ierr);
		ierr = VecScatterEnd(ctx, x, All, INSERT_VALUES, SCATTER_FORWARD);
		//CHKERRQ(ierr);
		ierr = VecScatterDestroy(&ctx);
		//CHKERRQ(ierr);
		//End solving system
		//---------------------------------------------------------------
			
		//---------------------------------------------------------------
		//Updating mesh
		Ione = 1;

		for (size_t i = 0; i < nodes_.size(); i++)
		{
		 	bounded_vector<double,2> deltaMesh;
		 	Idof = 2 * i;
		 	ierr = VecGetValues(All, Ione, &Idof, &val);
		 	//CHKERRQ(ierr);
		 	deltaMesh(0) = val;
			//std::cout<<i<<" dxmesh "<<val<<std::endl;
		 	Idof = 2 * i + 1;
		 	ierr = VecGetValues(All, Ione, &Idof, &val);
		 	//CHKERRQ(ierr);
		 	deltaMesh(1) = val;	
			//std::cout<<i<<" dymesh "<<val<<std::endl;

			nodes_[i]->setDeltaMesh(deltaMesh);			
			
		}

		/*for (Node* n : nodes_)
		{
			n->updateMesh();
			
 
		}
		*/
		// //End updating momentum
		// //---------------------------------------------------------------
		// ierr = KSPDestroy(&ksp); //CHKERRQ(ierr);
		
		ierr = VecDestroy(&b); //CHKERRQ(ierr);
		ierr = VecDestroy(&x); //CHKERRQ(ierr);
		ierr = VecDestroy(&All); //CHKERRQ(ierr);
		ierr = MatDestroy(&A); //CHKERRQ(ierr);

		MPI_Barrier(PETSC_COMM_WORLD);
	exportToParaview(itimeStep+1);
	

}

}



void FluidDomain::exportToParaview(const int& timestep)
{
	std::stringstream text;
	text << "results/" << "FluidOutput" << timestep << ".vtu";
	std::ofstream file(text.str());
	file.precision(16);

	//header
	file << "<?xml version=\"1.0\"?>" << "\n"
         << "<VTKFile type=\"UnstructuredGrid\">" << "\n"
		 << "  <UnstructuredGrid>" << "\n"
         << "  <Piece NumberOfPoints=\"" << nodes_.size()
         << "\"  NumberOfCells=\"" << elements_.size()
         << "\">" << "\n";
	//nodal coordinates
	file << "    <Points>" << "\n"
         << "      <DataArray type=\"Float64\" "
         << "NumberOfComponents=\"3\" format=\"ascii\">" << "\n";
    for (Node* n : nodes_)
	{
		file << n->getCurrentCoordinate()(0) << " " << n->getCurrentCoordinate()(1) << " " << 0.0 << "\n";
	}
    file << "      </DataArray>" << "\n"
         << "    </Points>" << "\n";
	//element connectivity
	file << "    <Cells>" << "\n"
         << "      <DataArray type=\"Int32\" "
         << "Name=\"connectivity\" format=\"ascii\">" << "\n";
    for (Element* e : elements_)
	{ 
		for (Node* n : e->getNodes())
		{
			file << n->getIndex() << " ";
		}
		file << "\n";
	}

	file << "      </DataArray>" << "\n";
	//offsets
	file << "      <DataArray type=\"Int32\""
		 << " Name=\"offsets\" format=\"ascii\">" << "\n";
	int aux = 0;
	for (Element* e : elements_)
	{
		aux += 3;
		file << aux << "\n";
	}
	file << "      </DataArray>" << "\n";
	//elements type
	file << "      <DataArray type=\"UInt8\" Name=\"types\" "
		 << "format=\"ascii\">" << "\n";

	for (Element* e : elements_)
	{
		file << 5 << "\n";
	}
	file << "      </DataArray>" << "\n"
		 << "    </Cells>" << "\n";
	//nodal results
	file << "    <PointData>" <<"\n";
	file << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
		<< "Name=\"Velocity\" format=\"ascii\">" << "\n";

	for (Node* n: nodes_)
	{
	 bounded_vector<double, 2> u_ = n->getCurrentVelocity();
		// std::cout << dx << " "
		// 	 << dy  << std::endl;
		file << u_(0) << " "
			 << u_(1)  << " " << 0.0 << "\n";
	}
	file << "      </DataArray> " << "\n";
	file << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
		<< "Name=\"MeshDisp\" format=\"ascii\">" << "\n";

	for (Node* n: nodes_)
	{
	 bounded_vector<double, 2> u_ = n->getDeltaMesh();
		// std::cout << dx << " "
		// 	 << dy  << std::endl;
		file << u_(0) << " "
			 << u_(1)  << " " << 0.0 << "\n";
	}
	file << "      </DataArray> " << "\n";

	file << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
		<< "Name=\"Momentum\" format=\"ascii\">" << "\n";

	for (Node* n: nodes_)
	{
	 bounded_vector<double, 2> u_ = n->getCurrentMomentum();
		// std::cout << dx << " "
		// 	 << dy  << std::endl;
		file << u_(0) << " "
			 << u_(1)  << " " << 0.0 << "\n";
	}
	file << "      </DataArray> " << "\n";
	file << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
		<< "Name=\"Pressure\" format=\"ascii\">" << "\n";

	for (Node* n: nodes_)
	{
	 double p = n->getCurrentPressure();
		// std::cout << dx << " "
		// 	 << dy  << std::endl;
		file << p << " "
			 << 0.0  << " " << 0.0 << "\n";
	}
	file << "      </DataArray> " << "\n";
	file << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
		<< "Name=\"Energy\" format=\"ascii\">" << "\n";

	for (Node* n: nodes_)
	{
	 double p = n->getCurrentInternalEnergy();
		// std::cout << dx << " "
		// 	 << dy  << std::endl;
		file << p << " "
			 << 0.0  << " " << 0.0 << "\n";
	}
	file << "      </DataArray> " << "\n";
		file << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
		<< "Name=\"Density\" format=\"ascii\">" << "\n";

	for (Node* n: nodes_)
	{
	 double p = n->getCurrentDensity();
		// std::cout << dx << " "
		// 	 << dy  << std::endl;
		file << p << " "
			 << 0.0  << " " << 0.0 << "\n";
	}
	file << "      </DataArray> " << "\n";

	file << "      <DataArray type=\"Float64\" NumberOfComponents=\"2\" "
		<< "Name=\"Temperature\" format=\"ascii\">" << "\n";

	for (Node* n: nodes_)
	{
	 double p = n->getCurrentTemperature();
		// std::cout << dx << " "
		// 	 << dy  << std::endl;
		file << p << " "
			 << 0.0  << "\n";
	}
	file << "      </DataArray> " << "\n";
	file << "      <DataArray type=\"Float64\" NumberOfComponents=\"2\" "
		<< "Name=\"LumpedMass\" format=\"ascii\">" << "\n";

	for (Node* n: nodes_)
	{
	 double p = n->getLumpedMass();
		// std::cout << dx << " "
		// 	 << dy  << std::endl;
		file << p << " "
			 << 0.0  << "\n";
	}
	file << "      </DataArray> " << "\n";
	file << "    </PointData>" << "\n";
	//elemental results
	file << "    <CellData>" << "\n";

	file << "    </CellData>" << "\n";
	//footnote
	file << "  </Piece>" << "\n"
		<< "  </UnstructuredGrid>" << "\n"
		<< "</VTKFile>" << "\n";
	file.close();

}

void FluidDomain::domainDecompositionMETIS() {
    
    std::string mirror2;
    mirror2 = "domain_decomposition.txt";
    std::ofstream mirrorData(mirror2.c_str());
    
    int size;

    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    idx_t objval;
    idx_t numEl = elements_.size();
    idx_t numNd = nodes_.size();
    idx_t ssize = size;
    idx_t one = 1;
    idx_t elem_start[numEl+1], elem_connec[3*numEl];
    elementPartition_ = new idx_t[numEl];
    nodePartition_ = new idx_t[numNd];


    for (idx_t i = 0; i < numEl+1; i++)
	{
        elem_start[i]=3*i;
    }
    for (idx_t jel = 0; jel < numEl; jel++)
	{   
        for (idx_t i=0; i < elements_[jel]->getNodes().size(); i++)
		{	
			int nodeIndex = elements_[jel]->getNodes()[i]->getIndex();
        	elem_connec[3*jel+i] = nodeIndex;
        }
    }

    //Performs the domain decomposition
	if (size == 1)
	{
		for (Node *node : nodes_)
			nodePartition_[node->getIndex()] = 0;
		for (Element *el : elements_)
			elementPartition_[el->getIndex()] = 0;
	}
	else
	{
		METIS_PartMeshDual(&numEl, &numNd, elem_start, elem_connec,
						NULL, NULL, &one, &ssize, NULL, NULL,
						&objval, elementPartition_, nodePartition_);
	}

    mirrorData << std::endl \
               << "FLUID MESH DOMAIN DECOMPOSITION - ELEMENTS" << std::endl;
    for (int i = 0; i < elements_.size(); i++)
	{
        mirrorData << "process = " << elementPartition_[i] \
                   << ", element = " << i << std::endl;
    }

    mirrorData << std::endl \
               << "FLUID MESH DOMAIN DECOMPOSITION - NODES" << std::endl;
    for (int i = 0; i < nodes_.size(); i++)
	{
        mirrorData << "process = " << nodePartition_[i] \
                   << ", node = " << i << std::endl;
    }

}
