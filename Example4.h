//generating fluidMesh
Geometry* fluid = new Geometry(0);

//Create base points
	Point* p1 = fluid->addPoint({ 0., 0. }, 0.01);
	Point* p2 = fluid->addPoint({ 1., 0. }, 0.01);
	Point* p3 = fluid->addPoint({ 1., 1. }, 0.01);
	Point* p4 = fluid->addPoint({ 0., 1. }, 0.01);

//Create base lines;
  	Line* l1 = fluid->addLine({ p1, p2});
  	Line* l2 = fluid->addLine({ p2, p3});
	Line* l3 = fluid->addLine({ p3, p4});
	Line* l4 = fluid->addLine({ p4, p1});

//Create surfaces. parameter: vector of lines
	PlaneSurface* s1 = fluid->addPlaneSurface({ l1, l2, l3, l4});

//setting number of line division for mesh creation: line, ndiv
	fluid->transfiniteLine({ l1 }, 20);
	fluid->transfiniteLine({ l2 }, 20); 
	fluid->transfiniteLine({ l3 }, 20);
	fluid->transfiniteLine({ l4 }, 20);
	

//adding boundary conditions
//  	fluid->addBoundaryCondition("NEUMANN", l1, {0.0}, {0.0}, "GLOBAL");
	fluid->addBoundaryCondition("DIRICHLET", l1, true, true, {0.}, {0.0}, "GLOBAL");
	fluid->addBoundaryCondition("DIRICHLET", l2, true, true, {0.0}, {0.}, "GLOBAL");
	fluid->addBoundaryCondition("DIRICHLET", l3, true, true, {0.0}, {0.}, "GLOBAL");
	fluid->addBoundaryCondition("DIRICHLET", l4, true, true, {0.0}, {0.}, "GLOBAL");
	//fluid->addBoundaryCondition("NEUMANN", l2, {0.0}, {0.0}, "GLOBAL");

	//creating fluid domain. parameter: fluid geometry
	FluidDomain* problem = new FluidDomain(fluid);
	
	//adding a material to the FLUID surface.
  	// undisturbed parameters: surface, viscosity, density,  Thermal condutivity, Isabaric specific heat, isochoric specific heat, undisturbed Temperature 

	problem->addSurfaceMaterial({ s1 }, .001, 1.0, 0., 1.4, 1., 0.);
	

	//setting FLUID analysis parameters. parameter: number of time steps, time step, gravityforce, time integrator parameter theta)
	problem->setAnalysisParameters(20000, 0.001, 0.0, 0.5);
	
	bounded_vector<double,2>uinf;

	uinf(0)=0.0;
	uinf(1)=0.0;
	problem->setInitialVelocities(uinf);

	//setting shock capturing parameter

	problem->setShockCapture(1.);
	
	problem->setSmallVelocityShockCapture(.0);

	problem->useSutherland(false);
	
	// FLUID generating mesh. parameter: element type, method, file .geo name, gmsh path, show mesh, mesh information
	problem->generateMesh("T3", "FRONT", "fluid", "", false, false);
	//solving transient problem:
	problem->solveCompressibleFlow();
