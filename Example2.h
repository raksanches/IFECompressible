//generating fluidMesh
Geometry* fluid = new Geometry(0);

//Create base points
	Point* p1 = fluid->addPoint({ 0., 0. }, 0.01);
	Point* p2 = fluid->addPoint({ .6, 0. }, 0.01);
	Point* p3 = fluid->addPoint({ .6, .2 }, 0.01);
	Point* p4 = fluid->addPoint({ 3., .2 }, 0.01);
	Point* p5 = fluid->addPoint({ 3., 1. }, 0.01);
	Point* p6 = fluid->addPoint({ 0., 1. }, 0.01);

//Create base lines;
  	Line* l1 = fluid->addLine({ p1, p2});
  	Line* l2 = fluid->addLine({ p2, p3});
	Line* l3 = fluid->addLine({ p3, p4});
	Line* l4 = fluid->addLine({ p4, p5});
	Line* l5 = fluid->addLine({ p5, p6});
	Line* l6 = fluid->addLine({ p6, p1});

//Create surfaces. parameter: vector of lines
	PlaneSurface* s1 = fluid->addPlaneSurface({ l1, l2, l3, l4, l5, l6});

//setting number of line division for mesh creation: line, ndiv
	fluid->transfiniteLine({ l1 }, 12);
	fluid->transfiniteLine({ l2 }, 6); 
	fluid->transfiniteLine({ l3 }, 48);
	fluid->transfiniteLine({ l4 }, 16);
	fluid->transfiniteLine({ l5 }, 60);
	fluid->transfiniteLine({ l6 }, 20);

//adding boundary conditions (fluid dirichlet, fluid neumann, mesh displacement)
//  	fluid->addBoundaryCondition("NEUMANN", l1, {0.0}, {0.0}, "GLOBAL");
// deixar: tipo, linha, restritox?, restritoy?, valorx, valory
	fluid->addBoundaryCondition("DIRICHLET", l2, true, false, 0., 0.0, "GLOBAL");
	fluid->addBoundaryCondition("DIRICHLET", l3, false, true, 0.0, 0., "GLOBAL");
	fluid->addBoundaryCondition("NEUMANN", l4, true, true, 0., 0.0, "GLOBAL");
	fluid->addBoundaryCondition("DIRICHLET", l1, false, true, 0.0, 0., "GLOBAL");
	fluid->addBoundaryCondition("DIRICHLET", l5, false, true, 0.0, 0., "GLOBAL");
	fluid->addBoundaryCondition("DIRICHLET", l6, true, true, 1.0, 0., "GLOBAL");
	//fluid->addBoundaryCondition("NEUMANN", l2, {0.0}, {0.0}, "GLOBAL");

//values: 0 for fixed and 1 for fluid/structure interaction
	fluid->addBoundaryCondition("MESHDISPLACEMENT", l2, true, true, 0., 0.0, "GLOBAL");
	fluid->addBoundaryCondition("MESHDISPLACEMENT", l3, false, true, 0.0, 0., "GLOBAL");
	fluid->addBoundaryCondition("MESHDISPLACEMENT", l4, true, true, 1., 0.0, "GLOBAL");
	fluid->addBoundaryCondition("MESHDISPLACEMENT", l1, true, true, 0.0, 0., "GLOBAL");
	fluid->addBoundaryCondition("MESHDISPLACEMENT", l5, false, true, 0.0, 0., "GLOBAL");
	fluid->addBoundaryCondition("MESHDISPLACEMENT", l6, true, true, 0.0, 0., "GLOBAL");


	//creating fluid domain. parameter: fluid geometry
	FluidDomain* problem = new FluidDomain(fluid);
	
	//adding a material to the FLUID surface.
  	// undisturbed parameters: surface, viscosity, density,  Thermal condutivity, Isabaric specific heat, isochoric specific heat, undisturbed Temperature 

	problem->addSurfaceMaterial({ s1 }, 0., 1.0, 0., 1.4, 1., 1.785714286);
	

	//setting FLUID analysis parameters. parameter: number of time steps, time step, gravityforce, time integrator parameter theta)
	problem->setAnalysisParameters(40000, 0.001, 0.0, .5);
	
	bounded_vector<double,2>uinf;

	uinf(0)=3.0;
	uinf(1)=0.0;
	problem->setInitialVelocities(uinf);

	bounded_vector<double,2>uMesh;

	uMesh(0)=0.0;
	uMesh(1)=0.0;
	problem->setInitialMeshVelocities(uMesh);

	//setting shock capturing parameter

	problem->setShockCapture(0.5);
	
	problem->setSmallVelocityShockCapture(0.);

	problem->useSutherland(false);
	
	// FLUID generating mesh. parameter: element type, method, file .geo name, gmsh path, show mesh, mesh information
	problem->generateMesh("T3", "FRONT", "fluid", "", false, false);
	//solving transient problem:
	problem->solveCompressibleFlowMoving();
