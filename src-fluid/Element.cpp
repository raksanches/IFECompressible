#include "Element.h"
#include "../libs/tensorlib/tensor4.h"
#include "../libs/tensorlib/operators.h"

Element::Element(const int& index, const std::vector<Node*>& nodes, Material* material, const double& thickness, const std::string& elementType)
{
	index_ = index;
	nodes_ = nodes;
	material_ = material;
	thickness_ = thickness;
	elementType_ = elementType;
	numberOfIntegrationPoints_ = 7;
	domainIntegrationPoints_ = hammerQuadrature(numberOfIntegrationPoints_);
	order_ = 1;
	meshStiffness_ = 1.;
	//neummanSides_ = setNeummanSides();


//	std::pair<vector<double>, vector<double> > pontosGauss = gaussLegendre(numberOfIntegrationPoints_);
	//pontosGauss_ = pontosGauss;
}

Element::~Element() {}

int Element::getIndex()
{
	return index_;
}

Node* Element::getNode(const int& index)
{
	return nodes_[index];
}

std::vector<Node*> Element::getNodes()
{
	return nodes_;
}

std::string Element::getType()
{
	return elementType_;
}

bounded_vector<double,6> Element::getRhsMomentum()
{
	bounded_vector<double,6> rhsMomentum;
	rhsMomentum[0] = nodes_[0]->getCurrentMomentum()(0);
	rhsMomentum[1] = nodes_[0]->getCurrentMomentum()(1);
	rhsMomentum[2] = nodes_[1]->getCurrentMomentum()(0);
	rhsMomentum[3] = nodes_[1]->getCurrentMomentum()(1);
	rhsMomentum[4] = nodes_[2]->getCurrentMomentum()(0);
	rhsMomentum[5] = nodes_[2]->getCurrentMomentum()(1);

	return rhsMomentum;

}

bounded_vector<double,3> Element::getRhsDensity()
{
	bounded_vector<double,3> rhsDensity;
	rhsDensity[0] = nodes_[0]->getCurrentDensity();
	rhsDensity[1] = nodes_[1]->getCurrentDensity();
	rhsDensity[2] = nodes_[2]->getCurrentDensity();

	return rhsDensity;

}

bounded_vector<double,3> Element::getRhsEnergy()
{
	bounded_vector<double,3> rhsEnergy;
	rhsEnergy[0] = nodes_[0]->getCurrentInternalEnergy();
	rhsEnergy[1] = nodes_[1]->getCurrentInternalEnergy();
	rhsEnergy[2] = nodes_[2]->getCurrentInternalEnergy();

	return rhsEnergy;

}

std::vector<Node*> Element::getSideNodes(const int& side)
{
	std::vector<Node*> nodes; nodes.reserve(2);
	if (side == 0)
	{
		nodes.push_back(nodes_[0]); nodes.push_back(nodes_[1]);
	}
	else if (side == 1)
	{
		nodes.push_back(nodes_[1]); nodes.push_back(nodes_[2]);
	}
	else if (side == 2)
	{
		nodes.push_back(nodes_[2]); nodes.push_back(nodes_[0]);
	}
	return nodes;
}



/*	std::vector<int> Element::setNeummanSides(){
		std::vector<int> bdsides; bdsides.reserve(3);
			if((nodes_[0]->getNeummanConstrain()==true)&&(nodes_[1]->getNeummanConstrain()==true))bdsides.push_back(0);
			if((nodes_[1]->getNeummanConstrain()==true)&&(nodes_[2]->getNeummanConstrain()==true))bdsides.push_back(1);
			if((nodes_[2]->getNeummanConstrain()==true)&&(nodes_[0]->getNeummanConstrain()==true))bdsides.push_back(3);
			return bdsides;
	}
*/

void Element::setAnalysisParameters(const double& deltat, const double& gravity, const double& beta, const double& gamma)
{
	deltat_ = deltat;
	gravity_ = gravity;
	beta_ = beta;
	gamma_ = gamma;
}


void Element::setTimeStep (double timeStep)
{
	timeStep_ = timeStep;
}

double Element::getTimeStep()
{
	return timeStep_;
}

std::pair<vector<double>, vector<double> > Element::gaussLegendre(const int& nga)
{
	vector<double> qsi(nga, 0.0);
	vector<double> weight(nga, 0.0);
	const double pi = std::acos(-1.0);
	double xmga, xlga, zga, p1ga, p2ga, p3ga, ppga, z1ga;
	int mga = (nga + 1) / 2;
	xmga = 0.0;
	xlga = 1.0;

	for (int iga = 1; iga <= mga; iga++)
	{
		zga = std::cos(pi*(double(iga) - 0.25) / (double(nga) + 0.5));
	g1:
		p1ga = 1.0;
		p2ga = 0.0;
		for (int jga = 1; jga <= nga; jga++)
		{
			p3ga = p2ga;
			p2ga = p1ga;
			p1ga = ((2.0 * double(jga) - 1.0) * zga * p2ga - (double(jga) - 1.0) * p3ga) / (double(jga));
		}

		ppga = nga * (zga * p1ga - p2ga) / (zga * zga - 1.0);
		z1ga = zga;
		zga = z1ga - p1ga / ppga;

		if (fabs(zga - z1ga) > 1.0e-15) goto g1;

		qsi(iga - 1) = xmga - xlga * zga;
		qsi(nga - iga) = xmga + xlga * zga;
		weight(iga - 1) = 2.0*xlga / ((1.0 - zga * zga)*ppga*ppga);
		weight(nga - iga) = weight(iga - 1);
	}

	return std::make_pair(qsi, weight);
}

matrix<double> Element:: hammerQuadrature(const int& nh)
{
	matrix<double> hammer(nh, 3, 0.0);

	if (nh == 1)
	{
		hammer(0,0) = 1.0 / 3.0; //xsi1
		hammer(0,1) = 1.0 / 3.0; //xsi2
		hammer(0,2) = 1.0 / 2.0; //weight
	}
	else if (nh == 3)
	{
		hammer(0,0) = 1.0 / 2.0;
		hammer(0,1) = 1.0 / 2.0;
		hammer(0,2) = 1.0 / 6.0;

		hammer(1,0) = 1.0 / 2.0;
		hammer(1,1) = 0.0;
		hammer(1,2) = 1.0 / 6.0;

		hammer(2,0) = 0.0;
		hammer(2,1) = 1.0 / 2.0;
		hammer(2,2) = 1.0 / 6.0;
	}
	else if (nh == 7)
	{
		hammer(0,0) = 1. / 3.;
        hammer(0,1) = 1. / 3.;
		hammer(0,2) = 0.11250;
            
        hammer(1,0) = (9. + 2. * sqrt(15.)) / 21.;
        hammer(1,1) = (6. - sqrt(15.)) / 21.;
		hammer(1,2) = (155. - sqrt(15.)) / 2400.;
          
        hammer(2,0) = (6. - sqrt(15.)) / 21.;
        hammer(2,1) = (9. + 2. * sqrt(15.)) / 21.;
		hammer(2,2) = (155. - sqrt(15.)) / 2400.;
          
        hammer(3,0) = (6. - sqrt(15.)) / 21.;
        hammer(3,1) = (6. - sqrt(15.)) / 21.;
		hammer(3,2) = (155. - sqrt(15.)) / 2400.;
          
        hammer(4,0) = (6. + sqrt(15.)) / 21.;
        hammer(4,1) = (6. + sqrt(15.)) / 21.;
		hammer(4,2) = (155. + sqrt(15.)) / 2400.;
          
        hammer(5,0) = (9. - 2. * sqrt(15.)) / 21.;
        hammer(5,1) = (6. + sqrt(15.)) / 21.;
		hammer(5,2) = (155. + sqrt(15.)) / 2400.;
          
        hammer(6,0) = (6. + sqrt(15.)) / 21.;
        hammer(6,1) = (9. - 2. * sqrt(15.)) / 21.;
		hammer(6,2) = (155. + sqrt(15.)) / 2400.;
	}
	return hammer;
}

vector<double> Element::domainShapeFunction(const double& xsi1, const double& xsi2)
{
	vector<double> phi(order_ + 2, 0.0);

	if (order_ == 1)
	{
		phi(0) = xsi1;
		phi(1) = xsi2;
		phi(2) = 1.0 - xsi1 - xsi2;
	}

	return phi;
}

matrix<double> Element::domainDerivativeShapeFunction(const double& xsi1, const double& xsi2)
{
	matrix<double> dphi_dxsi(2, order_ + 2, 0.0);

	if (order_ == 1)
	{
		dphi_dxsi(0,0) = 1.0;
		dphi_dxsi(0,1) = 0.0;
		dphi_dxsi(0,2) = -1.0;

		dphi_dxsi(1,0) = 0.0;
		dphi_dxsi(1,1) = 1.0;
		dphi_dxsi(1,2) = -1.0;
	}

	return dphi_dxsi;

}

bounded_vector<double, 2> Element::normalVector(const double& xsi)
{
	vector<double> dphi_dxsi = boundaryDerivativeShapeFunction(xsi);
	bounded_vector<double, 2> tangent;
	bounded_vector<double, 2> normal;

	for (size_t i = 0; i < 2; i++)
	{
		tangent(i) = 0.0;
		normal(i) = 0.0;
	}

	for (size_t j = 0; j < order_ + 1; j++)
	{
		tangent(0) += dphi_dxsi(j) * nodes_[j]->getCurrentCoordinate()(0);
		tangent(1) += dphi_dxsi(j) * nodes_[j]->getCurrentCoordinate()(1);
	}

	double jacobian = sqrt(pow(tangent(0), 2) + pow(tangent(1), 2));

	normal(0) = tangent(1) / jacobian;
	normal(1) = -tangent(0) / jacobian;

	return normal;
}

vector<double> Element::boundaryShapeFunction(const double& xsi)
{
	vector<double> phi(order_ + 1, 0.0);
	vector<double> nodesAdimensional(order_ + 1, 0.0);

	for (size_t i = 0; i < order_ + 1; i++)
	{
		nodesAdimensional(i) = ((2.0 * static_cast<double>(i + 1) - 2.0) / (static_cast<double>(order_)) - 1.0);
	}

	for (int j = 0; j < order_ + 1; j++)
	{
		phi(j) = 1.0;
		for (int i = 0; i < order_ + 1; i++)
		{
			double aux = 1.0;
				if (i != j)
				{
					phi(j) = phi(j) * (xsi - nodesAdimensional(i)) / (nodesAdimensional(j) - nodesAdimensional(i));
					for (int k = 0; k < order_ + 1; k++)
					{
						if ((i != k) && (j != k))
							aux *= (xsi - nodesAdimensional(k));
					}
				}
		}
	}

	return phi;
}

vector<double> Element::boundaryDerivativeShapeFunction(const double& xsi)
{
	vector<double> dphi_dxsi(order_ + 1, 0.0);
	vector<double> nodesAdimensional(order_ + 1, 0.0);

	for (size_t i = 0; i < order_ + 1; i++)
	{
		nodesAdimensional(i) = ((2.0 * static_cast<double>(i + 1) - 2.0) / (static_cast<double>(order_)) - 1.0);
	}

	for (int j = 0; j < order_ + 1; j++)
	{
		for (int i = 0; i < order_ + 1; i++)
		{
			double aux = 1.0;
			if (i != j)
			{
				for (int k = 0; k < order_ + 1; k++)
				{
					if ((i != k) && (j != k))
						aux *= (xsi - nodesAdimensional(k));
				}
				dphi_dxsi(j) += aux;
			}
		}
	}

	for (int i = 0; i < order_ + 1; i++)
	{
		for (int j = 0; j < order_ + 1; j++)
		{
			if (i != j)
				dphi_dxsi(i) /= (nodesAdimensional(i) - nodesAdimensional(j));
		}
	}

	return dphi_dxsi;
}

bounded_matrix<double, 2, 2> Element::currentJacobianMatrix(const double& xsi1, const double& xsi2)
{
	matrix<double> dphi_dxsi = domainDerivativeShapeFunction(xsi1, xsi2);
	double dx1_dxsi1 = 0.0, dx1_dxsi2 = 0.0, dx2_dxsi1 = 0.0, dx2_dxsi2 = 0.0;

	for (size_t i = 0; i < nodes_.size(); i++)
	{
		bounded_vector<double, 2> currentCoord = nodes_[i]->getCurrentCoordinate();
		dx1_dxsi1 += currentCoord(0) * dphi_dxsi(0,i);
		dx1_dxsi2 += currentCoord(0) * dphi_dxsi(1,i);
		dx2_dxsi1 += currentCoord(1) * dphi_dxsi(0,i);
		dx2_dxsi2 += currentCoord(1) * dphi_dxsi(1,i);
	}

	bounded_matrix<double, 2, 2> currentJacobianMatrix;
	currentJacobianMatrix(0,0) = dx1_dxsi1;
	currentJacobianMatrix(1,0) = dx2_dxsi1;
	currentJacobianMatrix(0,1) = dx1_dxsi2;
	currentJacobianMatrix(1,1) = dx2_dxsi2;

	return currentJacobianMatrix;
}

double Element::getArea()
{
	double area = 0.;
	bounded_vector<double, 2> currentCoord1 = nodes_[0]->getCurrentCoordinate();
	bounded_vector<double, 2> currentCoord2 = nodes_[1]->getCurrentCoordinate();
	bounded_vector<double, 2> currentCoord3 = nodes_[2]->getCurrentCoordinate();
	double D = 0.;

	D = fabs(currentCoord1(0)*currentCoord2(1)+currentCoord2(0)*currentCoord3(1)+currentCoord3(0)*currentCoord1(1)-(currentCoord1(1)*currentCoord2(0)+currentCoord2(1)*currentCoord3(0)+currentCoord3(1)*currentCoord1(0)));
	area = D/2.;

	return area;

}

double Element:: getH()
{
	double area = getArea();
	double h = 0.;
	double b =0.;
	double b1 =0.;
	double b2 =0.;
	double b3 =0.;
	bounded_vector<double, 2> currentCoord1 = nodes_[0]->getCurrentCoordinate();
	bounded_vector<double, 2> currentCoord2 = nodes_[1]->getCurrentCoordinate();
	bounded_vector<double, 2> currentCoord3 = nodes_[2]->getCurrentCoordinate();
	
	b1 = sqrt(pow(currentCoord2(0)-currentCoord1(0),2)+pow(currentCoord2(1)-currentCoord1(1),2));
	b2 = sqrt(pow(currentCoord3(0)-currentCoord2(0),2)+pow(currentCoord3(1)-currentCoord2(1),2));
	b3 = sqrt(pow(currentCoord1(0)-currentCoord3(0),2)+pow(currentCoord1(1)-currentCoord3(1),2));

	if (b1>b2 && b1>b3){
		b = b1;
	}
	else if (b2>b1 && b2>b3){
		b = b2;
	}
	else if (b3>b1 && b3>b2){
		b = b3;
	}
	h = 2.*area/b;

	return h;
}

std::pair<bounded_vector<double, 6>, bounded_matrix<double, 6, 6> > Element::elementContributionsS1(const double& qdif, const double& teta_, const bool& sutherland)
{	

	bounded_vector<double, 6> rhs;
	bounded_matrix<double, 6, 6> tangent;     // 2 por nó
	rhs.clear(); tangent.clear();
	double jac=0.;

	// Artificial viscosity
	double p1,p2,p3,s1,s2,s3,Se,ux1,ux2,ux3,ux,uy1,uy2,uy3,uy,uModule,T,calorv,calorp,c,pmed,rhomed_;

	p1 = nodes_[0]->getCurrentPressure();
	p2 = nodes_[1]->getCurrentPressure();
	p3 = nodes_[2]->getCurrentPressure();

	pmed=(p1+p2+p3)/3.;

	double rho1_ = nodes_[1]->getCurrentDensity();
	double rho2_ = nodes_[1]->getCurrentDensity();
	double rho3_ = nodes_[1]->getCurrentDensity();

	rhomed_ = (rho1_+rho2_+rho3_)/3.;

	s1 = fabs(p1-p2+p1-p3)/(fabs(p1-p2)+fabs(p1-p3)+.00000001);
	s2 = fabs(p2-p1+p2-p3)/(fabs(p2-p1)+fabs(p2-p3)+.00000001);
	s3 = fabs(p3-p1+p3-p2)/(fabs(p3-p1)+fabs(p3-p2)+.00000001);

	Se = (s1+s2+s3)/3.;

	ux1 = nodes_[0]->getCurrentVelocity()(0);
	ux2 = nodes_[1]->getCurrentVelocity()(0);
	ux3 = nodes_[2]->getCurrentVelocity()(0);
	ux = (ux1+ux2+ux3)/3.;
	uy1 = nodes_[0]->getCurrentVelocity()(1);
	uy2 = nodes_[1]->getCurrentVelocity()(1);
	uy3 = nodes_[2]->getCurrentVelocity()(1);
	uy = (uy1+uy2+uy3)/3.;
	uModule = sqrt(ux*ux+uy*uy);

	//std::cout<<h<<" hhhh "<<Se<<" Se "<<uModule<<std::endl;

	calorv = material_ -> getSpecificHeatv();
	calorp = material_ -> getSpecificHeatp();
	//T = (1./calorv)*(e-(u_(0)*u_(0)+u_(1)*u_(1))/2.);
	c = sqrt((calorp/calorv)*(pmed/rhomed_));

	bounded_vector<double, 2> rhou1_ = nodes_[0]->getCurrentMomentum();
	bounded_vector<double, 2> rhou2_ = nodes_[1]->getCurrentMomentum();
	bounded_vector<double, 2> rhou3_ = nodes_[2]->getCurrentMomentum();
		
	double xsi1 = 1./3.;
	double xsi2 = 1./3.;
	vector<double> phi = domainShapeFunction(xsi1, xsi2);
	matrix<double> dphi_dxsi = domainDerivativeShapeFunction(xsi1, xsi2); //row = direction, column = node

	//dy_dxsi at time t+1 = Aint
	double dy1_dxsi1 = 0.0, dy1_dxsi2 = 0.0, dy2_dxsi1 = 0.0, dy2_dxsi2 = 0.0;
	for (size_t i = 0; i < nodes_.size(); i++)
	{
 		bounded_vector<double, 2> currentCoord = nodes_[i]->getCurrentCoordinate();
 		dy1_dxsi1 += currentCoord(0) * dphi_dxsi(0,i);
 		dy1_dxsi2 += currentCoord(0) * dphi_dxsi(1,i);
 		dy2_dxsi1 += currentCoord(1) * dphi_dxsi(0,i);
 		dy2_dxsi2 += currentCoord(1) * dphi_dxsi(1,i);
	}

	//current jacobian
	jac = dy1_dxsi1 * dy2_dxsi2 - dy1_dxsi2 * dy2_dxsi1; 
	double h = sqrt(jac);
	double A = jac/2.;

	double shockVisc_ = qdif * h * (uModule + c) * Se;


	 
	for (int ih = 0; ih < numberOfIntegrationPoints_; ih++) // verificar o number
	{
	 	double xsi1 = domainIntegrationPoints_(ih,0);
	 	double xsi2 = domainIntegrationPoints_(ih,1);
	 	double weight = domainIntegrationPoints_(ih,2);

	 	vector<double> phi = domainShapeFunction(xsi1, xsi2);
	 	matrix<double> dphi_dxsi = domainDerivativeShapeFunction(xsi1, xsi2); //row = direction, column = node

	 	//dy_dxsi at time t+1 = Aint
	 	double dy1_dxsi1 = 0.0, dy1_dxsi2 = 0.0, dy2_dxsi1 = 0.0, dy2_dxsi2 = 0.0;
	 	for (size_t i = 0; i < nodes_.size(); i++)
	 	{
	 		bounded_vector<double, 2> currentCoord = nodes_[i]->getCurrentCoordinate();
	 		dy1_dxsi1 += currentCoord(0) * dphi_dxsi(0,i);
	 		dy1_dxsi2 += currentCoord(0) * dphi_dxsi(1,i);
	 		dy2_dxsi1 += currentCoord(1) * dphi_dxsi(0,i);
	 		dy2_dxsi2 += currentCoord(1) * dphi_dxsi(1,i);
	 	}

	 	//current jacobian
	 	jac = dy1_dxsi1 * dy2_dxsi2 - dy1_dxsi2 * dy2_dxsi1; 

		matrix<double> dxsi_dy(2, 2);
		dxsi_dy(0,0)= dy2_dxsi2/jac;
		dxsi_dy(0,1)= -dy1_dxsi2/jac;
		dxsi_dy(1,0)= -dy2_dxsi1/jac;
		dxsi_dy(1,1)= dy1_dxsi1/jac;

		matrix<double> dphi_dy(2, nodes_.size());
		for (size_t i = 0; i < nodes_.size(); i++)
		{
			dphi_dy(0,i)= dphi_dxsi(0,i)*dxsi_dy(0,0)+dphi_dxsi(1,i)*dxsi_dy(1,0);
			dphi_dy(1,i)= dphi_dxsi(0,i)*dxsi_dy(0,1)+dphi_dxsi(1,i)*dxsi_dy(1,1);
		}

	 	//dv_dxsi at time t+1
	 	/* double dv1_dxsi1 = 0.0, dv1_dxsi2 = 0.0, dv2_dxsi1 = 0.0, dv2_dxsi2 = 0.0;
	 	for (size_t i = 0; i < nodes_.size(); i++)
	 	{
	 		bounded_vector<double, 2> currentVel = nodes_[i]->getCurrentVelocity(); // definir essa função em nodes!!!!!
	 		dv1_dxsi1 += currentVel(0) * dphi_dxsi(0,i);
	 		dv1_dxsi2 += currentVel(0) * dphi_dxsi(1,i);
	 		dv2_dxsi1 += currentVel(1) * dphi_dxsi(0,i);
	 		dv2_dxsi2 += currentVel(1) * dphi_dxsi(1,i);
	 	} */
		//cálculo das variaveis nos pontos de hammer
				
		// cálculo das variaveis nos pontos de integração
		bounded_vector<double, 2> rhou_;
		rhou_(0) = 0.;
		rhou_(1) = 0.;
		double rho_ = 0.;
		double rhoe_ = 0.;
		double e =0.;
		double gamma_= (material_ -> getSpecificHeatp())/(material_ -> getSpecificHeatv());

		//		 material_ -> getUndisturbedTemperature
		
	 	for (size_t i = 0; i < nodes_.size(); i++)
		 {
			 //p_ += nodes_[i]->getCurrentPressure()*phi(i);
			 rho_ += nodes_[i]->getCurrentDensity()*phi(i);
			 rhou_(0) += nodes_[i]->getCurrentMomentum()(0)*phi(i);
			 rhou_(1) += nodes_[i]->getCurrentMomentum()(1)*phi(i);
			 rhoe_ += nodes_[i]->getCurrentInternalEnergy()*phi(i);
		 } 
				
		bounded_vector<double, 2> u_;
		u_(0) = rhou_(0)/rho_;
		u_(1) = rhou_(1)/rho_;
		e = rhoe_/rho_;
		double p_ = (gamma_-1.)*(rhoe_-.5*(rhou_(0)*u_(0) + rhou_(1)*u_(1)));

		// Calculo das derivadas de rhou

		matrix<double> drhou_dy(2, 2);
		drhou_dy.clear();
		bounded_vector<double, 2> drho_dy;
		drho_dy.clear();
		bounded_vector<double, 2> drhoe_dy;
		drhoe_dy.clear();

		for (size_t i = 0; i < nodes_.size(); i++)
		{	
			for (size_t k = 0; k < 2; k++)
			{
				drhoe_dy(k) += nodes_[i]->getCurrentInternalEnergy()*dphi_dy(k,i);

				drho_dy(k) += nodes_[i]->getCurrentDensity()*dphi_dy(k,i);
				
				drhou_dy(0,k) += nodes_[i]->getCurrentMomentum()(0)*dphi_dy(k,i);
				drhou_dy(1,k) += nodes_[i]->getCurrentMomentum()(1)*dphi_dy(k,i);
			}
		}	

		matrix<double> du_dy(2, 2);
		du_dy.clear();

		du_dy(0,0) = (drhou_dy(0,0)*rho_-rhou_(0)*drho_dy(0))/(rho_*rho_);
		du_dy(0,1) = (drhou_dy(0,1)*rho_-rhou_(0)*drho_dy(1))/(rho_*rho_);
		du_dy(1,0) = (drhou_dy(1,0)*rho_-rhou_(1)*drho_dy(0))/(rho_*rho_);
		du_dy(1,1) = (drhou_dy(1,1)*rho_-rhou_(1)*drho_dy(1))/(rho_*rho_);

		bounded_vector<double, 2> dp_dy;
		dp_dy(0) = (gamma_-1.)*(drhoe_dy(0)-.5*(drhou_dy(0,0)*u_(0) + rhou_(0)*du_dy(0,0)+drhou_dy(1,0)*u_(1) + rhou_(1)*du_dy(1,0)));
		dp_dy(1) = (gamma_-1.)*(drhoe_dy(1)-.5*(drhou_dy(0,1)*u_(0) + rhou_(0)*du_dy(0,1)+drhou_dy(1,1)*u_(1) + rhou_(1)*du_dy(1,1)));

		//mesh velocity
		bounded_vector<double, 2> meshU_;
		for (size_t i = 0; i < nodes_.size(); i++)
				{
				  meshU_(0) = 0.; // += nodes_[i]->getCurrentMeshVelocity()(0)*phi(i);
				  meshU_(1)= 0.; // += nodes_[i]->getCurrentMeshVelocity()(1)*phi(i);
				
			 	} 
		/*bounded_vector<double, 2> meshU_;
		meshU_(0) = 0.;
		meshU_(1) = 0.;*/

		//Get viscosity
		double mi_0 = material_->getViscosity();
		double mi;
		double tempInf = material_->getUndisturbedTemperature(); //initial temperature
		double T; // current temperature
		double calorv = material_->getSpecificHeatv();

		//Temperatura

		T=(1./calorv)*(e-(u_(0)*u_(0)+u_(1)*u_(1))/2.);

		if (sutherland==true)
		{
			//mi = mi_0;
			mi = mi_0*pow(((T+273.15)/(tempInf+273.15)),(3./2.))*((tempInf+273.15+110.4)/(T+273.15+110.4));
		}

		else if (sutherland==false)
		{
				mi = mi_0;
		}
		
		// shear stress
		matrix<double> tau_(2, 2);
		tau_.clear();
		matrix<double> dtau_dy(2, 2);
		dtau_dy.clear();

		tau_(0,0) = -2./3. * mi * (du_dy(0,0) + du_dy(1,1)) +2. * mi * du_dy(0,0);
		tau_(0,1) = mi*(du_dy(0,1)+du_dy(1,0));
		tau_(1,0) = mi*(du_dy(1,0)+du_dy(0,1));
		tau_(1,1) = -2./3. * mi * (du_dy(0,0) + du_dy(1,1)) +2. * mi * du_dy(1,1);

		// field forces
		bounded_vector<double, 2> g_;
		g_.clear();
		g_(0) = 0.;
		g_(1) = gravity_;

		//element tangent matrix
	 	for (int i = 0; i < nodes_.size(); i++)
		{
			//element rhs vector
			//momentum
			rhs(2*i) += -deltat_ * phi(i)*(du_dy(0,0)*rhou_(0)+u_(0)*drhou_dy(0,0)+du_dy(1,1)*rhou_(0)+u_(1)*drhou_dy(0,1))* jac* weight;
			rhs(2*i+1) += -deltat_ * phi(i)*(du_dy(0,0)*rhou_(1)+u_(0)*drhou_dy(1,0)+du_dy(1,1)*rhou_(1)+u_(1)*drhou_dy(1,1))* jac * weight;

			//mesh velocity
			rhs(2*i) += deltat_ * phi(i)*(meshU_(0)*drhou_dy(0,0)+meshU_(1)*drhou_dy(0,1)) * jac * weight;
			rhs(2*i+1) += deltat_ * phi(i)*(meshU_(0)*drhou_dy(1,0)+meshU_(1)*drhou_dy(1,1)) * jac * weight;

			rhs(2*i) += -deltat_ * (dphi_dy(0,i) * tau_(0,0)  + dphi_dy(1,i) * tau_(0,1)) * jac * weight;
			rhs(2*i+1) += -deltat_ * (dphi_dy(1,i) * tau_(1,1) + dphi_dy(0,i) * tau_(0,1)) * jac * weight;

			rhs(2*i) += deltat_ * dphi_dy(0,i) * p_ * jac * weight; 
			rhs(2*i+1) += deltat_ * dphi_dy(1,i) * p_ * jac * weight;

			//rhs(2*i) += -deltat_ * phi(i) * dp_dy(0) * jac * weight; //mudar a pressão posteriormente para aplicar no contorno
			//rhs(2*i+1) += -deltat_ * phi(i) * dp_dy(1) * jac * weight;

			
			rhs(2*i)  -= deltat_ * shockVisc_* ( dphi_dy(0,i) * drhou_dy(0,0) +dphi_dy(1,i) * drhou_dy(0,1))* jac * weight; 
			rhs(2*i+1)-= deltat_ * shockVisc_* ( dphi_dy(0,i) * drhou_dy(1,0) +dphi_dy(1,i) * drhou_dy(1,1))* jac * weight; 

			// field force (gravity)
			rhs(2*i) += deltat_*phi(i) * rho_* g_(0) * jac * weight;
			rhs(2*i+1) += deltat_*phi(i) * rho_ * g_(1) * jac * weight;

			// deltat²/2

			rhs(2*i) += -deltat_*deltat_/2. *(phi(i)*(du_dy(0,0)+du_dy(1,1))+dphi_dy(0,i)*u_(0)+dphi_dy(1,i)*u_(1))*
			((du_dy(0,0)*rhou_(0)+u_(0)*drhou_dy(0,0)+du_dy(1,1)*rhou_(0)+u_(1)*drhou_dy(0,1))-(meshU_(0)*drhou_dy(0,0) + meshU_(1) * drhou_dy(0,1)) + dp_dy(0)-rho_* g_(0))* jac * weight;
			
			rhs(2*i+1) += -deltat_*deltat_/2. *(phi(i)*(du_dy(0,0)+du_dy(1,1))+dphi_dy(0,i)*u_(0)+dphi_dy(1,i)*u_(1))*
			((du_dy(0,0)*rhou_(1)+u_(0)*drhou_dy(1,0)+du_dy(1,1)*rhou_(1)+u_(1)*drhou_dy(1,1))-(meshU_(0)*drhou_dy(1,0)+meshU_(1)*drhou_dy(1,1))+dp_dy(1)-rho_* g_(1))* jac * weight;



			for (int j = 0; j < nodes_.size(); j++)
			{
			
				tangent(2*i,2*j) += phi(i) * phi(j) * jac * weight;
				tangent(2*i+1,2*j+1) += phi(i) * phi(j) * jac * weight;


			
			}
		
	 	}
	}



	/*rhs(0) +=  deltat_ *(qdif * Se * (uModule + c) / h)* A*(-rhou1_(0)/6. + rhou2_(0)/12. + rhou3_(0)/12.);
	rhs(1) +=  deltat_ *(qdif * Se * (uModule + c) / h)* A*(-rhou1_(1)/6. + rhou2_(1)/12. + rhou3_(1)/12.);
	rhs(2) +=  deltat_ *(qdif * Se * (uModule + c) / h)* A*(rhou1_(0)/12. - rhou2_(0)/6.  + rhou3_(0)/12.);
	rhs(3) +=  deltat_ *(qdif * Se * (uModule + c) / h)* A*(rhou1_(1)/12. - rhou2_(1)/6.  + rhou3_(1)/12.);
	rhs(4) +=  deltat_ *(qdif * Se * (uModule + c) / h)* A*(rhou1_(0)/12. + rhou2_(0)/12. - rhou3_(0)/6.);
	rhs(5) +=  deltat_ *(qdif * Se * (uModule + c) / h)* A*(rhou1_(1)/12. + rhou2_(1)/12. - rhou3_(1)/6.);
*/
		
	
		std::pair<vector<double>,vector<double>>GaussQuad = gaussLegendre(numberOfIntegrationPoints_);

		//Lado 1
		if((nodes_[0]->getNeummanConstrain()==true)&&(nodes_[1]->getNeummanConstrain()==true)){

			bounded_vector<double, 2> currentCoord0 = nodes_[0]->getCurrentCoordinate();
			bounded_vector<double, 2> currentCoord1 = nodes_[1]->getCurrentCoordinate();

			double sideLenght = sqrt((currentCoord1(0)-currentCoord0(0))*(currentCoord1(0)-currentCoord0(0)) + (currentCoord1(1)-currentCoord0(1))*(currentCoord1(1)-currentCoord0(1)));

			bounded_vector<double, 2> Normal;
			Normal(0)=(currentCoord1(1)-currentCoord0(1))/sideLenght;
			Normal(1)=-(currentCoord1(0)-currentCoord0(0))/sideLenght;

	 		for (int ig = 0; ig < numberOfIntegrationPoints_; ig++) {
				double xsi1 = -GaussQuad.first(ig)/2. + .5;
	 			double xsi2 = GaussQuad.first(ig)/2. + .5;
	 			double weight = GaussQuad.second(ig);

	 			vector<double> phi = domainShapeFunction(xsi1, xsi2);
	 			matrix<double> dphi_dxsi = domainDerivativeShapeFunction(xsi1, xsi2); //row = direction, column = node

	 			//Computing deformation gradient. Below, y refers to current position and x to reference position

		 		//dy_dxsi at time t+1 = Aint
		 		double dy1_dxsi1 = 0.0, dy1_dxsi2 = 0.0, dy2_dxsi1 = 0.0, dy2_dxsi2 = 0.0;
	 				for (size_t i = 0; i < nodes_.size(); i++)
	 				{
	 					bounded_vector<double, 2> currentCoord = nodes_[i]->getCurrentCoordinate();
	 					dy1_dxsi1 += currentCoord(0) * dphi_dxsi(0,i);
	 					dy1_dxsi2 += currentCoord(0) * dphi_dxsi(1,i);
	 					dy2_dxsi1 += currentCoord(1) * dphi_dxsi(0,i);
	 					dy2_dxsi2 += currentCoord(1) * dphi_dxsi(1,i);
	 				}

			 	//current jacobian
			 	double jac = dy1_dxsi1 * dy2_dxsi2 - dy1_dxsi2 * dy2_dxsi1; 

				matrix<double> dxsi_dy(2, 2);
				dxsi_dy(0,0)= dy2_dxsi2/jac;
				dxsi_dy(0,1)= -dy1_dxsi2/jac;
				dxsi_dy(1,0)= -dy2_dxsi1/jac;
				dxsi_dy(1,1)= dy1_dxsi1/jac;

				matrix<double> dphi_dy(2, nodes_.size());
				for (size_t i = 0; i < nodes_.size(); i++)
				{
					dphi_dy(0,i)= dphi_dxsi(0,i)*dxsi_dy(0,0)+dphi_dxsi(1,i)*dxsi_dy(1,0);
					dphi_dy(1,i)= dphi_dxsi(0,i)*dxsi_dy(0,1)+dphi_dxsi(1,i)*dxsi_dy(1,1);
				}

				// cálculo das variaveis nos pontos de integração
				bounded_vector<double, 2> rhou_;
				rhou_(0) = 0.;
				rhou_(1) = 0.;
				double rho_ = 0.;
				double rhoe_ = 0.;
				double temperature_ =0.;
				double gamma_= (material_ -> getSpecificHeatp())/(material_ -> getSpecificHeatv());
		
			 	for (size_t i = 0; i < nodes_.size(); i++)
				{
					//p_ += nodes_[i]->getCurrentPressure()*phi(i);
					rho_ += nodes_[i]->getCurrentDensity()*phi(i);
					rhou_(0) += nodes_[i]->getCurrentMomentum()(0)*phi(i);
				 	rhou_(1) += nodes_[i]->getCurrentMomentum()(1)*phi(i);
				 	rhoe_ += nodes_[i]->getCurrentInternalEnergy()*phi(i);
					temperature_ += nodes_[i]->getCurrentTemperature()*phi(i);
			 	} 
			bounded_vector<double, 2> u_;
			u_(0) = rhou_(0)/rho_;
			u_(1) = rhou_(1)/rho_;
			//double p_ = material_->getDensity()*material_->getUndisturbedTemperature()*(material_ -> getSpecificHeatp()-material_ -> getSpecificHeatv());
			double p_ = material_->getDensity()*temperature_*(material_ -> getSpecificHeatp()-material_ -> getSpecificHeatv());
			//(gamma_-1.)*(rhoe_-.5*(rhou_(0)*u_(0) + rhou_(1)*u_(1)));

			// Calculo das derivadas de rhou

			matrix<double> drhou_dy(2, 2);
			drhou_dy.clear();
			bounded_vector<double, 2> drho_dy;
			drho_dy.clear();

			for (size_t i = 0; i < nodes_.size(); i++)
			{
						
				for (size_t k = 0; k < 2; k++)
				{
					drho_dy(k) += nodes_[i]->getCurrentDensity()*dphi_dy(k,i);

					drhou_dy(0,k) += nodes_[i]->getCurrentMomentum()(0)*dphi_dy(k,i);
					drhou_dy(1,k) += nodes_[i]->getCurrentMomentum()(1)*dphi_dy(k,i);
				}
			}	

			matrix<double> du_dy(2, 2);
			du_dy.clear();

			du_dy(0,0) = (drhou_dy(0,0)*rho_-rhou_(0)*drho_dy(0))/(rho_*rho_);
			du_dy(0,1) = (drhou_dy(0,1)*rho_-rhou_(0)*drho_dy(1))/(rho_*rho_);
			du_dy(1,0) = (drhou_dy(1,0)*rho_-rhou_(1)*drho_dy(0))/(rho_*rho_);
			du_dy(1,1) = (drhou_dy(1,1)*rho_-rhou_(1)*drho_dy(1))/(rho_*rho_);

			///Get viscosity
			double mi_0 = material_->getViscosity();
			double mi;
			double tempInf = material_->getUndisturbedTemperature(); //initial temperature
			double T; // current temperature
			double calorv = material_->getSpecificHeatv();

			//Temperatura
			double e = rhoe_/rho_;
			T=(1./calorv)*(e-(u_(0)*u_(0)+u_(1)*u_(1))/2.);
		
			if (sutherland==true)
			{
				//mi = mi_0;
				mi = mi_0*pow(((T+273.15)/(tempInf+273.15)),(3./2.))*((tempInf+273.15+110.4)/(T+273.15+110.4));
			}
		
			else if (sutherland==false)
			{
				mi = mi_0;
			}
			// shear stress
			matrix<double> tau_(2, 2);
			tau_.clear();
			matrix<double> dtau_dy(2, 2);
			dtau_dy.clear();

			tau_(0,0) = -2./3. * mi * (du_dy(0,0) + du_dy(1,1)) +2. * mi * du_dy(0,0);
			tau_(0,1) = mi*(du_dy(0,1)+du_dy(1,0));
			tau_(1,0) = mi*(du_dy(1,0)+du_dy(0,1));
			tau_(1,1) = -2./3. * mi * (du_dy(0,0) + du_dy(1,1)) +2. * mi * du_dy(1,1);

			for(int i = 0; i < nodes_.size()-1; i++){
			//	rhs(2*i) += deltat_ *phi(i) * ((0.*tau_(0,0)-p_)*Normal(0)+0.*tau_(0,1)*Normal(1))*weight*sideLenght/2.;
			//	rhs(2*i+1) += deltat_ *phi(i) * ((0.*tau_(1,1)-p_)*Normal(1)+0.*tau_(1,0)*Normal(0))*weight*sideLenght/2.;
			}

	 	} // loop over gauss points
	}//if side in Neumann boundary
			//Lado 2
	if((nodes_[1]->getNeummanConstrain()==true)&&(nodes_[2]->getNeummanConstrain()==true)){

	 	
		bounded_vector<double, 2> currentCoord0 = nodes_[1]->getCurrentCoordinate();
		bounded_vector<double, 2> currentCoord1 = nodes_[2]->getCurrentCoordinate();

		double sideLenght = sqrt((currentCoord1(0)-currentCoord0(0))*(currentCoord1(0)-currentCoord0(0)) + (currentCoord1(1)-currentCoord0(1))*(currentCoord1(1)-currentCoord0(1)));

		bounded_vector<double, 2> Normal;
		Normal(0)=(currentCoord1(1)-currentCoord0(1))/sideLenght;
		Normal(1)=-(currentCoord1(0)-currentCoord0(0))/sideLenght;

	 	for (int ig = 0; ig < numberOfIntegrationPoints_; ig++) {

			double xsi1 = 0.; 
	 		double xsi2 = -GaussQuad.first(ig)/2. + .5;

	 		double weight = GaussQuad.second(ig);

	 		vector<double> phi = domainShapeFunction(xsi1, xsi2);
	 		matrix<double> dphi_dxsi = domainDerivativeShapeFunction(xsi1, xsi2); //row = direction, column = node

	 		//Computing deformation gradient. Below, y refers to current position and x to reference position

		 	//dy_dxsi at time t+1 = Aint
		 	double dy1_dxsi1 = 0.0, dy1_dxsi2 = 0.0, dy2_dxsi1 = 0.0, dy2_dxsi2 = 0.0;
	 		for (size_t i = 0; i < nodes_.size(); i++)
	 		{
	 			bounded_vector<double, 2> currentCoord = nodes_[i]->getCurrentCoordinate();
	 			dy1_dxsi1 += currentCoord(0) * dphi_dxsi(0,i);
	 			dy1_dxsi2 += currentCoord(0) * dphi_dxsi(1,i);
	 			dy2_dxsi1 += currentCoord(1) * dphi_dxsi(0,i);
	 			dy2_dxsi2 += currentCoord(1) * dphi_dxsi(1,i);
	 		}

		 	//current jacobian
		 	double jac = dy1_dxsi1 * dy2_dxsi2 - dy1_dxsi2 * dy2_dxsi1; 

			matrix<double> dxsi_dy(2, 2);
			dxsi_dy(0,0)= dy2_dxsi2/jac;
			dxsi_dy(0,1)= -dy1_dxsi2/jac;
			dxsi_dy(1,0)= -dy2_dxsi1/jac;
			dxsi_dy(1,1)= dy1_dxsi1/jac;

			matrix<double> dphi_dy(2, nodes_.size());
			for (size_t i = 0; i < nodes_.size(); i++)
			{
					dphi_dy(0,i)= dphi_dxsi(0,i)*dxsi_dy(0,0)+dphi_dxsi(1,i)*dxsi_dy(1,0);
					dphi_dy(1,i)= dphi_dxsi(0,i)*dxsi_dy(0,1)+dphi_dxsi(1,i)*dxsi_dy(1,1);
			}

			// cálculo das variaveis nos pontos de integração
			bounded_vector<double, 2> rhou_;
			rhou_(0) = 0.;
			rhou_(1) = 0.;
			double rho_ = 0.;
			double rhoe_ = 0.;
			double temperature_ =0.;
			double gamma_= (material_ -> getSpecificHeatp())/(material_ -> getSpecificHeatv());
		
		 	for (size_t i = 0; i < nodes_.size(); i++)
			{
				//p_ += nodes_[i]->getCurrentPressure()*phi(i);
				rho_ += nodes_[i]->getCurrentDensity()*phi(i);
				rhou_(0) += nodes_[i]->getCurrentMomentum()(0)*phi(i);
			 	rhou_(1) += nodes_[i]->getCurrentMomentum()(1)*phi(i);
			 	rhoe_ += nodes_[i]->getCurrentInternalEnergy()*phi(i);
				temperature_ += nodes_[i]->getCurrentTemperature()*phi(i);
		 	} 
			bounded_vector<double, 2> u_;
			u_(0) = rhou_(0)/rho_;
			u_(1) = rhou_(1)/rho_;
			//double p_ = material_->getDensity()*material_->getUndisturbedTemperature()*(material_ -> getSpecificHeatp()-material_ -> getSpecificHeatv());
			double p_ = material_->getDensity()*temperature_*(material_ -> getSpecificHeatp()-material_ -> getSpecificHeatv());
			//double p_ = (gamma_-1.)*(rhoe_-.5*(rhou_(0)*u_(0) + rhou_(1)*u_(1)));

			// Calculo das derivadas de rhou

			matrix<double> drhou_dy(2, 2);
			drhou_dy.clear();
			bounded_vector<double, 2> drho_dy;
			drho_dy.clear();

			for (size_t i = 0; i < nodes_.size(); i++)
			{
						
				for (size_t k = 0; k < 2; k++)
				{
					drho_dy(k) += nodes_[i]->getCurrentDensity()*dphi_dy(k,i);

					drhou_dy(0,k) += nodes_[i]->getCurrentMomentum()(0)*dphi_dy(k,i);
					drhou_dy(1,k) += nodes_[i]->getCurrentMomentum()(1)*dphi_dy(k,i);
				}
			}	

			matrix<double> du_dy(2, 2);
			du_dy.clear();

			du_dy(0,0) = (drhou_dy(0,0)*rho_-rhou_(0)*drho_dy(0))/(rho_*rho_);
			du_dy(0,1) = (drhou_dy(0,1)*rho_-rhou_(0)*drho_dy(1))/(rho_*rho_);
			du_dy(1,0) = (drhou_dy(1,0)*rho_-rhou_(1)*drho_dy(0))/(rho_*rho_);
			du_dy(1,1) = (drhou_dy(1,1)*rho_-rhou_(1)*drho_dy(1))/(rho_*rho_);

			///Get viscosity
			double mi_0 = material_->getViscosity();
			double mi;
			double tempInf = material_->getUndisturbedTemperature(); //initial temperature
			double T; // current temperature
			double calorv = material_->getSpecificHeatv();

			//Temperatura
			double e = rhoe_/rho_;
			T=(1./calorv)*(e-(u_(0)*u_(0)+u_(1)*u_(1))/2.);

			if (sutherland==true)
			{
				//mi = mi_0;
				mi = mi_0*pow(((T+273.15)/(tempInf+273.15)),(3./2.))*((tempInf+273.15+110.4)/(T+273.15+110.4));
			}
		
			else if (sutherland==false)
			{
				mi = mi_0;
			}
			// shear stress
			matrix<double> tau_(2, 2);
			tau_.clear();
			matrix<double> dtau_dy(2, 2);
			dtau_dy.clear();

			tau_(0,0) = -2./3. * mi * (du_dy(0,0) + du_dy(1,1)) +2. * mi * du_dy(0,0);
			tau_(0,1) = mi*(du_dy(0,1)+du_dy(1,0));
			tau_(1,0) = mi*(du_dy(1,0)+du_dy(0,1));
			tau_(1,1) = -2./3. * mi * (du_dy(0,0) + du_dy(1,1)) +2. * mi * du_dy(1,1);

			//for(int i = 1; i < nodes_.size(); i++){

			//	rhs(2*i) += deltat_ *phi(i) * ((0.*tau_(0,0)-p_)*Normal(0)+0.*tau_(0,1)*Normal(1))*weight*sideLenght/2.;
			//	rhs(2*i+1) += deltat_ *phi(i) * ((0.*tau_(1,1)-p_)*Normal(1)+0.*tau_(1,0)*Normal(0))*weight*sideLenght/2.;
			//}
		}
	}
			//Lado 3
	if((nodes_[2]->getNeummanConstrain()==true)&&(nodes_[0]->getNeummanConstrain()==true)){

	 	bounded_vector<double, 2> currentCoord0 = nodes_[2]->getCurrentCoordinate();
		bounded_vector<double, 2> currentCoord1 = nodes_[0]->getCurrentCoordinate();

		double sideLenght = sqrt((currentCoord1(0)-currentCoord0(0))*(currentCoord1(0)-currentCoord0(0)) + (currentCoord1(1)-currentCoord0(1))*(currentCoord1(1)-currentCoord0(1)));

		bounded_vector<double, 2> Normal;
		Normal(0)=(currentCoord1(1)-currentCoord0(1))/sideLenght;
		Normal(1)=-(currentCoord1(0)-currentCoord0(0))/sideLenght;

	 	for (int ig = 0; ig < numberOfIntegrationPoints_; ig++) {

			double xsi1 = GaussQuad.first(ig)/2. + .5;
	 		double xsi2 = 0.;
	 		double weight = GaussQuad.second(ig);

	 		vector<double> phi = domainShapeFunction(xsi1, xsi2);
	 		matrix<double> dphi_dxsi = domainDerivativeShapeFunction(xsi1, xsi2); //row = direction, column = node

	 		//Computing deformation gradient. Below, y refers to current position and x to reference position

		 	//dy_dxsi at time t+1 = Aint
		 	double dy1_dxsi1 = 0.0, dy1_dxsi2 = 0.0, dy2_dxsi1 = 0.0, dy2_dxsi2 = 0.0;
	 		for (size_t i = 0; i < nodes_.size(); i++)
	 		{
	 			bounded_vector<double, 2> currentCoord = nodes_[i]->getCurrentCoordinate();
	 			dy1_dxsi1 += currentCoord(0) * dphi_dxsi(0,i);
	 			dy1_dxsi2 += currentCoord(0) * dphi_dxsi(1,i);
	 			dy2_dxsi1 += currentCoord(1) * dphi_dxsi(0,i);
	 			dy2_dxsi2 += currentCoord(1) * dphi_dxsi(1,i);
	 		}

		 	//current jacobian
		 	double jac = dy1_dxsi1 * dy2_dxsi2 - dy1_dxsi2 * dy2_dxsi1; 

			matrix<double> dxsi_dy(2, 2);
			dxsi_dy(0,0)= dy2_dxsi2/jac;
			dxsi_dy(0,1)= -dy1_dxsi2/jac;
			dxsi_dy(1,0)= -dy2_dxsi1/jac;
			dxsi_dy(1,1)= dy1_dxsi1/jac;

			matrix<double> dphi_dy(2, nodes_.size());
			for (size_t i = 0; i < nodes_.size(); i++)
			{
					dphi_dy(0,i)= dphi_dxsi(0,i)*dxsi_dy(0,0)+dphi_dxsi(1,i)*dxsi_dy(1,0);
					dphi_dy(1,i)= dphi_dxsi(0,i)*dxsi_dy(0,1)+dphi_dxsi(1,i)*dxsi_dy(1,1);
			}

			// cálculo das variaveis nos pontos de integração
			bounded_vector<double, 2> rhou_;
			rhou_(0) = 0.;
			rhou_(1) = 0.;
			double rho_ = 0.;
			double rhoe_ = 0.;
			double temperature_ =0.;
			double gamma_= (material_ -> getSpecificHeatp())/(material_ -> getSpecificHeatv());
		
		 	for (size_t i = 0; i < nodes_.size(); i++)
			{
				//p_ += nodes_[i]->getCurrentPressure()*phi(i);
				rho_ += nodes_[i]->getCurrentDensity()*phi(i);
				rhou_(0) += nodes_[i]->getCurrentMomentum()(0)*phi(i);
			 	rhou_(1) += nodes_[i]->getCurrentMomentum()(1)*phi(i);
			 	rhoe_ += nodes_[i]->getCurrentInternalEnergy()*phi(i);
				temperature_ += nodes_[i]->getCurrentTemperature()*phi(i);
		 	} 
			bounded_vector<double, 2> u_;
			u_(0) = rhou_(0)/rho_;
			u_(1) = rhou_(1)/rho_;
			//double p_ = material_->getDensity()*material_->getUndisturbedTemperature()*(material_ -> getSpecificHeatp()-material_ -> getSpecificHeatv());
			double p_ = material_->getDensity()*temperature_*(material_ -> getSpecificHeatp()-material_ -> getSpecificHeatv());
		 	
			// Calculo das derivadas de rhou

			matrix<double> drhou_dy(2, 2);
			drhou_dy.clear();
			bounded_vector<double, 2> drho_dy;
			drho_dy.clear();

			for (size_t i = 0; i < nodes_.size(); i++)
			{
						
				for (size_t k = 0; k < 2; k++)
				{
					drho_dy(k) += nodes_[i]->getCurrentDensity()*dphi_dy(k,i);

					drhou_dy(0,k) += nodes_[i]->getCurrentMomentum()(0)*dphi_dy(k,i);
					drhou_dy(1,k) += nodes_[i]->getCurrentMomentum()(1)*dphi_dy(k,i);
				}
			}	

			matrix<double> du_dy(2, 2);
			du_dy.clear();

			du_dy(0,0) = (drhou_dy(0,0)*rho_-rhou_(0)*drho_dy(0))/(rho_*rho_);
			du_dy(0,1) = (drhou_dy(0,1)*rho_-rhou_(0)*drho_dy(1))/(rho_*rho_);
			du_dy(1,0) = (drhou_dy(1,0)*rho_-rhou_(1)*drho_dy(0))/(rho_*rho_);
			du_dy(1,1) = (drhou_dy(1,1)*rho_-rhou_(1)*drho_dy(1))/(rho_*rho_);

			///Get viscosity
			double mi_0 = material_->getViscosity();
			double mi;
			double tempInf = material_->getUndisturbedTemperature(); //initial temperature
			double T; // current temperature
			double calorv = material_->getSpecificHeatv();

			//Temperatura
			double e = rhoe_/rho_;
			T=(1./calorv)*(e-(u_(0)*u_(0)+u_(1)*u_(1))/2.);

			if (sutherland==true)
			{
				//mi = mi_0;
				mi = mi_0*pow(((T+273.15)/(tempInf+273.15)),(3./2.))*((tempInf+273.15+110.4)/(T+273.15+110.4));
			}
		
			else if (sutherland==false)
			{
				mi = mi_0;
			}
			// shear stress
			matrix<double> tau_(2, 2);
			tau_.clear();
			matrix<double> dtau_dy(2, 2);
			dtau_dy.clear();

			tau_(0,0) = -2./3. * mi * (du_dy(0,0) + du_dy(1,1)) +2. * mi * du_dy(0,0);
			tau_(0,1) = mi*(du_dy(0,1)+du_dy(1,0));
			tau_(1,0) = mi*(du_dy(1,0)+du_dy(0,1));
			tau_(1,1) = -2./3. * mi * (du_dy(0,0) + du_dy(1,1)) +2. * mi * du_dy(1,1);
	
			//rhs(4) += deltat_ * phi(2) * ((0.*tau_(0,0)-p_)*Normal(0)+0.*tau_(0,1)*Normal(1))*weight*sideLenght/2.;
			//rhs(5) += deltat_ * phi(2) * ((0.*tau_(1,1)-p_)*Normal(1)+0.*tau_(1,0)*Normal(0))*weight*sideLenght/2.;
			//rhs(0) += deltat_ * phi(0) * ((0.*tau_(0,0)-p_)*Normal(0)+0.*tau_(0,1)*Normal(1))*weight*sideLenght/2.;
			//rhs(1) += deltat_ * phi(0) * ((0.*tau_(1,1)-p_)*Normal(1)+0.*tau_(1,0)*Normal(0))*weight*sideLenght/2.;

		}
	} 
return std::make_pair(rhs, tangent);
}

std::pair<bounded_vector<double, 3>, bounded_matrix<double, 3, 3> > Element::elementContributionsS2(const double& qdif, const double& teta_, const bool& sutherland)
{	

	bounded_vector<double, 3> rhs;
	bounded_matrix<double, 3, 3> tangent;     
	rhs.clear(); tangent.clear();

	// auto param = computePSPG();
	// double t_PSPG = param.first; t_PSPG = 0.0;
	// double t_LSIC = param.second;
	double jac = 0.;

	// Artificial viscosity
	double p1,p2,p3,s1,s2,s3,Se,ux1,ux2,ux3,ux,uy1,uy2,uy3,uy,uModule,T,calorv,calorp,c,pmed,rhomed_;

	p1 = nodes_[0]->getCurrentPressure();
	p2 = nodes_[1]->getCurrentPressure();
	p3 = nodes_[2]->getCurrentPressure();

	pmed=(p1+p2+p3)/3.;

	double rho1_ = nodes_[1]->getCurrentDensity();
	double rho2_ = nodes_[1]->getCurrentDensity();
	double rho3_ = nodes_[1]->getCurrentDensity();

	rhomed_ = (rho1_+rho2_+rho3_)/3.;

	s1 = fabs(p1-p2+p1-p3)/(fabs(p1-p2)+fabs(p1-p3)+.00000001);
	s2 = fabs(p2-p1+p2-p3)/(fabs(p2-p1)+fabs(p2-p3)+.00000001);
	s3 = fabs(p3-p1+p3-p2)/(fabs(p3-p1)+fabs(p3-p2)+.00000001);

	Se = (s1+s2+s3)/3.;

	ux1 = nodes_[0]->getCurrentVelocity()(0);
	ux2 = nodes_[1]->getCurrentVelocity()(0);
	ux3 = nodes_[2]->getCurrentVelocity()(0);
	ux = (ux1+ux2+ux3)/3.;
	uy1 = nodes_[0]->getCurrentVelocity()(1);
	uy2 = nodes_[1]->getCurrentVelocity()(1);
	uy3 = nodes_[2]->getCurrentVelocity()(1);
	uy = (uy1+uy2+uy3)/3.;
	uModule = sqrt(ux*ux+uy*uy);

	//std::cout<<h<<" hhhh "<<Se<<" Se "<<uModule<<std::endl;

	calorv = material_ -> getSpecificHeatv();
	calorp = material_ -> getSpecificHeatp();
	//T = (1./calorv)*(e-(u_(0)*u_(0)+u_(1)*u_(1))/2.);
	c = sqrt((calorp/calorv)*(pmed/rhomed_));

	bounded_vector<double, 2> rhou1_ = nodes_[0]->getCurrentMomentum();
	bounded_vector<double, 2> rhou2_ = nodes_[1]->getCurrentMomentum();
	bounded_vector<double, 2> rhou3_ = nodes_[2]->getCurrentMomentum();
		
	double xsi1 = 1./3.;
	double xsi2 = 1./3.;
	vector<double> phi = domainShapeFunction(xsi1, xsi2);
	matrix<double> dphi_dxsi = domainDerivativeShapeFunction(xsi1, xsi2); //row = direction, column = node

	//dy_dxsi at time t+1 = Aint
	double dy1_dxsi1 = 0.0, dy1_dxsi2 = 0.0, dy2_dxsi1 = 0.0, dy2_dxsi2 = 0.0;
	for (size_t i = 0; i < nodes_.size(); i++)
	{
 		bounded_vector<double, 2> currentCoord = nodes_[i]->getCurrentCoordinate();
 		dy1_dxsi1 += currentCoord(0) * dphi_dxsi(0,i);
 		dy1_dxsi2 += currentCoord(0) * dphi_dxsi(1,i);
 		dy2_dxsi1 += currentCoord(1) * dphi_dxsi(0,i);
 		dy2_dxsi2 += currentCoord(1) * dphi_dxsi(1,i);
	}

	//current jacobian
	jac = dy1_dxsi1 * dy2_dxsi2 - dy1_dxsi2 * dy2_dxsi1; 
	double h = sqrt(jac);
	double A = jac/2.;

	double shockVisc_ = qdif * h * (uModule + c) * Se;

	for (int ih = 0; ih < numberOfIntegrationPoints_; ih++)
	{
	
		double xsi1 = domainIntegrationPoints_(ih,0);
		double xsi2 = domainIntegrationPoints_(ih,1);
		double weight = domainIntegrationPoints_(ih,2);

		vector<double> phi = domainShapeFunction(xsi1, xsi2);
		matrix<double> dphi_dxsi = domainDerivativeShapeFunction(xsi1, xsi2); //row = direction, column = node

		//Computing deformation gradient. Below, y refers to current position and x to reference position

		//dy_dxsi at time t+1 = Aint
		double dy1_dxsi1 = 0.0, dy1_dxsi2 = 0.0, dy2_dxsi1 = 0.0, dy2_dxsi2 = 0.0;
		for (size_t i = 0; i < nodes_.size(); i++)
		{
			bounded_vector<double, 2> currentCoord = nodes_[i]->getCurrentCoordinate();
			dy1_dxsi1 += currentCoord(0) * dphi_dxsi(0,i);
			dy1_dxsi2 += currentCoord(0) * dphi_dxsi(1,i);
			dy2_dxsi1 += currentCoord(1) * dphi_dxsi(0,i);
			dy2_dxsi2 += currentCoord(1) * dphi_dxsi(1,i);
		}

		//current jacobian
		jac = dy1_dxsi1 * dy2_dxsi2 - dy1_dxsi2 * dy2_dxsi1; 

	    matrix<double> dxsi_dy(2, 2);
		dxsi_dy(0,0)= dy2_dxsi2/jac;
		dxsi_dy(0,1)= -dy1_dxsi2/jac;
		dxsi_dy(1,0)= -dy2_dxsi1/jac;
		dxsi_dy(1,1)= dy1_dxsi1/jac;


		matrix<double> dphi_dy(2, nodes_.size());
		for (int i = 0; i < nodes_.size(); i++)
		{
			dphi_dy(0,i) = dxsi_dy(0,0) * dphi_dxsi(0,i) + dxsi_dy(1,0) * dphi_dxsi(1,i);
			dphi_dy(1,i) = dxsi_dy(0,1) * dphi_dxsi(0,i) + dxsi_dy(1,1) * dphi_dxsi(1,i);
		}


		//Pressure gradient dp_dy
		bounded_vector<double, 2> dp_dy; dp_dy.clear();
		for (int i = 0; i < nodes_.size(); i++)
		{
			dp_dy(0) += dphi_dy(0,i) * nodes_[i]->getCurrentPressure();   // definir as funções de pressão em nodes!!!!!
			dp_dy(1) += dphi_dy(1,i) * nodes_[i]->getCurrentPressure();
		}

		//std::cout<< "Passo 1" << std::endl;

		
		// cálculo das variaveis nos pontos de integração
		bounded_vector<double, 2> rhou_;
		bounded_vector<double, 2> rhou2_;
		bounded_vector<double, 2> deltaRhou_;
		rhou_(0) = 0.;
		rhou_(1) = 0.;
		rhou2_(0) = 0.;
		rhou2_(1) = 0.;
		deltaRhou_(0) = 0.;
		deltaRhou_(1) = 0.;
		double rho_ = 0.;
		double rhoe_ = 0.;
		double e = 0.;

		//std::cout<< "Passo 2" << std::endl;

	 	for (size_t i = 0; i < nodes_.size(); i++)
		 {
			 rho_ += nodes_[i]->getCurrentDensity()*phi(i);
			 rhou_(0) += nodes_[i]->getCurrentMomentum()(0)*phi(i);
			 rhou_(1) += nodes_[i]->getCurrentMomentum()(1)*phi(i);
			 //rhou2_(0) += (nodes_[i]->getCurrentMomentum()(0)+nodes_[i]->getDeltaMomentum()(0))*phi(i);
			 //rhou2_(1) += (nodes_[i]->getCurrentMomentum()(1)+nodes_[i]->getDeltaMomentum()(1))*phi(i);
			 deltaRhou_(0) += nodes_[i]->getDeltaMomentum()(0)*phi(i);
			 deltaRhou_(1) += nodes_[i]->getDeltaMomentum()(1)*phi(i);
			 rhoe_ += nodes_[i]->getCurrentInternalEnergy()*phi(i);
		 } 
		bounded_vector<double, 2> u_;
		u_(0) = rhou_(0)/rho_;
		u_(1) = rhou_(1)/rho_;
		e = rhoe_/rho_;

		// Calculo das derivadas de rhou
		matrix<double> drhou_dy(2, 2);
		drhou_dy.clear();
		bounded_vector<double, 2> drho_dy;
		drho_dy.clear();

		for (size_t i = 0; i < nodes_.size(); i++)
		{	
			for (size_t k = 0; k < 2; k++)
			{
				drho_dy(k) += nodes_[i]->getCurrentDensity()*dphi_dy(k,i);
				
				drhou_dy(0,k) += nodes_[i]->getCurrentMomentum()(0)*dphi_dy(k,i);
				drhou_dy(1,k) += nodes_[i]->getCurrentMomentum()(1)*dphi_dy(k,i);
			}
		}	

		matrix<double> du_dy(2, 2);
		du_dy.clear();

		du_dy(0,0) = (drhou_dy(0,0)*rho_-rhou_(0)*drho_dy(0))/(rho_*rho_);
		du_dy(0,1) = (drhou_dy(0,1)*rho_-rhou_(0)*drho_dy(1))/(rho_*rho_);
		du_dy(1,0) = (drhou_dy(1,0)*rho_-rhou_(1)*drho_dy(0))/(rho_*rho_);
		du_dy(1,1) = (drhou_dy(1,1)*rho_-rhou_(1)*drho_dy(1))/(rho_*rho_);

		//mesh velocity
		bounded_vector<double, 2> meshU_;
		bounded_vector<double, 2> rhoMeshU_;
		meshU_(0) = 0.;
		meshU_(1) = 0.;
		for (size_t i = 0; i < nodes_.size(); i++)
				{
				  meshU_(0) += nodes_[i]->getCurrentMeshVelocity()(0)*phi(i);
				  meshU_(1) += nodes_[i]->getCurrentMeshVelocity()(1)*phi(i);
				
			 	} 
		
		matrix<double> dmeshU_dy(2, 2);
		dmeshU_dy(0,0)=0.; dmeshU_dy(0,1)=0.; dmeshU_dy(1,0)=0.; dmeshU_dy(1,1)=0.;
		for (size_t i = 0; i < nodes_.size(); i++)
		{
						
			for (size_t k = 0; k < 2; k++)
			{
			  dmeshU_dy(0,k)+=  nodes_[i]->getCurrentMeshVelocity()(0)*dphi_dy(k,i);
			  dmeshU_dy(1,k)+=  nodes_[i]->getCurrentMeshVelocity()(1)*dphi_dy(k,i);
			}
		}	
		rhoMeshU_(0) = rho_*meshU_(0);
		rhoMeshU_(1) = rho_*meshU_(1);
		
		/*
		bounded_vector<double, 2> meshU_;
		bounded_vector<double, 2> rhoMeshU_;
		matrix<double> dmeshU_dy(2,2);
		meshU_(0) = 0.;
		meshU_(1) = 0.;
		dmeshU_dy(0,0) = 0.;
		dmeshU_dy(0,1) = 0.;
		dmeshU_dy(1,0) = 0.;
		dmeshU_dy(1,1) = 0.;
		rhoMeshU_(0) = 0.;
		rhoMeshU_(1) = 0.;*/

		//dv_dxsi at time t+1
	 	double dv1_dxsi1 = 0.0, dv1_dxsi2 = 0.0, dv2_dxsi1 = 0.0, dv2_dxsi2 = 0.0;
	 	for (size_t i = 0; i < nodes_.size(); i++)
	 	{
	 		bounded_vector<double, 2> currentVel = nodes_[i]->getCurrentVelocity(); // definir essa função em nodes!!!!!
	 		dv1_dxsi1 += currentVel(0) * dphi_dxsi(0,i);
	 		dv1_dxsi2 += currentVel(0) * dphi_dxsi(1,i);
	 		dv2_dxsi1 += currentVel(1) * dphi_dxsi(0,i);
	 		dv2_dxsi2 += currentVel(1) * dphi_dxsi(1,i);
	 	}

		// field forces
		bounded_vector<double, 2> g_;
		g_.clear();
		g_(0) = 0.;
		g_(1) = gravity_;

		//element rhs vector
		for (size_t i = 0; i < nodes_.size(); i++)
		{
			// 1 parcela			
			rhs (i) += deltat_*(dphi_dy(0,i)*(rhou_(0) + teta_*deltaRhou_(0)) + dphi_dy(1,i)*(rhou_(1) + teta_*deltaRhou_(1))) * jac * weight;

			// 2 parcela

			rhs(i) -= deltat_*rho_*(dphi_dy(0,i)*meshU_(0) + phi(i)*dmeshU_dy(0,0)+dphi_dy(1,i)*meshU_(1) + phi(i)*dmeshU_dy(1,1))* jac* weight; //confirmar a mudança

			rhs(i) -= deltat_ * shockVisc_ * (dphi_dy(0,i) * drho_dy(0) + dphi_dy(1,i) * drho_dy(1)) * jac * weight;

			// deltat²/2

			//rhs(i) += deltat_*deltat_/2. *((dphi_dy(0,i)*u_(0)+dphi_dy(1,i)*u_(1)+phi(i)*(du_dy(0,0)+du_dy(1,1)))*meshU_(0)*drho_dy(0)+
			//(dphi_dy(0,i)*u_(0)+dphi_dy(1,i)*u_(1)+phi(i)*(du_dy(0,0)+du_dy(1,1)))*meshU_(1)*drho_dy(1));


			for (int j = 0; j < nodes_.size(); j++)
			{
			
				tangent(i,j) += phi(i) * phi(j) * jac * weight;
			
			}
		}
	}	


		//artificialMi = qdif * h * (uModule + c) * Se;
//		rhs(0) +=  deltat_ *(qdif * Se * (uModule + c) / h) * A * (-rho1_/6. + rho2_/12. + rho3_/12.);
//		rhs(1) +=  deltat_ *(qdif * Se * (uModule + c) / h) * A * (rho1_/12. - rho2_/6.  + rho3_/12.);
//		rhs(2) +=  deltat_ *(qdif * Se * (uModule + c) / h) * A * (rho1_/12. + rho2_/12. - rho3_/6.);
		

	std::pair<vector<double>,vector<double>>GaussQuad = gaussLegendre(numberOfIntegrationPoints_);

	//Lado 1
	if(((nodes_[0]->getNeummanConstrain()==true)&&(nodes_[1]->getNeummanConstrain()==true))||((nodes_[0]->getConstrain()==true)&&(nodes_[1]->getConstrain()==true))){

		bounded_vector<double, 2> currentCoord0 = nodes_[0]->getCurrentCoordinate();
		bounded_vector<double, 2> currentCoord1 = nodes_[1]->getCurrentCoordinate();

		double sideLenght = sqrt((currentCoord1(0)-currentCoord0(0))*(currentCoord1(0)-currentCoord0(0)) + (currentCoord1(1)-currentCoord0(1))*(currentCoord1(1)-currentCoord0(1)));

		bounded_vector<double, 2> Normal;
		Normal(0)=(currentCoord1(1)-currentCoord0(1))/sideLenght;
		Normal(1)=-(currentCoord1(0)-currentCoord0(0))/sideLenght;

	 	for (int ig = 0; ig < numberOfIntegrationPoints_; ig++) {

			double xsi1 = -GaussQuad.first(ig)/2. + .5;
	 		double xsi2 = GaussQuad.first(ig)/2. + .5;
	 		double weight = GaussQuad.second(ig);
			
	 		vector<double> phi = domainShapeFunction(xsi1, xsi2);
	 		matrix<double> dphi_dxsi = domainDerivativeShapeFunction(xsi1, xsi2); //row = direction, column = node

	 		//Computing deformation gradient. Below, y refers to current position and x to reference position
		 	//dy_dxsi at time t+1 = Aint
		 	double dy1_dxsi1 = 0.0, dy1_dxsi2 = 0.0, dy2_dxsi1 = 0.0, dy2_dxsi2 = 0.0;
	 		for (size_t i = 0; i < nodes_.size(); i++)
	 		{
	 			bounded_vector<double, 2> currentCoord = nodes_[i]->getCurrentCoordinate();
	 			dy1_dxsi1 += currentCoord(0) * dphi_dxsi(0,i);
	 			dy1_dxsi2 += currentCoord(0) * dphi_dxsi(1,i);
	 			dy2_dxsi1 += currentCoord(1) * dphi_dxsi(0,i);
	 			dy2_dxsi2 += currentCoord(1) * dphi_dxsi(1,i);
	 		}

		 	//current jacobian
		 	double jac = dy1_dxsi1 * dy2_dxsi2 - dy1_dxsi2 * dy2_dxsi1; 

			matrix<double> dxsi_dy(2, 2);
			dxsi_dy(0,0)= dy2_dxsi2/jac;
			dxsi_dy(0,1)= -dy1_dxsi2/jac;
			dxsi_dy(1,0)= -dy2_dxsi1/jac;
			dxsi_dy(1,1)= dy1_dxsi1/jac;

			matrix<double> dphi_dy(2, nodes_.size());
			for (size_t i = 0; i < nodes_.size(); i++)
			{
					dphi_dy(0,i)= dphi_dxsi(0,i)*dxsi_dy(0,0)+dphi_dxsi(1,i)*dxsi_dy(1,0);
					dphi_dy(1,i)= dphi_dxsi(0,i)*dxsi_dy(0,1)+dphi_dxsi(1,i)*dxsi_dy(1,1);
			}

			// cálculo das variaveis nos pontos de integração
			// cálculo das variaveis nos pontos de integração
			bounded_vector<double, 2> rhou_;
			bounded_vector<double, 2> rhou2_;
			bounded_vector<double, 2> deltaRhou_;
			rhou_(0) = 0.;
			rhou_(1) = 0.;
			rhou2_(0) = 0.;
			rhou2_(1) = 0.;
			deltaRhou_(0) = 0.;
			deltaRhou_(1) = 0.;
			double rho_ = 0.;
		
		 	for (size_t i = 0; i < nodes_.size(); i++)
			{

				rho_ += nodes_[i]->getCurrentDensity()*phi(i);
				rhou_(0) += nodes_[i]->getCurrentMomentum()(0)*phi(i);
			 	rhou_(1) += nodes_[i]->getCurrentMomentum()(1)*phi(i);
				deltaRhou_(0) += nodes_[i]->getDeltaMomentum()(0)*phi(i);
			 	deltaRhou_(1) += nodes_[i]->getDeltaMomentum()(1)*phi(i);
		 	} 
			bounded_vector<double, 2> u_;
			u_(0) = rhou_(0)/rho_;
			u_(1) = rhou_(1)/rho_;

			// Calculo das derivadas de rhou
			matrix<double> drhou_dy(2, 2);
			drhou_dy.clear();
			bounded_vector<double, 2> drho_dy;
			drho_dy.clear();

			for (size_t i = 0; i < nodes_.size(); i++)
			{						
				for (size_t k = 0; k < 2; k++)
				{
					drho_dy(k) += nodes_[i]->getCurrentDensity()*dphi_dy(k,i);

					drhou_dy(0,k) += nodes_[i]->getCurrentMomentum()(0)*dphi_dy(k,i);
					drhou_dy(1,k) += nodes_[i]->getCurrentMomentum()(1)*dphi_dy(k,i);
				}
			}	

			//mesh velocity

			bounded_vector<double, 2> meshU_;
			bounded_vector<double, 2> rhoMeshU_;
			meshU_(0) = 0.;
			meshU_(1) = 0.;
			for (size_t i = 0; i < nodes_.size(); i++)
				{
					meshU_(0) += nodes_[i]->getCurrentMeshVelocity()(0)*phi(i);
				 	meshU_(1) += nodes_[i]->getCurrentMeshVelocity()(1)*phi(i);
				
			 	} 
			matrix<double> dmeshU_dy(2, 2);
			for (size_t i = 0; i < nodes_.size(); i++)
			{
						
				for (size_t k = 0; k < 2; k++)
				{
					dmeshU_dy(0,k) += meshU_(0)*dphi_dy(k,i);
					dmeshU_dy(1,k) += meshU_(1)*dphi_dy(k,i);
				}
			}	
			rhoMeshU_(0) = rho_*meshU_(0);
			rhoMeshU_(1) = rho_*meshU_(1);
		
			/*bounded_vector<double, 2> meshU_;
			matrix<double> dmeshU_dy(2,2);
			bounded_vector<double, 2> rhoMeshU_;
			meshU_(0) = 0.;
			meshU_(1) = 0.;
			dmeshU_dy(0,0) = 0.;
			dmeshU_dy(0,1) = 0.;
			dmeshU_dy(1,0) = 0.;
			dmeshU_dy(1,1) = 0.;
			rhoMeshU_(0) = rho_*meshU_(0);
			rhoMeshU_(1) = rho_*meshU_(1);*/


			for(int i = 0; i < nodes_.size()-1; i++){
				rhs(i) += -deltat_ *phi(i) * (Normal(0)*(rhou_(0)+teta_*deltaRhou_(0)-rhoMeshU_(0))+Normal(1)*(rhou_(1) + teta_*deltaRhou_(1)-rhoMeshU_(1)))*weight*sideLenght/2.;
				
				// parcela deltat²/2

			//	rhs(i) += -deltat_*deltat_/2. *phi(i) * (Normal(0)*meshU_(0)*drho_dy(0)+Normal(1)*meshU_(1)*drho_dy(1))*weight*sideLenght/2.;
			}

	 	} // loop over gauss points
	}//if side in Neumann boundary
			//Lado 2
	if(((nodes_[1]->getNeummanConstrain()==true)&&(nodes_[2]->getNeummanConstrain()==true))||((nodes_[1]->getConstrain()==true)&&(nodes_[2]->getConstrain()==true))){

		bounded_vector<double, 2> currentCoord0 = nodes_[1]->getCurrentCoordinate();
		bounded_vector<double, 2> currentCoord1 = nodes_[2]->getCurrentCoordinate();

		double sideLenght = sqrt((currentCoord1(0)-currentCoord0(0))*(currentCoord1(0)-currentCoord0(0)) + (currentCoord1(1)-currentCoord0(1))*(currentCoord1(1)-currentCoord0(1)));

		bounded_vector<double, 2> Normal;
		Normal(0)=(currentCoord1(1)-currentCoord0(1))/sideLenght;
		Normal(1)=-(currentCoord1(0)-currentCoord0(0))/sideLenght;

	 	for (int ig = 0; ig < numberOfIntegrationPoints_; ig++) {

			double xsi1 = 0.; 
	 		double xsi2 = -GaussQuad.first(ig)/2. + .5;
			 double weight = GaussQuad.second(ig);

	 		vector<double> phi = domainShapeFunction(xsi1, xsi2);
	 		matrix<double> dphi_dxsi = domainDerivativeShapeFunction(xsi1, xsi2); //row = direction, column = node

	 		//Computing deformation gradient. Below, y refers to current position and x to reference position

		 	//dy_dxsi at time t+1 = Aint
		 	double dy1_dxsi1 = 0.0, dy1_dxsi2 = 0.0, dy2_dxsi1 = 0.0, dy2_dxsi2 = 0.0;
	 		for (size_t i = 0; i < nodes_.size(); i++)
	 		{
	 			bounded_vector<double, 2> currentCoord = nodes_[i]->getCurrentCoordinate();
	 			dy1_dxsi1 += currentCoord(0) * dphi_dxsi(0,i);
	 			dy1_dxsi2 += currentCoord(0) * dphi_dxsi(1,i);
	 			dy2_dxsi1 += currentCoord(1) * dphi_dxsi(0,i);
	 			dy2_dxsi2 += currentCoord(1) * dphi_dxsi(1,i);
	 		}

		 	//current jacobian
		 	double jac = dy1_dxsi1 * dy2_dxsi2 - dy1_dxsi2 * dy2_dxsi1; 

			matrix<double> dxsi_dy(2, 2);
			dxsi_dy(0,0)= dy2_dxsi2/jac;
			dxsi_dy(0,1)= -dy1_dxsi2/jac;
			dxsi_dy(1,0)= -dy2_dxsi1/jac;
			dxsi_dy(1,1)= dy1_dxsi1/jac;

			matrix<double> dphi_dy(2, nodes_.size());
			for (size_t i = 0; i < nodes_.size(); i++)
			{
					dphi_dy(0,i)= dphi_dxsi(0,i)*dxsi_dy(0,0)+dphi_dxsi(1,i)*dxsi_dy(1,0);
					dphi_dy(1,i)= dphi_dxsi(0,i)*dxsi_dy(0,1)+dphi_dxsi(1,i)*dxsi_dy(1,1);
			}

			// cálculo das variaveis nos pontos de integração
			// cálculo das variaveis nos pontos de integração
			bounded_vector<double, 2> rhou_;
			bounded_vector<double, 2> rhou2_;
			bounded_vector<double, 2> deltaRhou_;
			rhou_(0) = 0.;
			rhou_(1) = 0.;
			rhou2_(0) = 0.;
			rhou2_(1) = 0.;
			deltaRhou_(0) = 0.;
			deltaRhou_(1) = 0.;
			double rho_ = 0.;
		
		 	for (size_t i = 0; i < nodes_.size(); i++)
			{

				rho_ += nodes_[i]->getCurrentDensity()*phi(i);
				rhou_(0) += nodes_[i]->getCurrentMomentum()(0)*phi(i);
			 	rhou_(1) += nodes_[i]->getCurrentMomentum()(1)*phi(i);
				deltaRhou_(0) += nodes_[i]->getDeltaMomentum()(0)*phi(i);
			 	deltaRhou_(1) += nodes_[i]->getDeltaMomentum()(1)*phi(i);
		 	} 
			bounded_vector<double, 2> u_;
			u_(0) = rhou_(0)/rho_;
			u_(1) = rhou_(1)/rho_;

			// Calculo das derivadas de rhou

			matrix<double> drhou_dy(2, 2);
			drhou_dy.clear();
			bounded_vector<double, 2> drho_dy;
			drho_dy.clear();

			for (size_t i = 0; i < nodes_.size(); i++)
			{
						
				for (size_t k = 0; k < 2; k++)
				{
					drho_dy(k) += nodes_[i]->getCurrentDensity()*dphi_dy(k,i);

					drhou_dy(0,k) += nodes_[i]->getCurrentMomentum()(0)*dphi_dy(k,i);
					drhou_dy(1,k) += nodes_[i]->getCurrentMomentum()(1)*dphi_dy(k,i);
				}
			}	


			//mesh velocity
			bounded_vector<double, 2> meshU_;
			bounded_vector<double, 2> rhoMeshU_;
			meshU_(0) = 0.;
			meshU_(1) = 0.;
			for (size_t i = 0; i < nodes_.size(); i++)
				{
					meshU_(0) += nodes_[i]->getCurrentMeshVelocity()(0)*phi(i);
				 	meshU_(1) += nodes_[i]->getCurrentMeshVelocity()(1)*phi(i);
				
			 	} 
			matrix<double> dmeshU_dy(2, 2);
			for (size_t i = 0; i < nodes_.size(); i++)
			{
						
				for (size_t k = 0; k < 2; k++)
				{
					dmeshU_dy(0,k) += meshU_(0)*dphi_dy(k,i);
					dmeshU_dy(1,k) += meshU_(1)*dphi_dy(k,i);
				}
			}	
			rhoMeshU_(0) = rho_*meshU_(0);
			rhoMeshU_(1) = rho_*meshU_(1);
			
			/*bounded_vector<double, 2> meshU_;
			matrix<double> dmeshU_dy(2,2);
			bounded_vector<double, 2> rhoMeshU_;
			meshU_(0) = 0.;
			meshU_(1) = 0.;
			dmeshU_dy(0,0) = 0.;
			dmeshU_dy(0,1) = 0.;
			dmeshU_dy(1,0) = 0.;
			dmeshU_dy(1,1) = 0.;
			rhoMeshU_(0) = rho_*meshU_(0);
			rhoMeshU_(1) = rho_*meshU_(1);*/
			

			for(int i = 1; i < nodes_.size(); i++){

				rhs(i) += -deltat_ *phi(i) * (Normal(0)*(rhou_(0)+teta_*deltaRhou_(0)-rhoMeshU_(0))+Normal(1)*(rhou_(1)+teta_*deltaRhou_(1)-rhoMeshU_(1)))*weight*sideLenght/2.;
				
				// parcela deltat²/2

				//rhs(i) += -deltat_*deltat_/2. *phi(i) * (Normal(0)*meshU_(0)*drho_dy(0)+Normal(1)*meshU_(1)*drho_dy(1))*weight*sideLenght/2.;
			
			}

	 	} // loop over gauss points
	}
			//Lado 3
	if(((nodes_[2]->getNeummanConstrain()==true)&&(nodes_[0]->getNeummanConstrain()==true))||((nodes_[2]->getConstrain()==true)&&(nodes_[0]->getConstrain()==true))){

	 	bounded_vector<double, 2> currentCoord0 = nodes_[2]->getCurrentCoordinate();
		bounded_vector<double, 2> currentCoord1 = nodes_[0]->getCurrentCoordinate();

		double sideLenght = sqrt((currentCoord1(0)-currentCoord0(0))*(currentCoord1(0)-currentCoord0(0)) + (currentCoord1(1)-currentCoord0(1))*(currentCoord1(1)-currentCoord0(1)));

		bounded_vector<double, 2> Normal;
		Normal(0)=(currentCoord1(1)-currentCoord0(1))/sideLenght;
		Normal(1)=-(currentCoord1(0)-currentCoord0(0))/sideLenght;

	 	for (int ig = 0; ig < numberOfIntegrationPoints_; ig++) {
			double xsi1 = GaussQuad.first(ig)/2. + .5;
	 		double xsi2 = 0.;
	 		double weight = GaussQuad.second(ig);

				 		vector<double> phi = domainShapeFunction(xsi1, xsi2);
	 		matrix<double> dphi_dxsi = domainDerivativeShapeFunction(xsi1, xsi2); //row = direction, column = node

	 		//Computing deformation gradient. Below, y refers to current position and x to reference position

		 	//dy_dxsi at time t+1 = Aint
		 	double dy1_dxsi1 = 0.0, dy1_dxsi2 = 0.0, dy2_dxsi1 = 0.0, dy2_dxsi2 = 0.0;
	 		for (size_t i = 0; i < nodes_.size(); i++)
	 		{
	 			bounded_vector<double, 2> currentCoord = nodes_[i]->getCurrentCoordinate();
	 			dy1_dxsi1 += currentCoord(0) * dphi_dxsi(0,i);
	 			dy1_dxsi2 += currentCoord(0) * dphi_dxsi(1,i);
	 			dy2_dxsi1 += currentCoord(1) * dphi_dxsi(0,i);
	 			dy2_dxsi2 += currentCoord(1) * dphi_dxsi(1,i);
	 		}

		 	//current jacobian
		 	double jac = dy1_dxsi1 * dy2_dxsi2 - dy1_dxsi2 * dy2_dxsi1; 

			matrix<double> dxsi_dy(2, 2);
			dxsi_dy(0,0)= dy2_dxsi2/jac;
			dxsi_dy(0,1)= -dy1_dxsi2/jac;
			dxsi_dy(1,0)= -dy2_dxsi1/jac;
			dxsi_dy(1,1)= dy1_dxsi1/jac;

			matrix<double> dphi_dy(2, nodes_.size());
			for (size_t i = 0; i < nodes_.size(); i++)
			{
					dphi_dy(0,i)= dphi_dxsi(0,i)*dxsi_dy(0,0)+dphi_dxsi(1,i)*dxsi_dy(1,0);
					dphi_dy(1,i)= dphi_dxsi(0,i)*dxsi_dy(0,1)+dphi_dxsi(1,i)*dxsi_dy(1,1);
			}

			// cálculo das variaveis nos pontos de integração
			// cálculo das variaveis nos pontos de integração
			bounded_vector<double, 2> rhou_;
			bounded_vector<double, 2> rhou2_;
			bounded_vector<double, 2> deltaRhou_;
			rhou_(0) = 0.;
			rhou_(1) = 0.;
			rhou2_(0) = 0.;
			rhou2_(1) = 0.;
			deltaRhou_(0) = 0.;
			deltaRhou_(1) = 0.;
			double rho_ = 0.;
		
		 	for (size_t i = 0; i < nodes_.size(); i++)
			{

				rho_ += nodes_[i]->getCurrentDensity()*phi(i);
				rhou_(0) += nodes_[i]->getCurrentMomentum()(0)*phi(i);
			 	rhou_(1) += nodes_[i]->getCurrentMomentum()(1)*phi(i);
				deltaRhou_(0) += nodes_[i]->getDeltaMomentum()(0)*phi(i);
			 	deltaRhou_(1) += nodes_[i]->getDeltaMomentum()(1)*phi(i);
		 	} 
			bounded_vector<double, 2> u_;
			u_(0) = rhou_(0)/rho_;
			u_(1) = rhou_(1)/rho_;

			// Calculo das derivadas de rhou

			matrix<double> drhou_dy(2, 2);
			drhou_dy.clear();
			bounded_vector<double, 2> drho_dy;
			drho_dy.clear();

			for (size_t i = 0; i < nodes_.size(); i++)
			{
						
				for (size_t k = 0; k < 2; k++)
				{
					drho_dy(k) += nodes_[i]->getCurrentDensity()*dphi_dy(k,i);

					drhou_dy(0,k) += nodes_[i]->getCurrentMomentum()(0)*dphi_dy(k,i);
					drhou_dy(1,k) += nodes_[i]->getCurrentMomentum()(1)*dphi_dy(k,i);
				}
			}	


			//mesh velocity
			bounded_vector<double, 2> meshU_;
			bounded_vector<double, 2> rhoMeshU_;
			meshU_(0) = 0.;
			meshU_(1) = 0.;
			for (size_t i = 0; i < nodes_.size(); i++)
				{
					meshU_(0) += nodes_[i]->getCurrentMeshVelocity()(0)*phi(i);
				 	meshU_(1) += nodes_[i]->getCurrentMeshVelocity()(1)*phi(i);
				
			 	} 
			matrix<double> dmeshU_dy(2, 2);
			for (size_t i = 0; i < nodes_.size(); i++)
			{
						
				for (size_t k = 0; k < 2; k++)
				{
					dmeshU_dy(0,k) += meshU_(0)*dphi_dy(k,i);
					dmeshU_dy(1,k) += meshU_(1)*dphi_dy(k,i);
				}
			}	
			rhoMeshU_(0) = rho_*meshU_(0);
			rhoMeshU_(1) = rho_*meshU_(1);
		
			/*bounded_vector<double, 2> meshU_;
			matrix<double> dmeshU_dy(2,2);
			bounded_vector<double, 2> rhoMeshU_;
			meshU_(0) = 0.;
			meshU_(1) = 0.;
			dmeshU_dy(0,0) = 0.;
			dmeshU_dy(0,1) = 0.;
			dmeshU_dy(1,0) = 0.;
			dmeshU_dy(1,1) = 0.;
			rhoMeshU_(0) = rho_*meshU_(0);
			rhoMeshU_(1) = rho_*meshU_(1);*/
		
			rhs(2) += -deltat_ *phi(2) * (Normal(0)*(rhou_(0)+teta_*deltaRhou_(0)-rhoMeshU_(0))+Normal(1)*(rhou_(1)+teta_*deltaRhou_(1)-rhoMeshU_(1)))*weight*sideLenght/2.;
			rhs(0) += -deltat_ *phi(0) * (Normal(0)*(rhou_(0)+teta_*deltaRhou_(0)-rhoMeshU_(0))+Normal(1)*(rhou_(1)+teta_*deltaRhou_(1)-rhoMeshU_(1)))*weight*sideLenght/2.;

			// parcela deltat²/2

			//rhs(2) += -deltat_*deltat_/2. *phi(2) * (Normal(0)*meshU_(0)*drho_dy(0)+Normal(1)*meshU_(1)*drho_dy(1))*weight*sideLenght/2.;
			//rhs(0) += -deltat_*deltat_/2. *phi(0) * (Normal(0)*meshU_(0)*drho_dy(0)+Normal(1)*meshU_(1)*drho_dy(1))*weight*sideLenght/2.;


	 	} // loop over gauss points
	}
	

	return std::make_pair(rhs, tangent);
}

std::pair<bounded_vector<double, 3>, bounded_matrix<double, 3, 3> > Element::elementContributionsS3(const double& qdif, const double& teta_, const bool& sutherland)
{	

	bounded_vector<double, 3> rhs;
	bounded_matrix<double, 3, 3> tangent;     
	rhs.clear(); tangent.clear();
	double jac=0.;

// Artificial viscosity
	double p1,p2,p3,s1,s2,s3,Se,ux1,ux2,ux3,ux,uy1,uy2,uy3,uy,uModule,T,calorv,calorp,c,pmed,rhomed_;

	p1 = nodes_[0]->getCurrentPressure();
	p2 = nodes_[1]->getCurrentPressure();
	p3 = nodes_[2]->getCurrentPressure();

	pmed=(p1+p2+p3)/3.;

	double rho1_ = nodes_[1]->getCurrentDensity();
	double rho2_ = nodes_[1]->getCurrentDensity();
	double rho3_ = nodes_[1]->getCurrentDensity();

	rhomed_ = (rho1_+rho2_+rho3_)/3.;

	s1 = fabs(p1-p2+p1-p3)/(fabs(p1-p2)+fabs(p1-p3)+.00000001);
	s2 = fabs(p2-p1+p2-p3)/(fabs(p2-p1)+fabs(p2-p3)+.00000001);
	s3 = fabs(p3-p1+p3-p2)/(fabs(p3-p1)+fabs(p3-p2)+.00000001);

	Se = (s1+s2+s3)/3.;

	ux1 = nodes_[0]->getCurrentVelocity()(0);
	ux2 = nodes_[1]->getCurrentVelocity()(0);
	ux3 = nodes_[2]->getCurrentVelocity()(0);
	ux = (ux1+ux2+ux3)/3.;
	uy1 = nodes_[0]->getCurrentVelocity()(1);
	uy2 = nodes_[1]->getCurrentVelocity()(1);
	uy3 = nodes_[2]->getCurrentVelocity()(1);
	uy = (uy1+uy2+uy3)/3.;
	uModule = sqrt(ux*ux+uy*uy);

	//std::cout<<h<<" hhhh "<<Se<<" Se "<<uModule<<std::endl;

	calorv = material_ -> getSpecificHeatv();
	calorp = material_ -> getSpecificHeatp();
	//T = (1./calorv)*(e-(u_(0)*u_(0)+u_(1)*u_(1))/2.);
	c = sqrt((calorp/calorv)*(pmed/rhomed_));

	bounded_vector<double, 2> rhou1_ = nodes_[0]->getCurrentMomentum();
	bounded_vector<double, 2> rhou2_ = nodes_[1]->getCurrentMomentum();
	bounded_vector<double, 2> rhou3_ = nodes_[2]->getCurrentMomentum();
		
	double xsi1 = 1./3.;
	double xsi2 = 1./3.;
	vector<double> phi = domainShapeFunction(xsi1, xsi2);
	matrix<double> dphi_dxsi = domainDerivativeShapeFunction(xsi1, xsi2); //row = direction, column = node

	//dy_dxsi at time t+1 = Aint
	double dy1_dxsi1 = 0.0, dy1_dxsi2 = 0.0, dy2_dxsi1 = 0.0, dy2_dxsi2 = 0.0;
	for (size_t i = 0; i < nodes_.size(); i++)
	{
 		bounded_vector<double, 2> currentCoord = nodes_[i]->getCurrentCoordinate();
 		dy1_dxsi1 += currentCoord(0) * dphi_dxsi(0,i);
 		dy1_dxsi2 += currentCoord(0) * dphi_dxsi(1,i);
 		dy2_dxsi1 += currentCoord(1) * dphi_dxsi(0,i);
 		dy2_dxsi2 += currentCoord(1) * dphi_dxsi(1,i);
	}

	//current jacobian
	jac = dy1_dxsi1 * dy2_dxsi2 - dy1_dxsi2 * dy2_dxsi1; 
	double h = sqrt(jac);
	double A = jac/2.;

	//artificial viscosity for shock capturing
	double shockVisc_ = qdif * h * (uModule + c) * Se;

	for (int ih = 0; ih < numberOfIntegrationPoints_; ih++)
	{
		double xsi1 = domainIntegrationPoints_(ih,0);
		double xsi2 = domainIntegrationPoints_(ih,1);
		double weight = domainIntegrationPoints_(ih,2);

		vector<double> phi = domainShapeFunction(xsi1, xsi2);
		matrix<double> dphi_dxsi = domainDerivativeShapeFunction(xsi1, xsi2); //row = direction, column = node

		//Computing deformation gradient. Below, y refers to current position and x to reference position

		//dy_dxsi at time t+1 = Aint
		double dy1_dxsi1 = 0.0, dy1_dxsi2 = 0.0, dy2_dxsi1 = 0.0, dy2_dxsi2 = 0.0;
		for (size_t i = 0; i < nodes_.size(); i++)
		{
			bounded_vector<double, 2> currentCoord = nodes_[i]->getCurrentCoordinate();
			dy1_dxsi1 += currentCoord(0) * dphi_dxsi(0,i);
			dy1_dxsi2 += currentCoord(0) * dphi_dxsi(1,i);
			dy2_dxsi1 += currentCoord(1) * dphi_dxsi(0,i);
			dy2_dxsi2 += currentCoord(1) * dphi_dxsi(1,i);
		}

		//current jacobian
		jac = dy1_dxsi1 * dy2_dxsi2 - dy1_dxsi2 * dy2_dxsi1; 

		matrix<double> dxsi_dy(2, 2);
		dxsi_dy(0,0)= dy2_dxsi2/jac;
		dxsi_dy(0,1)= -dy1_dxsi2/jac;
		dxsi_dy(1,0)= -dy2_dxsi1/jac;
		dxsi_dy(1,1)= dy1_dxsi1/jac;

		matrix<double> dphi_dy(2, nodes_.size());
		for (int i = 0; i < nodes_.size(); i++)
		{
			dphi_dy(0,i) = dxsi_dy(0,0) * dphi_dxsi(0,i) + dxsi_dy(1,0) * dphi_dxsi(1,i);
			dphi_dy(1,i) = dxsi_dy(0,1) * dphi_dxsi(0,i) + dxsi_dy(1,1) * dphi_dxsi(1,i);
		}

		bounded_vector<double, 2> rhou_;
		rhou_(0) = 0.;
		rhou_(1) = 0.;
		double rho_ = 0.;
		double p_ = 0.;
		double rhoe_ = 0.;
		bounded_vector<double, 2> drhoe_dy;
		drhoe_dy.clear();
		bounded_vector<double, 2> drho_dy;
		drho_dy.clear();
		matrix<double> drhou_dy(2, 2);
		drhou_dy.clear();
		bounded_vector<double, 2> dp_dy;
		dp_dy.clear();

		//std::cout<< "Passo 2" << std::endl;

		//Interpolating main variables and calculating its derivatives
	 	for (size_t i = 0; i < nodes_.size(); i++)
		 {
			p_ += nodes_[i]->getCurrentPressure()*phi(i);
			rho_ += nodes_[i]->getCurrentDensity()*phi(i);
			rhou_(0) += nodes_[i]->getCurrentMomentum()(0)*phi(i);
			rhou_(1) += nodes_[i]->getCurrentMomentum()(1)*phi(i);
			rhoe_ += nodes_[i]->getCurrentInternalEnergy()*phi(i);
			for (int k=0;k<2;k++){
				drho_dy(k) += nodes_[i]->getCurrentDensity()*dphi_dy(k,i);
				drhou_dy(0,k) += nodes_[i]->getCurrentMomentum()(0)*dphi_dy(k,i);
				drhou_dy(1,k) += nodes_[i]->getCurrentMomentum()(1)*dphi_dy(k,i);
				drhoe_dy(k) += nodes_[i]->getCurrentInternalEnergy()*dphi_dy(k,i);
				dp_dy(k) += dphi_dy(k,i) * nodes_[i]->getCurrentPressure();
			}
		} 

		//Calculating secondary variables

		//Get viscosity
		double mi_0 = material_->getViscosity();
		double mi;
		bounded_vector<double, 2> u_; //velocity
		double tempInf = material_->getUndisturbedTemperature(); //previous temperature
		double T; // current temperature
		double e; //specific energy
		matrix<double> du_dy(2, 2); //velocity gradient
		matrix<double> tau_(2, 2); // shear stress
		bounded_vector<double, 2> de_dy; //energy gradient
		bounded_vector<double, 2> dT_dy; //temperature gradient
		
		//Temperatura
		e=rhoe_/rho_;
		for (size_t i = 0; i < 2; i++){
			de_dy(i)=(drhoe_dy(i)*rho_-drho_dy(i)*rhoe_)/(rho_*rho_);
		}
		double calorv = material_->getSpecificHeatv();

		T=(1./calorv)*(e-(u_(0)*u_(0)+u_(1)*u_(1))/2.);
		//std::cout<< true <<std::endl;
		if (sutherland==true)
		{
			//mi = mi_0;
			mi = mi_0*pow(((T+273.15)/(tempInf+273.15)),(3./2.))*((tempInf+273.15+110.4)/(T+273.15+110.4));
		}
		
		else if (sutherland==false)
		{
			mi = mi_0;
		}
	
		u_(0) = rhou_(0)/rho_;
		u_(1) = rhou_(1)/rho_;
		du_dy(0,0) = (drhou_dy(0,0)*rho_-rhou_(0)*drho_dy(0))/(rho_*rho_);
		du_dy(0,1) = (drhou_dy(0,1)*rho_-rhou_(0)*drho_dy(1))/(rho_*rho_);
		du_dy(1,0) = (drhou_dy(1,0)*rho_-rhou_(1)*drho_dy(0))/(rho_*rho_);
		du_dy(1,1) = (drhou_dy(1,1)*rho_-rhou_(1)*drho_dy(1))/(rho_*rho_);
	
		tau_(0,0) = -2./3. * mi * (du_dy(0,0) + du_dy(1,1)) +2. * mi * du_dy(0,0);
		tau_(0,1) = mi*(du_dy(0,1)+du_dy(1,0));
		tau_(1,0) = mi*(du_dy(1,0)+du_dy(0,1));
		tau_(1,1) = -2./3. * mi * (du_dy(0,0) + du_dy(1,1)) +2. * mi * du_dy(1,1);

		for (size_t i = 0; i < 2; i++){
			dT_dy(i)=(1./calorv)*(de_dy(i)-u_(0)*du_dy(0,i)-u_(1)*du_dy(1,i));
		}
		
		// determinação do Kt
		double Kt_ = material_->getThermalCond();

		//mesh velocity
		bounded_vector<double, 2> meshU_;
		bounded_vector<double, 2> rhoMeshU_;
		meshU_(0) = 0.;
		meshU_(1) = 0.;
		for (size_t i = 0; i < nodes_.size(); i++)
				{
					meshU_(0) += nodes_[i]->getCurrentMeshVelocity()(0)*phi(i);
				 	meshU_(1) += nodes_[i]->getCurrentMeshVelocity()(1)*phi(i);
				
			 	} 
		matrix<double> dmeshU_dy(2, 2);
		for (size_t i = 0; i < nodes_.size(); i++)
		{
					
			for (size_t k = 0; k < 2; k++)
			{
				dmeshU_dy(0,k) += meshU_(0)*dphi_dy(k,i);
				dmeshU_dy(1,k) += meshU_(1)*dphi_dy(k,i);
			}
		}	
		rhoMeshU_(0) = rho_*meshU_(0);
		rhoMeshU_(1) = rho_*meshU_(1);
		
		/*bounded_vector<double, 2> meshU_;
		matrix<double> dmeshU_dy(2,2);
		meshU_(0) = 0.;
		meshU_(1) = 0.;
		dmeshU_dy(0,0) = 0.;
		dmeshU_dy(0,1) = 0.;
		dmeshU_dy(1,0) = 0.;
		dmeshU_dy(1,1) = 0.;*/
		
		
		//element rhs vector
		for (size_t i = 0; i < nodes_.size(); i++)
		{
			
			// 1 parcela
			rhs (i) += deltat_ * phi(i) * (-(du_dy(0,0)*(rhoe_+p_) + u_(0)*(drhoe_dy(0)+dp_dy(0))+du_dy(1,1)*(rhoe_+p_)
 			+ u_(1)*(drhoe_dy(1)+dp_dy(1))) + meshU_(0)*drhoe_dy(0) + meshU_(1)*drhoe_dy(1))* jac * weight;

			// 2 parcela
			rhs (i) -= deltat_ * (dphi_dy(0,i)*(tau_(0,0)*u_(0)+tau_(0,1)*u_(1)+ Kt_*dT_dy(0))+dphi_dy(1,i)*(tau_(1,0)*u_(0)+tau_(1,1)*u_(1)+Kt_*dT_dy(1)))* jac * weight;
		
			// 3 parcela deltat²/2

			rhs (i) += deltat_*deltat_/2. *((du_dy(0,0)+du_dy(1,1))*phi(i)+dphi_dy(0,i)*u_(0)+dphi_dy(1,i)*u_(1))*(-(du_dy(0,0)*(rhoe_+p_) + u_(0)*(drhoe_dy(0)+dp_dy(0))+du_dy(1,1)*(rhoe_+p_)
 			+ u_(1)*(drhoe_dy(1)+dp_dy(1))) + meshU_(0)*drhoe_dy(0) + meshU_(1)*drhoe_dy(1))* jac * weight;

			//shock capturing contribution
			rhs (i) -= deltat_ * shockVisc_ * (dphi_dy(0,i) * drhoe_dy(0) + dphi_dy(1,i) * drhoe_dy(1))* jac * weight;

			for (int j = 0; j < nodes_.size(); j++)
			{
			
				tangent(i,j) += phi(i) * phi(j) * jac * weight;
			
			}
		}
	}	

		//artificialMi = qdif * h * (uModule + c) * Se;
/*
		rhs(0) +=  deltat_ *(qdif * Se * (uModule + c) / h) * A * (-rhoe1_/6.  + rhoe2_/12. + rhoe3_/12.);
		rhs(1) +=  deltat_ *(qdif * Se * (uModule + c) / h) * A * ( rhoe1_/12. - rhoe2_/6.  + rhoe3_/12.);
		rhs(2) +=  deltat_ *(qdif * Se * (uModule + c) / h) * A * ( rhoe1_/12. + rhoe2_/12. - rhoe3_/6.);
*/

	std::pair<vector<double>,vector<double>>GaussQuad = gaussLegendre(numberOfIntegrationPoints_);

	//Lado 1
	if(((nodes_[0]->getNeummanConstrain()==true)&&(nodes_[1]->getNeummanConstrain()==true))||((nodes_[0]->getConstrain()==true)&&(nodes_[1]->getConstrain()==true))){


		bounded_vector<double, 2> currentCoord0 = nodes_[0]->getCurrentCoordinate();
		bounded_vector<double, 2> currentCoord1 = nodes_[1]->getCurrentCoordinate();

		double sideLenght = sqrt((currentCoord1(0)-currentCoord0(0))*(currentCoord1(0)-currentCoord0(0)) + (currentCoord1(1)-currentCoord0(1))*(currentCoord1(1)-currentCoord0(1)));

		bounded_vector<double, 2> Normal;
		Normal(0)=(currentCoord1(1)-currentCoord0(1))/sideLenght;
		Normal(1)=-(currentCoord1(0)-currentCoord0(0))/sideLenght;

	 	for (int ig = 0; ig < numberOfIntegrationPoints_; ig++) {

			double xsi1 = -GaussQuad.first(ig)/2. + .5;
	 		double xsi2 = GaussQuad.first(ig)/2. + .5;
			
	 		double weight = GaussQuad.second(ig);

	 		vector<double> phi = domainShapeFunction(xsi1, xsi2);
	 		matrix<double> dphi_dxsi = domainDerivativeShapeFunction(xsi1, xsi2); //row = direction, column = node

	 		//Computing deformation gradient. Below, y refers to current position and x to reference position

		 	//dy_dxsi at time t+1 = Aint
		 	double dy1_dxsi1 = 0.0, dy1_dxsi2 = 0.0, dy2_dxsi1 = 0.0, dy2_dxsi2 = 0.0;
	 		for (size_t i = 0; i < nodes_.size(); i++)
	 		{
	 			bounded_vector<double, 2> currentCoord = nodes_[i]->getCurrentCoordinate();
	 			dy1_dxsi1 += currentCoord(0) * dphi_dxsi(0,i);
	 			dy1_dxsi2 += currentCoord(0) * dphi_dxsi(1,i);
	 			dy2_dxsi1 += currentCoord(1) * dphi_dxsi(0,i);
	 			dy2_dxsi2 += currentCoord(1) * dphi_dxsi(1,i);
	 		}

		 	//current jacobian
		 	double jac = dy1_dxsi1 * dy2_dxsi2 - dy1_dxsi2 * dy2_dxsi1; 

			matrix<double> dxsi_dy(2, 2);
			dxsi_dy(0,0)= dy2_dxsi2/jac;
			dxsi_dy(0,1)= -dy1_dxsi2/jac;
			dxsi_dy(1,0)= -dy2_dxsi1/jac;
			dxsi_dy(1,1)= dy1_dxsi1/jac;

			matrix<double> dphi_dy(2, nodes_.size());
			for (size_t i = 0; i < nodes_.size(); i++)
			{
					dphi_dy(0,i)= dphi_dxsi(0,i)*dxsi_dy(0,0)+dphi_dxsi(1,i)*dxsi_dy(1,0);
					dphi_dy(1,i)= dphi_dxsi(0,i)*dxsi_dy(0,1)+dphi_dxsi(1,i)*dxsi_dy(1,1);
			}

			// cálculo das variaveis nos pontos de integração
			bounded_vector<double, 2> rhou_;
			rhou_(0) = 0.;
			rhou_(1) = 0.;
			double rho_ = 0.;
			double rhoe_ = 0.;
			bounded_vector<double, 2> drhoe_dy;
			drhoe_dy.clear();
			bounded_vector<double, 2> drho_dy;
			drho_dy.clear();
			matrix<double> drhou_dy(2, 2);
			drhou_dy.clear();
			double gamma_= (material_ -> getSpecificHeatp())/(material_ -> getSpecificHeatv());
		
		 	for (size_t i = 0; i < nodes_.size(); i++)
			{

				rho_ += nodes_[i]->getCurrentDensity()*phi(i);
				rhou_(0) += nodes_[i]->getCurrentMomentum()(0)*phi(i);
				rhou_(1) += nodes_[i]->getCurrentMomentum()(1)*phi(i);
				rhoe_ += nodes_[i]->getCurrentInternalEnergy()*phi(i);
				for (int k=0;k<2;k++){
					drho_dy(k) += nodes_[i]->getCurrentDensity()*dphi_dy(k,i);
					drhou_dy(0,k) += nodes_[i]->getCurrentMomentum()(0)*dphi_dy(k,i);
					drhou_dy(1,k) += nodes_[i]->getCurrentMomentum()(1)*dphi_dy(k,i);
					drhoe_dy(k) += nodes_[i]->getCurrentInternalEnergy()*dphi_dy(k,i);
				}
		 	} 
			
		///Get viscosity
			double mi_0 = material_->getViscosity();
			double mi;
			double tempInf = material_->getUndisturbedTemperature(); //initial temperature
			double T; // current temperature
			double calorv = material_->getSpecificHeatv();
			bounded_vector<double, 2> u_; //velocity
			u_(0) = rhou_(0)/rho_;
			u_(1) = rhou_(1)/rho_;

			//Temperatura
			double e = rhoe_/rho_;
			T=(1./calorv)*(e-(u_(0)*u_(0)+u_(1)*u_(1))/2.);

			if (sutherland==true)
			{
				//mi = mi_0;
				mi = mi_0*pow(((T+273.15)/(tempInf+273.15)),(3./2.))*((tempInf+273.15+110.4)/(T+273.15+110.4));
			}
		
			else if (sutherland==false)
			{
				mi = mi_0;
			}
			
		matrix<double> du_dy(2, 2); //velocity gradient
		matrix<double> tau_(2, 2); // shear stress
		bounded_vector<double, 2> de_dy; //energy gradient
		bounded_vector<double, 2> dT_dy; //temperature gradient

		du_dy(0,0) = (drhou_dy(0,0)*rho_-rhou_(0)*drho_dy(0))/(rho_*rho_);
		du_dy(0,1) = (drhou_dy(0,1)*rho_-rhou_(0)*drho_dy(1))/(rho_*rho_);
		du_dy(1,0) = (drhou_dy(1,0)*rho_-rhou_(1)*drho_dy(0))/(rho_*rho_);
		du_dy(1,1) = (drhou_dy(1,1)*rho_-rhou_(1)*drho_dy(1))/(rho_*rho_);
	
		tau_(0,0) = -2./3. * mi * (du_dy(0,0) + du_dy(1,1)) +2. * mi * du_dy(0,0);
		tau_(0,1) = mi*(du_dy(0,1)+du_dy(1,0));
		tau_(1,0) = mi*(du_dy(1,0)+du_dy(0,1));
		tau_(1,1) = -2./3. * mi * (du_dy(0,0) + du_dy(1,1)) +2. * mi * du_dy(1,1);
		//Temperatura

		for (size_t i = 0; i < 2; i++){
			de_dy(i)=(drhoe_dy(i)*rho_-drho_dy(i)*rhoe_)/(rho_*rho_);
		}

		for (size_t i = 0; i < 2; i++){
			dT_dy(i)=(1./calorv)*(de_dy(i)-u_(0)*du_dy(0,i)-u_(1)*du_dy(1,i));
		}
		
		// determinação do Kt
		double Kt_ = material_->getThermalCond();

			for(int i = 0; i < nodes_.size()-1; i++){

				//rhs(i) += deltat_ *phi(i) * (Normal(0)*(tau_(0,0)*u_(0)+tau_(0,1)*u_(1)+ Kt_*dT_dy(0))+Normal(1)*(tau_(1,0)*u_(0)+tau_(1,1)*u_(1)+Kt_*dT_dy(1)))*weight*sideLenght/2.;
				
			}

	 	} // loop over gauss points
	}//if side in Neumann boundary
			//Lado 2
	if(((nodes_[1]->getNeummanConstrain()==true)&&(nodes_[2]->getNeummanConstrain()==true))||((nodes_[1]->getConstrain()==true)&&(nodes_[2]->getConstrain()==true))){

	 	
		bounded_vector<double, 2> currentCoord0 = nodes_[1]->getCurrentCoordinate();
		bounded_vector<double, 2> currentCoord1 = nodes_[2]->getCurrentCoordinate();

		double sideLenght = sqrt((currentCoord1(0)-currentCoord0(0))*(currentCoord1(0)-currentCoord0(0)) + (currentCoord1(1)-currentCoord0(1))*(currentCoord1(1)-currentCoord0(1)));

		bounded_vector<double, 2> Normal;
		Normal(0)=(currentCoord1(1)-currentCoord0(1))/sideLenght;
		Normal(1)=-(currentCoord1(0)-currentCoord0(0))/sideLenght;

	 	for (int ig = 0; ig < numberOfIntegrationPoints_; ig++) {

			double xsi1 = 0.; 
	 		double xsi2 = -GaussQuad.first(ig)/2. + .5;
	 		double weight = GaussQuad.second(ig);

	 		vector<double> phi = domainShapeFunction(xsi1, xsi2);
	 		matrix<double> dphi_dxsi = domainDerivativeShapeFunction(xsi1, xsi2); //row = direction, column = node

	 		//Computing deformation gradient. Below, y refers to current position and x to reference position

		 	//dy_dxsi at time t+1 = Aint
		 	double dy1_dxsi1 = 0.0, dy1_dxsi2 = 0.0, dy2_dxsi1 = 0.0, dy2_dxsi2 = 0.0;
	 		for (size_t i = 0; i < nodes_.size(); i++)
	 		{
	 			bounded_vector<double, 2> currentCoord = nodes_[i]->getCurrentCoordinate();
	 			dy1_dxsi1 += currentCoord(0) * dphi_dxsi(0,i);
	 			dy1_dxsi2 += currentCoord(0) * dphi_dxsi(1,i);
	 			dy2_dxsi1 += currentCoord(1) * dphi_dxsi(0,i);
	 			dy2_dxsi2 += currentCoord(1) * dphi_dxsi(1,i);
	 		}

		 	//current jacobian
		 	double jac = dy1_dxsi1 * dy2_dxsi2 - dy1_dxsi2 * dy2_dxsi1; 

			matrix<double> dxsi_dy(2, 2);
			dxsi_dy(0,0)= dy2_dxsi2/jac;
			dxsi_dy(0,1)= -dy1_dxsi2/jac;
			dxsi_dy(1,0)= -dy2_dxsi1/jac;
			dxsi_dy(1,1)= dy1_dxsi1/jac;

			matrix<double> dphi_dy(2, nodes_.size());
			for (size_t i = 0; i < nodes_.size(); i++)
			{
					dphi_dy(0,i)= dphi_dxsi(0,i)*dxsi_dy(0,0)+dphi_dxsi(1,i)*dxsi_dy(1,0);
					dphi_dy(1,i)= dphi_dxsi(0,i)*dxsi_dy(0,1)+dphi_dxsi(1,i)*dxsi_dy(1,1);
			}

			// cálculo das variaveis nos pontos de integração
			bounded_vector<double, 2> rhou_;
			rhou_(0) = 0.;
			rhou_(1) = 0.;
			double rho_ = 0.;
			double rhoe_ = 0.;
			bounded_vector<double, 2> drhoe_dy;
			drhoe_dy.clear();
			bounded_vector<double, 2> drho_dy;
			drho_dy.clear();
			matrix<double> drhou_dy(2, 2);
			drhou_dy.clear();
			double gamma_= (material_ -> getSpecificHeatp())/(material_ -> getSpecificHeatv());
		
		 	for (size_t i = 0; i < nodes_.size(); i++)
			{

				rho_ += nodes_[i]->getCurrentDensity()*phi(i);
				rhou_(0) += nodes_[i]->getCurrentMomentum()(0)*phi(i);
				rhou_(1) += nodes_[i]->getCurrentMomentum()(1)*phi(i);
				rhoe_ += nodes_[i]->getCurrentInternalEnergy()*phi(i);
				for (int k=0;k<2;k++){
					drho_dy(k) += nodes_[i]->getCurrentDensity()*dphi_dy(k,i);
					drhou_dy(0,k) += nodes_[i]->getCurrentMomentum()(0)*dphi_dy(k,i);
					drhou_dy(1,k) += nodes_[i]->getCurrentMomentum()(1)*dphi_dy(k,i);
					drhoe_dy(k) += nodes_[i]->getCurrentInternalEnergy()*dphi_dy(k,i);
				}
		 	} 
			
		//Get viscosity
		double mi_0 = material_->getViscosity();
		double mi;
		double tempInf = material_->getUndisturbedTemperature(); //initial temperature
		double T; // current temperature
		double calorv = material_->getSpecificHeatv();
		bounded_vector<double, 2> u_; //velocity
		u_(0) = rhou_(0)/rho_;
		u_(1) = rhou_(1)/rho_;

		//Temperatura
		double e = rhoe_/rho_;
		T=(1./calorv)*(e-(u_(0)*u_(0)+u_(1)*u_(1))/2.);

		if (sutherland==true)
		{
			//mi = mi_0;
			mi = mi_0*pow(((T+273.15)/(tempInf+273.15)),(3./2.))*((tempInf+273.15+110.4)/(T+273.15+110.4));
		}
		
		else if (sutherland==false)
		{
				mi = mi_0;
		}

		matrix<double> du_dy(2, 2); //velocity gradient
		matrix<double> tau_(2, 2); // shear stress
		bounded_vector<double, 2> de_dy; //energy gradient
		bounded_vector<double, 2> dT_dy; //temperature gradient


		du_dy(0,0) = (drhou_dy(0,0)*rho_-rhou_(0)*drho_dy(0))/(rho_*rho_);
		du_dy(0,1) = (drhou_dy(0,1)*rho_-rhou_(0)*drho_dy(1))/(rho_*rho_);
		du_dy(1,0) = (drhou_dy(1,0)*rho_-rhou_(1)*drho_dy(0))/(rho_*rho_);
		du_dy(1,1) = (drhou_dy(1,1)*rho_-rhou_(1)*drho_dy(1))/(rho_*rho_);
	
		tau_(0,0) = -2./3. * mi * (du_dy(0,0) + du_dy(1,1)) +2. * mi * du_dy(0,0);
		tau_(0,1) = mi*(du_dy(0,1)+du_dy(1,0));
		tau_(1,0) = mi*(du_dy(1,0)+du_dy(0,1));
		tau_(1,1) = -2./3. * mi * (du_dy(0,0) + du_dy(1,1)) +2. * mi * du_dy(1,1);

		//Temperatura

		for (size_t i = 0; i < 2; i++){
			de_dy(i)=(drhoe_dy(i)*rho_-drho_dy(i)*rhoe_)/(rho_*rho_);
		}

		for (size_t i = 0; i < 2; i++){
			dT_dy(i)=(1./calorv)*(de_dy(i)-u_(0)*du_dy(0,i)-u_(1)*du_dy(1,i));
		}
		
		// determinação do Kt
		double Kt_ = material_->getThermalCond();

			for(int i = 1; i < nodes_.size(); i++){

				//rhs(i) += deltat_ *phi(i) * (Normal(0)*(tau_(0,0)*u_(0)+tau_(0,1)*u_(1)+ Kt_*dT_dy(0))+Normal(1)*(tau_(1,0)*u_(0)+tau_(1,1)*u_(1)+Kt_*dT_dy(1)))*weight*sideLenght/2.;
				
			}

	 	} // loop over gauss points
	}
			//Lado 3
	if(((nodes_[2]->getNeummanConstrain()==true)&&(nodes_[0]->getNeummanConstrain()==true))||((nodes_[2]->getConstrain()==true)&&(nodes_[0]->getConstrain()==true))){

	 	bounded_vector<double, 2> currentCoord0 = nodes_[2]->getCurrentCoordinate();
		bounded_vector<double, 2> currentCoord1 = nodes_[0]->getCurrentCoordinate();

		double sideLenght = sqrt((currentCoord1(0)-currentCoord0(0))*(currentCoord1(0)-currentCoord0(0)) + (currentCoord1(1)-currentCoord0(1))*(currentCoord1(1)-currentCoord0(1)));

		bounded_vector<double, 2> Normal;
		Normal(0)=(currentCoord1(1)-currentCoord0(1))/sideLenght;
		Normal(1)=-(currentCoord1(0)-currentCoord0(0))/sideLenght;

	 	for (int ig = 0; ig < numberOfIntegrationPoints_; ig++) {
			double xsi1 = GaussQuad.first(ig)/2. + .5;
	 		double xsi2 = 0.;
			double weight = GaussQuad.second(ig);

	 		vector<double> phi = domainShapeFunction(xsi1, xsi2);
	 		matrix<double> dphi_dxsi = domainDerivativeShapeFunction(xsi1, xsi2); //row = direction, column = node

	 		//Computing deformation gradient. Below, y refers to current position and x to reference position

		 	//dy_dxsi at time t+1 = Aint
		 	double dy1_dxsi1 = 0.0, dy1_dxsi2 = 0.0, dy2_dxsi1 = 0.0, dy2_dxsi2 = 0.0;
	 		for (size_t i = 0; i < nodes_.size(); i++)
	 		{
	 			bounded_vector<double, 2> currentCoord = nodes_[i]->getCurrentCoordinate();
	 			dy1_dxsi1 += currentCoord(0) * dphi_dxsi(0,i);
	 			dy1_dxsi2 += currentCoord(0) * dphi_dxsi(1,i);
	 			dy2_dxsi1 += currentCoord(1) * dphi_dxsi(0,i);
	 			dy2_dxsi2 += currentCoord(1) * dphi_dxsi(1,i);
	 		}

		 	//current jacobian
		 	double jac = dy1_dxsi1 * dy2_dxsi2 - dy1_dxsi2 * dy2_dxsi1; 

			matrix<double> dxsi_dy(2, 2);
			dxsi_dy(0,0)= dy2_dxsi2/jac;
			dxsi_dy(0,1)= -dy1_dxsi2/jac;
			dxsi_dy(1,0)= -dy2_dxsi1/jac;
			dxsi_dy(1,1)= dy1_dxsi1/jac;

			matrix<double> dphi_dy(2, nodes_.size());
			for (size_t i = 0; i < nodes_.size(); i++)
			{
					dphi_dy(0,i)= dphi_dxsi(0,i)*dxsi_dy(0,0)+dphi_dxsi(1,i)*dxsi_dy(1,0);
					dphi_dy(1,i)= dphi_dxsi(0,i)*dxsi_dy(0,1)+dphi_dxsi(1,i)*dxsi_dy(1,1);
			}

			// cálculo das variaveis nos pontos de integração
			bounded_vector<double, 2> rhou_;
			rhou_(0) = 0.;
			rhou_(1) = 0.;
			double rho_ = 0.;
			double rhoe_ = 0.;
			bounded_vector<double, 2> drhoe_dy;
			drhoe_dy.clear();
			bounded_vector<double, 2> drho_dy;
			drho_dy.clear();
			matrix<double> drhou_dy(2, 2);
			drhou_dy.clear();
			double gamma_= (material_ -> getSpecificHeatp())/(material_ -> getSpecificHeatv());
		
		 	for (size_t i = 0; i < nodes_.size(); i++)
			{

				rho_ += nodes_[i]->getCurrentDensity()*phi(i);
				rhou_(0) += nodes_[i]->getCurrentMomentum()(0)*phi(i);
				rhou_(1) += nodes_[i]->getCurrentMomentum()(1)*phi(i);
				rhoe_ += nodes_[i]->getCurrentInternalEnergy()*phi(i);
				for (int k=0;k<2;k++){
					drho_dy(k) += nodes_[i]->getCurrentDensity()*dphi_dy(k,i);
					drhou_dy(0,k) += nodes_[i]->getCurrentMomentum()(0)*dphi_dy(k,i);
					drhou_dy(1,k) += nodes_[i]->getCurrentMomentum()(1)*dphi_dy(k,i);
					drhoe_dy(k) += nodes_[i]->getCurrentInternalEnergy()*dphi_dy(k,i);
				}
		 	} 
			
		//Get viscosity
		double mi_0 = material_->getViscosity(); // initial viscosity
		double mi; // current viscosity
		double tempInf = material_->getUndisturbedTemperature(); //initial temperature
		double T; // current temperature
		double calorv = material_->getSpecificHeatv();
		bounded_vector<double, 2> u_; //velocity
		u_(0) = rhou_(0)/rho_;
		u_(1) = rhou_(1)/rho_;

		//Temperatura
		double e = rhoe_/rho_;
		T=(1./calorv)*(e-(u_(0)*u_(0)+u_(1)*u_(1))/2.);

		if (sutherland==true)
		{
			//mi = mi_0;
			mi = mi_0*pow(((T+273.15)/(tempInf+273.15)),(3./2.))*((tempInf+273.15+110.4)/(T+273.15+110.4));
		}
		
		else if (sutherland==false)
		{
				mi = mi_0;
		}
	
		matrix<double> du_dy(2, 2); //velocity gradient
		matrix<double> tau_(2, 2); // shear stress
		bounded_vector<double, 2> de_dy; //energy gradient
		bounded_vector<double, 2> dT_dy; //temperature gradient

	

		du_dy(0,0) = (drhou_dy(0,0)*rho_-rhou_(0)*drho_dy(0))/(rho_*rho_);
		du_dy(0,1) = (drhou_dy(0,1)*rho_-rhou_(0)*drho_dy(1))/(rho_*rho_);
		du_dy(1,0) = (drhou_dy(1,0)*rho_-rhou_(1)*drho_dy(0))/(rho_*rho_);
		du_dy(1,1) = (drhou_dy(1,1)*rho_-rhou_(1)*drho_dy(1))/(rho_*rho_);
	
		tau_(0,0) = -2./3. * mi * (du_dy(0,0) + du_dy(1,1)) +2. * mi * du_dy(0,0);
		tau_(0,1) = mi*(du_dy(0,1)+du_dy(1,0));
		tau_(1,0) = mi*(du_dy(1,0)+du_dy(0,1));
		tau_(1,1) = -2./3. * mi * (du_dy(0,0) + du_dy(1,1)) +2. * mi * du_dy(1,1);
		
		//Temperatura
		for (size_t i = 0; i < 2; i++){
			de_dy(i)=(drhoe_dy(i)*rho_-drho_dy(i)*rhoe_)/(rho_*rho_);
		}

		for (size_t i = 0; i < 2; i++){
			dT_dy(i)=(1./calorv)*(de_dy(i)-u_(0)*du_dy(0,i)-u_(1)*du_dy(1,i));
		}
		
		// determinação do Kt
		double Kt_ = material_->getThermalCond();

	
		//rhs(2) += deltat_ *phi(2) * (Normal(0)*(tau_(0,0)*u_(0)+tau_(0,1)*u_(1)+ Kt_*dT_dy(0))+Normal(1)*(tau_(1,0)*u_(0)+tau_(1,1)*u_(1)+Kt_*dT_dy(1)))*weight*sideLenght/2.;
		//rhs(0) += deltat_ *phi(0) * (Normal(0)*(tau_(0,0)*u_(0)+tau_(0,1)*u_(1)+ Kt_*dT_dy(0))+Normal(1)*(tau_(1,0)*u_(0)+tau_(1,1)*u_(1)+Kt_*dT_dy(1)))*weight*sideLenght/2.;

	 	} // loop over gauss points
	}
	
	return std::make_pair(rhs, tangent);
}

bounded_vector<double, 6> Element::elementShockLowSpeedS1()
{
	bounded_vector<double, 6> rhs;
	bounded_vector<double, 6> rhsaux;
	bounded_matrix<double, 6, 6> tangent;     // 2 por nó
	rhs.clear(); tangent.clear();
	double jac=0.;

	 for (int ih = 0; ih < numberOfIntegrationPoints_; ih++) // verificar o number
	{
	 	double xsi1 = domainIntegrationPoints_(ih,0);
	 	double xsi2 = domainIntegrationPoints_(ih,1);
	 	double weight = domainIntegrationPoints_(ih,2);

	 	vector<double> phi = domainShapeFunction(xsi1, xsi2);
	 	matrix<double> dphi_dxsi = domainDerivativeShapeFunction(xsi1, xsi2); //row = direction, column = node

		double dy1_dxsi1 = 0.0, dy1_dxsi2 = 0.0, dy2_dxsi1 = 0.0, dy2_dxsi2 = 0.0;
	 	for (size_t i = 0; i < nodes_.size(); i++)
	 	{
	 		bounded_vector<double, 2> currentCoord = nodes_[i]->getCurrentCoordinate();
	 		dy1_dxsi1 += currentCoord(0) * dphi_dxsi(0,i);
	 		dy1_dxsi2 += currentCoord(0) * dphi_dxsi(1,i);
	 		dy2_dxsi1 += currentCoord(1) * dphi_dxsi(0,i);
	 		dy2_dxsi2 += currentCoord(1) * dphi_dxsi(1,i);
	 	}

	 	//current jacobian
	 	jac = dy1_dxsi1 * dy2_dxsi2 - dy1_dxsi2 * dy2_dxsi1; 
		
		for (int i = 0; i < nodes_.size(); i++)
		{
			for (int j = 0; j < nodes_.size(); j++)
			{
			if(i!=j){
				tangent(2*i,2*j) += phi(i) * phi(j) * jac * weight;
				tangent(2*i+1,2*j+1) += phi(i) * phi(j) * jac * weight;
			}

			
			}
		
	 	}
	}
	
	for (int i = 0; i < nodes_.size(); i++)
	{
		rhsaux(2*i) = nodes_[i]->getDeltaMomentum()(0);
		rhsaux(2*i+1) = nodes_[i]->getDeltaMomentum()(1);
	}

	for (int i = 0; i < 2*nodes_.size(); i++)
	{
		for (int j = 0; j < 2*nodes_.size(); j++)
		{				
			rhs(i) += tangent(i,j)* rhsaux(j);
		}
	}

	return rhs;
}

bounded_vector<double, 3> Element::elementShockLowSpeedS2()
{
	bounded_vector<double, 3> rhs;
	bounded_vector<double, 3> rhsaux;
	bounded_matrix<double, 3, 3> tangent;     // 2 por nó
	rhs.clear(); tangent.clear();
	double jac=0.;

	 for (int ih = 0; ih < numberOfIntegrationPoints_; ih++) // verificar o number
	 {
	 	double xsi1 = domainIntegrationPoints_(ih,0);
	 	double xsi2 = domainIntegrationPoints_(ih,1);
	 	double weight = domainIntegrationPoints_(ih,2);

	 	vector<double> phi = domainShapeFunction(xsi1, xsi2);
	 	matrix<double> dphi_dxsi = domainDerivativeShapeFunction(xsi1, xsi2); //row = direction, column = node

		double dy1_dxsi1 = 0.0, dy1_dxsi2 = 0.0, dy2_dxsi1 = 0.0, dy2_dxsi2 = 0.0;
	 	for (size_t i = 0; i < nodes_.size(); i++)
	 	{
	 		bounded_vector<double, 2> currentCoord = nodes_[i]->getCurrentCoordinate();
	 		dy1_dxsi1 += currentCoord(0) * dphi_dxsi(0,i);
	 		dy1_dxsi2 += currentCoord(0) * dphi_dxsi(1,i);
	 		dy2_dxsi1 += currentCoord(1) * dphi_dxsi(0,i);
	 		dy2_dxsi2 += currentCoord(1) * dphi_dxsi(1,i);
	 	}

	 	//current jacobian
	 	jac = dy1_dxsi1 * dy2_dxsi2 - dy1_dxsi2 * dy2_dxsi1; 
		
		for (int i = 0; i < nodes_.size(); i++)
		{
			for (int j = 0; j < nodes_.size(); j++)
			{
				if(i!=j){
				tangent(i,j) += phi(i) * phi(j) * jac * weight;
				}			
			}
		
	 	}
	 }
	
	for (int i = 0; i < nodes_.size(); i++)
	{
		rhsaux(i) = nodes_[i]->getDeltaDensity();
			
	}

	for (int i = 0; i < nodes_.size(); i++)
	{
		for (int j = 0; j < nodes_.size(); j++)
		{
			
			rhs(i) += tangent(i,j)* rhsaux(j);
		}
	}

	return rhs;
}

bounded_vector<double, 3> Element::elementShockLowSpeedS3()
{
	bounded_vector<double, 3> rhs;
	bounded_vector<double, 3> rhsaux;
	bounded_matrix<double, 3, 3> tangent;     // 2 por nó
	rhs.clear(); tangent.clear();
	double jac=0.;

	 for (int ih = 0; ih < numberOfIntegrationPoints_; ih++) // verificar o number
	 {
	 	double xsi1 = domainIntegrationPoints_(ih,0);
	 	double xsi2 = domainIntegrationPoints_(ih,1);
	 	double weight = domainIntegrationPoints_(ih,2);

	 	vector<double> phi = domainShapeFunction(xsi1, xsi2);
	 	matrix<double> dphi_dxsi = domainDerivativeShapeFunction(xsi1, xsi2); //row = direction, column = node

		double dy1_dxsi1 = 0.0, dy1_dxsi2 = 0.0, dy2_dxsi1 = 0.0, dy2_dxsi2 = 0.0;
	 	for (size_t i = 0; i < nodes_.size(); i++)
	 	{
	 		bounded_vector<double, 2> currentCoord = nodes_[i]->getCurrentCoordinate();
	 		dy1_dxsi1 += currentCoord(0) * dphi_dxsi(0,i);
	 		dy1_dxsi2 += currentCoord(0) * dphi_dxsi(1,i);
	 		dy2_dxsi1 += currentCoord(1) * dphi_dxsi(0,i);
	 		dy2_dxsi2 += currentCoord(1) * dphi_dxsi(1,i);
	 	}

	 	//current jacobian
	 	jac = dy1_dxsi1 * dy2_dxsi2 - dy1_dxsi2 * dy2_dxsi1; 
		
		for (int i = 0; i < nodes_.size(); i++)
		{
			for (int j = 0; j < nodes_.size(); j++)
			{
				if(i!=j){
				tangent(i,j) += phi(i) * phi(j) * jac * weight;			
				}
			}
		
	 	}
	 }
	for (int i = 0; i < nodes_.size(); i++)
	{
		rhsaux(i) = nodes_[i]->getDeltaEnergy();
		
	}

	for (int i = 0; i < nodes_.size(); i++)
	{
		for (int j = 0; j < nodes_.size(); j++)
		{
				
			rhs(i) += tangent(i,j)* rhsaux(j);
		}
	}

	return rhs;
}

	bounded_matrix<double, 3, 3> Element:: elementContributionsMesh(){
	bounded_matrix<double, 3, 3> tangent;     // 2 por nó
	tangent.clear();
	
	double jac=0.;

	
	 
	for (int ih = 0; ih < numberOfIntegrationPoints_; ih++) // verificar o number
	{
	 	double xsi1 = domainIntegrationPoints_(ih,0);
	 	double xsi2 = domainIntegrationPoints_(ih,1);
	 	double weight = domainIntegrationPoints_(ih,2);

		vector<double> phi = domainShapeFunction(xsi1, xsi2);
	 	matrix<double> dphi_dxsi = domainDerivativeShapeFunction(xsi1, xsi2); //row = direction, column = node

	 	//dy_dxsi at time t+1 = Aint
	 	double dy1_dxsi1 = 0.0, dy1_dxsi2 = 0.0, dy2_dxsi1 = 0.0, dy2_dxsi2 = 0.0;
	 	for (size_t i = 0; i < nodes_.size(); i++)
	 	{
	 		bounded_vector<double, 2> currentCoord = nodes_[i]->getCurrentCoordinate();
	 		dy1_dxsi1 += currentCoord(0) * dphi_dxsi(0,i);
	 		dy1_dxsi2 += currentCoord(0) * dphi_dxsi(1,i);
	 		dy2_dxsi1 += currentCoord(1) * dphi_dxsi(0,i);
	 		dy2_dxsi2 += currentCoord(1) * dphi_dxsi(1,i);
	 	}

	 	//current jacobian
	 	jac = dy1_dxsi1 * dy2_dxsi2 - dy1_dxsi2 * dy2_dxsi1; 

		matrix<double> dxsi_dy(2, 2);
		dxsi_dy(0,0)= dy2_dxsi2/jac;
		dxsi_dy(0,1)= -dy1_dxsi2/jac;
		dxsi_dy(1,0)= -dy2_dxsi1/jac;
		dxsi_dy(1,1)= dy1_dxsi1/jac;

		matrix<double> dphi_dy(2, nodes_.size());
		for (size_t i = 0; i < nodes_.size(); i++)
		{
			dphi_dy(0,i)= dphi_dxsi(0,i)*dxsi_dy(0,0)+dphi_dxsi(1,i)*dxsi_dy(1,0);
			dphi_dy(1,i)= dphi_dxsi(0,i)*dxsi_dy(0,1)+dphi_dxsi(1,i)*dxsi_dy(1,1);
		}

	    for (int i = 0; i < 3; i++){
        	for (int j = 0; j < 3; j++){        
            	tangent(i  ,j  ) += (dphi_dy(0,i) * dphi_dy(0,j) +
                                              dphi_dy(1,i) * dphi_dy(1,j)) * weight * jac * meshStiffness_;
	        };
    	};

	}

	return tangent;

}
