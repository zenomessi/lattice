// system_general_cpp.cpp
//  
//
//  Created by Zeno Messi on 28.01.19
//

#include<gsl/gsl_math.h>
#include "system_general_cpp.hpp"
#include<limits>

double epsilon_s(std::numeric_limits<double>::epsilon());

System_General::System_General(const Parseur& par): Lattice(par), spring_coef(par.get(Spring)), bend_coef(par.get(Bend)){
	
	this->param=new Parseur(par);
	cerr << "A new System_General was created" << endl;
}

double System_General::get_state(const gsl_vector * var, const unsigned int& i, unsigned int j) const{
	
	if(j==0) return gsl_vector_get(var, i*2)+this->points[i].first;
	return gsl_vector_get(var, i*2+1)+this->points[i].second;
}

double System_General::get_var(const gsl_vector * var, const unsigned int& i, unsigned int j) const{
	
	if(j==0) return gsl_vector_get(var, i*2);
	return gsl_vector_get(var, i*2+1);
}

double System_General::get_spring_E(const gsl_vector * var) const{//one call for the whole system
	//a more human intelligible version can be found in the older version system_general.cpp
	double res(0.);
	
	for(unsigned int i(0); i<this->lattice_hinges_size; i++){
		for(unsigned int j(0); j<2; j++){
			double v_x(this->get_var(var, this->get_hinge_point(i,j+1),0)-this->get_var(var, this->get_hinge_point(i,j), 0));
			double v_y(this->get_var(var, this->get_hinge_point(i,j+1),1)-this->get_var(var, this->get_hinge_point(i,j), 1));
			double r_x(this->get_point(this->get_hinge_point(i,j+1), 0)-this->get_point(this->get_hinge_point(i,j), 0));
			double r_y(this->get_point(this->get_hinge_point(i,j+1), 1)-this->get_point(this->get_hinge_point(i,j), 1));
			
			//res+= pow((v_x*(2.*r_x+v_x)+v_y*(2.*r_y+v_y))/(sqrt((r_x+v_x)*(r_x+v_x)+(r_y+v_y)*(r_y+v_y))+1.),2.);//rest_length==1.
			res+= pow((v_x*(2.*r_x+v_x)+v_y*(2.*r_y+v_y))/(sqrt((r_x+v_x)*(r_x+v_x)+(r_y+v_y)*(r_y+v_y))+this->rest_length),2.);//rest_length variable
		}
	}
	//cerr << "Le problème n'est pas dans spring" << endl;
	//cerr << "Spring E " << res*spring_coef/2. << endl;
	return res*this->spring_coef/2.;
}

//double System_General::get_spring_dE_x(const gsl_vector* var, const unsigned int & i) const{
	////for x axis
	//double res(0.);

	//for(unsigned int j(0); j<this->hinges_size; j++){
		//for(unsigned int k(0); k<3; k+=2){
			//if(this->get_hinge_point(j,k)==i){
				
				//double v_x(this->get_var(var, this->get_hinge_point(j,1),0)-this->get_var(var, i, 0));
				//double v_y(this->get_var(var, this->get_hinge_point(j,1),1)-this->get_var(var, i, 1));
				//double r_x(this->get_point(this->get_hinge_point(j,1), 0)-this->get_point(i, 0));
				//double r_y(this->get_point(this->get_hinge_point(j,1), 1)-this->get_point(i, 1));
				
				//res+= (v_x*(2.*r_x+v_x)+v_y*(2.*r_y+v_y))/(sqrt((r_x+v_x)*(r_x+v_x)+(r_y+v_y)*(r_y+v_y))+1.)*(v_x+r_x)/sqrt((r_x+v_x)*(r_x+v_x)+(r_y+v_y)*(r_y+v_y));
				////res+= (v_x*(2.*r_x+v_x)+v_y*(2.*r_y+v_y))/(sqrt((r_x+v_x)*(r_x+v_x)+(r_y+v_y)*(r_y+v_y))+this->rest_length)*(v_x*r_x)/((r_x+v_x)*(r_x+v_x)+(r_y+v_y)*(r_y+v_y))
			//}
		//}
	//}
	//return res*this->spring_coef;
//}
//double System_General::get_spring_dE_y(const gsl_vector* var, const unsigned int & i) const{
	////for y axis
	//double res(0.);
	
	//for(unsigned int j(0); j<this->hinges_size; j++){
		//for(unsigned int k(0); k<3; k+=2){
			//if(this->get_hinge_point(j,k)==i){
				
				//double v_x(this->get_var(var, this->get_hinge_point(j,1),0)-this->get_var(var, i, 0));
				//double v_y(this->get_var(var, this->get_hinge_point(j,1),1)-this->get_var(var, i, 1));
				//double r_x(this->get_point(this->get_hinge_point(j,1), 0)-this->get_point(i, 0));
				//double r_y(this->get_point(this->get_hinge_point(j,1), 1)-this->get_point(i, 1));
				
				//res+= (v_x*(2.*r_x+v_x)+v_y*(2.*r_y+v_y))/(sqrt((r_x+v_x)*(r_x+v_x)+(r_y+v_y)*(r_y+v_y))+1.)*(v_y+r_y)/sqrt((r_x+v_x)*(r_x+v_x)+(r_y+v_y)*(r_y+v_y));
				////res+= (v_x*(2.*r_x+v_x)+v_y*(2.*r_y+v_y))/(sqrt((r_x+v_x)*(r_x+v_x)+(r_y+v_y)*(r_y+v_y))+this->rest_length)*(v_y*r_y)/((r_x+v_x)*(r_x+v_x)+(r_y+v_y)*(r_y+v_y))
			//}
		//}
	//}
	//return res*this->spring_coef;
//}

double System_General::get_spring_dE_x(const gsl_vector* var, const unsigned int & i) const{
	
	double res(0.);
	for(auto it(edges_first[i].begin()); it!=edges_first[i].end(); ++it){
		double r_x(-this->points[i].first+this->points[*it].first);
		double r_y(-this->points[i].second+this->points[*it].second);
		double v_x(-this->get_var(var, i,0)+this->get_var(var, *it,0));
		double v_y(-this->get_var(var, i,1)+this->get_var(var, *it,1));
		
		//res-= (v_x*(2.*r_x+v_x)+v_y*(2.*r_y+v_y))/(sqrt((r_x+v_x)*(r_x+v_x)+(r_y+v_y)*(r_y+v_y))+1.)*(v_x+r_x)/sqrt((r_x+v_x)*(r_x+v_x)+(r_y+v_y)*(r_y+v_y));//rest_length==1.
		res-= (v_x*(2.*r_x+v_x)+v_y*(2.*r_y+v_y))/(sqrt((r_x+v_x)*(r_x+v_x)+(r_y+v_y)*(r_y+v_y))+this->rest_length)*(v_x+r_x)/sqrt((r_x+v_x)*(r_x+v_x)+(r_y+v_y)*(r_y+v_y));//rest_length variable
	}
	for(auto it(edges_second[i].begin()); it!=edges_second[i].end(); ++it){
		double r_x(this->points[i].first-this->points[*it].first);
		double r_y(this->points[i].second-this->points[*it].second);
		double v_x(this->get_var(var, i,0)-this->get_var(var, *it,0));
		double v_y(this->get_var(var, i,1)-this->get_var(var, *it,1));
		
		//res+= (v_x*(2.*r_x+v_x)+v_y*(2.*r_y+v_y))/(sqrt((r_x+v_x)*(r_x+v_x)+(r_y+v_y)*(r_y+v_y))+1.)*(v_x+r_x)/sqrt((r_x+v_x)*(r_x+v_x)+(r_y+v_y)*(r_y+v_y));//rest_length==1.
		res+= (v_x*(2.*r_x+v_x)+v_y*(2.*r_y+v_y))/(sqrt((r_x+v_x)*(r_x+v_x)+(r_y+v_y)*(r_y+v_y))+this->rest_length)*(v_x+r_x)/sqrt((r_x+v_x)*(r_x+v_x)+(r_y+v_y)*(r_y+v_y));//rest_length variable
	}
	return res*this->spring_coef;
}

double System_General::get_spring_dE_y(const gsl_vector* var, const unsigned int & i) const{
	
	double res(0.);
	for(auto it(edges_first[i].begin()); it!=edges_first[i].end(); ++it){
		double r_x(-this->points[i].first+this->points[*it].first);
		double r_y(-this->points[i].second+this->points[*it].second);
		double v_x(-this->get_var(var, i,0)+this->get_var(var, *it,0));
		double v_y(-this->get_var(var, i,1)+this->get_var(var, *it,1));
		
		//res-= (v_x*(2.*r_x+v_x)+v_y*(2.*r_y+v_y))/(sqrt((r_x+v_x)*(r_x+v_x)+(r_y+v_y)*(r_y+v_y))+1.)*(v_y+r_y)/sqrt((r_x+v_x)*(r_x+v_x)+(r_y+v_y)*(r_y+v_y));//rest_length==1.
		res-= (v_x*(2.*r_x+v_x)+v_y*(2.*r_y+v_y))/(sqrt((r_x+v_x)*(r_x+v_x)+(r_y+v_y)*(r_y+v_y))+this->rest_length)*(v_y+r_y)/sqrt((r_x+v_x)*(r_x+v_x)+(r_y+v_y)*(r_y+v_y));//rest_length variable
	}
	for(auto it(edges_second[i].begin()); it!=edges_second[i].end(); ++it){
		double r_x(this->points[i].first-this->points[*it].first);
		double r_y(this->points[i].second-this->points[*it].second);
		double v_x(this->get_var(var, i,0)-this->get_var(var, *it,0));
		double v_y(this->get_var(var, i,1)-this->get_var(var, *it,1));
		
		//res+= (v_x*(2.*r_x+v_x)+v_y*(2.*r_y+v_y))/(sqrt((r_x+v_x)*(r_x+v_x)+(r_y+v_y)*(r_y+v_y))+1.)*(v_y+r_y)/sqrt((r_x+v_x)*(r_x+v_x)+(r_y+v_y)*(r_y+v_y));//rest_length==1.
		res+= (v_x*(2.*r_x+v_x)+v_y*(2.*r_y+v_y))/(sqrt((r_x+v_x)*(r_x+v_x)+(r_y+v_y)*(r_y+v_y))+this->rest_length)*(v_y+r_y)/sqrt((r_x+v_x)*(r_x+v_x)+(r_y+v_y)*(r_y+v_y));//rest_length variable
	}
	return res*this->spring_coef;
}

double System_General::get_bend_E(const gsl_vector* var) const{
	
	double res(0.);
	for(auto it(this->hinges.begin()); it!=this->hinges.end(); ++it){
	//for(unsigned int i(0); i<){
		
		//double r_x(this->get_point(this->get_hinge_point(
		double r_x(this->points[get<1>(*it)].first-this->points[get<0>(*it)].first);
		double r_y(this->points[get<1>(*it)].second-this->points[get<0>(*it)].second);
		double v_1_x(this->get_var(var, get<1>(*it),0)-this->get_var(var, get<0>(*it),0));
		double v_2_x(this->get_var(var, get<2>(*it),0)-this->get_var(var, get<1>(*it),0));
		double v_1_y(this->get_var(var, get<1>(*it),1)-this->get_var(var, get<0>(*it),1));
		double v_2_y(this->get_var(var, get<2>(*it),1)-this->get_var(var, get<1>(*it),1));
		
		double norm_u_1(sqrt((v_1_x+r_x)*(v_1_x+r_x)+(v_1_y+r_y)*(v_1_y+r_y)));
		double norm_u_2(sqrt((v_2_x+r_x)*(v_2_x+r_x)+(v_2_y+r_y)*(v_2_y+r_y)));
		
		if((r_x+v_1_x)*(r_x+v_2_x)+(r_y+v_1_y)*(r_y+v_2_y)+norm_u_1*norm_u_2<0.1) res+= 1- ((r_x+v_1_x)*(r_x+v_2_x)+(r_y+v_1_y)*(r_y+v_2_y))/(norm_u_1*norm_u_2);
		else res+= pow((v_1_x-v_2_x)*r_y-(v_1_y-v_2_y)*r_x+v_1_x*v_2_y-v_1_y*v_2_x,2.)/((r_x+v_1_x)*(r_x+v_2_x)+(r_y+v_1_y)*(r_y+v_2_y)+norm_u_1*norm_u_2)/(norm_u_1*norm_u_2);
	}
	//cerr << "Le problème n'est pas dans bend" << endl;
	//cerr << "Bend E " << res*this->bend_coef << endl;
	return res*this->bend_coef;
}

double System_General::get_bend_dE_x(const gsl_vector* var, const unsigned int& i) const{
	
	double res(0.);
	
	//first loop on the hinges where it's the first node
	for(auto it(this->hinges_first[i].begin()); it!= hinges_first[i].end(); ++it){
		
		double r_x(this->points[it->first].first-this->points[i].first);
		double r_y(this->points[it->first].second-this->points[i].second);
		double v_1_x(this->get_var(var, it->first,0)-this->get_var(var, i,0));
		double v_1_y(this->get_var(var, it->first,1)-this->get_var(var, i,1));
		double v_2_x(this->get_var(var, it->second,0)-this->get_var(var, it->first,0));
		double v_2_y(this->get_var(var, it->second,1)-this->get_var(var, it->first,1));
		
		res-=this->compute_bend_dE_x(r_x,r_y,v_1_x,v_1_y,v_2_x,v_2_y);
	}
	
	//then where it's the second...
	for(auto it(this->hinges_mid[i].begin()); it!= hinges_mid[i].end(); ++it){
		
		double r_x(this->points[i].first-this->points[it->first].first);
		double r_y(this->points[i].second-this->points[it->first].second);
		double v_1_x(this->get_var(var, i,0)-this->get_var(var, it->first,0));
		double v_1_y(this->get_var(var, i,1)-this->get_var(var, it->first,1));
		double v_2_x(this->get_var(var, it->second,0)-this->get_var(var, i,0));
		double v_2_y(this->get_var(var, it->second,1)-this->get_var(var, i,1));
		
		res+=this->compute_bend_dE_x(r_x,r_y,v_1_x,v_1_y,v_2_x,v_2_y);
		res-=this->compute_bend_dE_x(r_x,r_y,v_2_x,v_2_y,v_1_x,v_1_y);
	}
	
	//And last but not least, where it's last
	for(auto it(this->hinges_last[i].begin()); it!= hinges_last[i].end(); ++it){
		
		double r_x(this->points[it->second].first-this->points[it->first].first);
		double r_y(this->points[it->second].second-this->points[it->first].second);
		double v_1_x(this->get_var(var, it->second,0)-this->get_var(var, it->first,0));
		double v_1_y(this->get_var(var, it->second,1)-this->get_var(var, it->first,1));
		double v_2_x(this->get_var(var, i,0)-this->get_var(var, it->second,0));
		double v_2_y(this->get_var(var, i,1)-this->get_var(var, it->second,1));
		
		res+=this->compute_bend_dE_x(r_x,r_y,v_2_x,v_2_y,v_1_x,v_1_y);
		//cerr << this->compute_bend_dE_x(r_x,r_y,v_2_x,v_2_y,v_1_x,v_1_y) << endl;
	}
	return res*bend_coef;
}

double System_General::get_bend_dE_y(const gsl_vector* var, const unsigned int& i) const{
	
	double res(0.);
	
	//first loop on the hinges where it's the first node
	for(auto it(this->hinges_first[i].begin()); it!= hinges_first[i].end(); ++it){
		
		double r_x(this->points[it->first].first-this->points[i].first);
		double r_y(this->points[it->first].second-this->points[i].second);
		double v_1_x(this->get_var(var, it->first,0)-this->get_var(var, i,0));
		double v_1_y(this->get_var(var, it->first,1)-this->get_var(var, i,1));
		double v_2_x(this->get_var(var, it->second,0)-this->get_var(var, it->first,0));
		double v_2_y(this->get_var(var, it->second,1)-this->get_var(var, it->first,1));
		
		res-=this->compute_bend_dE_y(r_x,r_y,v_1_x,v_1_y,v_2_x,v_2_y);
	}
	
	//then where it's the second...
	for(auto it(this->hinges_mid[i].begin()); it!= hinges_mid[i].end(); ++it){
		
		double r_x(this->points[i].first-this->points[it->first].first);
		double r_y(this->points[i].second-this->points[it->first].second);
		double v_1_x(this->get_var(var, i,0)-this->get_var(var, it->first,0));
		double v_1_y(this->get_var(var, i,1)-this->get_var(var, it->first,1));
		double v_2_x(this->get_var(var, it->second,0)-this->get_var(var, i,0));
		double v_2_y(this->get_var(var, it->second,1)-this->get_var(var, i,1));
		
		res+=this->compute_bend_dE_y(r_x,r_y,v_1_x,v_1_y,v_2_x,v_2_y);
		res-=this->compute_bend_dE_y(r_x,r_y,v_2_x,v_2_y,v_1_x,v_1_y);
	}
	
	//And last but not least, where it's last
	for(auto it(this->hinges_last[i].begin()); it!= hinges_last[i].end(); ++it){
		
		double r_x(this->points[it->second].first-this->points[it->first].first);
		double r_y(this->points[it->second].second-this->points[it->first].second);
		double v_1_x(this->get_var(var, it->second,0)-this->get_var(var, it->first,0));
		double v_1_y(this->get_var(var, it->second,1)-this->get_var(var, it->first,1));
		double v_2_x(this->get_var(var, i,0)-this->get_var(var, it->second,0));
		double v_2_y(this->get_var(var, i,1)-this->get_var(var, it->second,1));
		
		res+=this->compute_bend_dE_y(r_x,r_y,v_2_x,v_2_y,v_1_x,v_1_y);
	}
	
	return res*bend_coef;
}

double System_General::compute_bend_dE_x(const double& r_x, const double& r_y, const double& v_1_x, const double& v_1_y, const double& v_2_x, const double v_2_y) const{
    
	double norm_u_1_s((v_1_x+r_x)*(v_1_x+r_x)+(v_1_y+r_y)*(v_1_y+r_y));
	double norm_u_2_s((v_2_x+r_x)*(v_2_x+r_x)+(v_2_y+r_y)*(v_2_y+r_y));

	double v_2_v_1_u_1((v_2_x-v_1_x)*(r_x+v_1_x)+(v_2_y-v_1_y)*(r_y+v_1_y));//LHT of the numerator
	double u_1_u_2((r_x+v_1_x)*(r_x+v_2_x)+(r_y+v_1_y)*(r_y+v_2_y));//RHT of the numerator
	
	return ((r_x*v_2_v_1_u_1 + v_1_x*u_1_u_2)/norm_u_1_s-v_2_x)/(sqrt(norm_u_1_s*norm_u_2_s));
}

double System_General::compute_bend_dE_y(const double& r_x, const double& r_y, const double& v_1_x, const double& v_1_y, const double& v_2_x, const double v_2_y) const{
    
	double norm_u_1_s((v_1_x+r_x)*(v_1_x+r_x)+(v_1_y+r_y)*(v_1_y+r_y));
	double norm_u_2_s((v_2_x+r_x)*(v_2_x+r_x)+(v_2_y+r_y)*(v_2_y+r_y));

	double v_2_v_1_u_1((v_2_x-v_1_x)*(r_x+v_1_x)+(v_2_y-v_1_y)*(r_y+v_1_y));//LHT of the numerator
	double u_1_u_2((r_x+v_1_x)*(r_x+v_2_x)+(r_y+v_1_y)*(r_y+v_2_y));//RHT of the numerator
	
	return ((r_y*v_2_v_1_u_1 + v_1_y*u_1_u_2)/norm_u_1_s-v_2_y)/(sqrt(norm_u_1_s*norm_u_2_s));
}

System_General::~System_General(){
	
	cerr << "System_General desctructor entered" << endl;
	delete this->param;
	cerr << "The System_General instance has been correctly destroyed" << endl;
}
