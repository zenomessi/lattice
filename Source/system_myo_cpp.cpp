//
//  system_myo_cpp.cpp
//  
//	
//  Created by Zeno Messi on 29.01.19
//

#include<gsl/gsl_math.h>
#include "system_myo_cpp.hpp"

System::System(Parseur& par): System_General(par), F_magn(par.get(FMagnitude)), dip_prob(par.get(Myo)), min_dist(0.0001){
	
	cerr << "creating the System instance" << endl;
	unsigned seed2 = std::chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator(seed2);
	//default_random_engine generator(3582721925);
	
	//default_random_engine generator(par.get(SeedMyo));
	//cerr << "Seedmyo " << par.get(SeedMyo) << endl;
	//cerr << "Seedlattice " << par.get_seed(SeedLat) << endl;
	this->param->set_seed(SeedMyo, seed2);
	//cerr << "Seedmyo " << seed2 << endl;
	//cerr << "Seedlattice " << this->par.get_seed(SeedLat) << endl;
	uniform_real_distribution<double> distribution(0.0,1.0);
	
	//create the list of dipoles
	//version with dipoles also without the link
	for(unsigned int i(0); i<latsize; i++){
		for(unsigned int j(i); j<latsize; j++){
			//if(this->points[i].first>0) continue;//only dipoles on the left side of the cell
			
			if(sqrt(pow(this->points[i].first-this->points[j].first,2.)+pow(this->points[i].second-this->points[j].second,2.))<2.5*this->rest_length){
				
				if(distribution(generator)<this->dip_prob) this->dip_list.push_back(make_pair(i,j));
			}
		}
	}
	
	//version with dipoles only on links
	//for(unsigned int i(0); i<this->lattice_hinges_size; i++){
		
		//if(distribution(generator)<this->dip_prob) this->dip_list.push_back(make_pair(get<0>(this->hinges[i]), get<2>(this->hinges[i])));
	//}
	
	//create the helping vector for dipoles
	this->dip_list_other=vector<vector<unsigned int> >(this->latsize);
	for(auto it(this->dip_list.begin()); it!=this->dip_list.end(); ++it){
		dip_list_other[it->first].push_back(it->second);
		dip_list_other[it->second].push_back(it->first);
	}
}

//this function is made to create the dipoles from a file of densities of the former final stage
void System::set_dipoles(const string& filepath){
	
	cerr << "creating dipoles from file " << filepath << endl;
	unsigned seed2 = std::chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator(seed2);
	
	uniform_real_distribution<double> distribution(0.0,1.0);
	
	vector<double> dip_density;
	ifstream densfile(filepath+"_dipoles_density.dat", ios::in);
	unsigned int side_out, size_out, max_distance_ui, dip_count;
	densfile.ignore(3,' ');
	densfile >> side_out >> size_out >> max_distance_ui >> dip_count;
	cerr << "size of a square of the grid in pixels " << side_out << endl;
	cerr << "size of the grid side in squares " << size_out << endl;
	cerr << "maximal distance to the origin " << max_distance_ui << endl;
	cerr << "number of dipoles " << dip_count << endl;
	double temp;
	densfile >> temp;
	while(!densfile.eof()){
		dip_density.push_back(temp);
		densfile >> temp;
	}
	
	this->dip_list.clear();
	
	Cell last_final_contour(filepath+"_final_outline.dat");
	last_final_contour.setcharac();
	
	//for the version with only dipoles on one side of the cell at the beginning, we change the dipoles density on the go
	//this->dip_prob=0.1;
	
	uniform_int_distribution<unsigned int> points_distrib(0, this->latsize-1);
	while(this->dip_list.size()<dip_count){
		
		unsigned int i(points_distrib(generator));
		unsigned int j(points_distrib(generator));
		if(i==j) continue;
		bool flag_deja(false);
		for(auto dipole : this->dip_list){
			if((dipole.first==i && dipole.second ==j) or (dipole.first==j && dipole.second==i)){
				flag_deja=true;
				break;
			}
		}
		if(flag_deja) continue;
		
	//create the list of dipoles
	//version with dipoles also without the link
	//for(unsigned int i(0); i<latsize; i++){
		//for(unsigned int j(i); j<latsize; j++){
			
			if(sqrt(pow(this->points[i].first-this->points[j].first,2.)+pow(this->points[i].second-this->points[j].second,2.))<2.5*this->rest_length){//version with the distance between the two points
			//int index_i(this->indices[j].first-this->indices[i].first);
			//int index_j(this->indices[j].second-this->indices[i].second);
			//if((index_i==2 and index_j==0) or (index_i==0 and index_j==2) or (index_i==2 and index_j==-2)//positive direction
			//or (index_i==-2 and index_j==0) or (index_i==0 and index_j==-2) or (index_i==-2 and index_j==2)){//negative direction
				
				//the middle of the dipole defines its position
				Point dip_mid(make_pair((this->points[i].first+this->points[j].first)/2., (this->points[i].second+this->points[j].second)/2.));
				//shift and scale this position
				Point dip_mid_shifted_scaled(make_pair((dip_mid.first+max_distance_ui)/side_out, (dip_mid.second+max_distance_ui)/side_out));
				//test if the range is ok
				if(last_final_contour.incell(dip_mid.first,dip_mid.second)){
					//if the range is ok use the dipoles density
					if(distribution(generator)<dip_density[(unsigned int)floor(dip_mid_shifted_scaled.first)+(unsigned int)floor(dip_mid_shifted_scaled.second)*size_out]) this->dip_list.push_back(make_pair(i,j));
				}
				//else use the normal method
				else{
					//cerr << "using default dipole generation" << endl;
					if(distribution(generator)<this->dip_prob) this->dip_list.push_back(make_pair(i,j));
				}
			}
		//}
	//}
	}
	cerr << "Previous stage #dipoles " << dip_count << " New stage " << this->dip_list.size() << endl;
	cerr << "Out of dipoles creation" << endl;
	//version with dipoles only on links
	//for(unsigned int i(0); i<this->lattice_hinges_size; i++){
		
		//if(distribution(generator)<this->dip_prob) this->dip_list.push_back(make_pair(get<0>(this->hinges[i]), get<2>(this->hinges[i])));
	//}
	
	//create the helping vector for dipoles
	this->dip_list_other=vector<vector<unsigned int> >(this->latsize);
	cerr << "Latsize "<< latsize << endl;
	cerr << "Entering the creation of helper vector" << endl;
	for(auto it(this->dip_list.begin()); it!=this->dip_list.end(); ++it){
		//cerr << it->first << " " << it->second << endl;
		dip_list_other[it->first].push_back(it->second);
		dip_list_other[it->second].push_back(it->first);
	}
	cerr << "Out of dipoles helper creation" << endl;
}

//System(Parseur& par, const unsigned int & i): F_magn(par.get(FMagnitude)), dip_prob(par.get(Myo)), min_dist(0.01){
	
	//this->param=new Parseur(par);
	//this->rest_length=1.;
	//this->link_prob=1.;
	
	
//}

unsigned int System::get_dip_list_size() const{
	return this->dip_list.size();
}

double System::get_force_E(const gsl_vector* var) const{
	
	double res(0.);
	double min(this->min_dist*2.*this->rest_length);//minimal distance for full force //variable rest_length
	for(auto it(this->dip_list.begin()); it!=dip_list.end(); ++it){
		
		double distance(sqrt(pow((this->points[it->first].first+this->get_var(var,it->first,0))-(this->points[it->second].first+this->get_var(var,it->second,0)),2.)
		+pow((this->points[it->first].second+this->get_var(var, it->first,1))-(this->points[it->second].second+this->get_var(var,it->second,1)),2.)));
		if(distance>min) res+= (distance-min) +min/4.;
		else{
			if(distance>min/2.) res+= distance*distance/min-distance+min/4.;
		}
	}
	//cerr << "Le problème n'est pas dans Force" << endl;
	//cerr << "Force E " << res*F_magn << endl;
	return res*this->F_magn;
}

double System::get_force_dE_x(const gsl_vector* var, const unsigned int& index) const{
	//cerr<< "entered Force_dE_x" << endl;
	double res(0.);
	double min(this->min_dist*2.*this->rest_length);//variable rest_length
	for(auto it(this->dip_list_other[index].begin()); it!=this->dip_list_other[index].end(); ++it){
		
		//double dir_vec(this->points[*it].first+this->get_var(var, *it,0)-this->points[index].first-this->get_var(var,index,0));
		double dir_vec(this->points[index].first+this->get_var(var,index,0)-this->points[*it].first-this->get_var(var, *it,0));
		
		double distance(sqrt(pow((this->points[index].first+this->get_var(var,index,0))-(this->points[*it].first+this->get_var(var,*it,0)),2.)
		+pow(this->points[index].second+this->get_var(var, index,1)-(this->points[*it].second+this->get_var(var,*it,1)),2.)));
		
		if(distance>min) res+=dir_vec/distance;
		//else if(distance>min/2.) res+=(2.*distance/min-1.)*dir_vec/distance;
		else if(distance>min/2.) res+=2.*dir_vec/min-dir_vec/distance;
	}
	//cerr << "Le problème n'est pas dans Force_dEx" << endl;
	return res*this->F_magn;
}

double System::get_force_dE_y(const gsl_vector* var, const unsigned int& index) const{
	//cerr<< "entered Force_dE_y" << endl;
	double res(0.);
	double min(this->min_dist*2.*this->rest_length);//variable rest_length
	for(auto it(this->dip_list_other[index].begin()); it!=this->dip_list_other[index].end(); ++it){
		
		//double dir_vec(this->points[*it].second+this->get_var(var, *it,1)-this->points[index].second-this->get_var(var,index,1));
		double dir_vec(this->points[index].second+this->get_var(var,index,1)-this->points[*it].second-this->get_var(var, *it,1));
		
		double distance(sqrt(pow((this->points[index].first+this->get_var(var,index,0))-(this->points[*it].first+this->get_var(var,*it,0)),2.)
		+pow(this->points[index].second+this->get_var(var, index,1)-(this->points[*it].second+this->get_var(var,*it,1)),2.)));
		
		if(distance>min) res+=dir_vec/distance;
		//else if(distance>min/2.) res+=(2.*distance/min-1.)*dir_vec/distance;
		else if(distance>min/2.) res+=2.*dir_vec/min-dir_vec/distance;
	}
	//cerr << "Le problème n'est pas dans Force_dEy" << endl;
	return res*this->F_magn;
}

void System::plot_dipoles_density(const gsl_vector* var, const unsigned int& side_out, const string& s) const{
	
	//there are two ways the density can be normalized.
	//One is with the maximal number of dipoles that could be present(count all nodes ,
	//the other with the maximal number of links in the simulation
	
	//side_out is the size of a square of the grid in which we average the dipoles density
	//size_out is the number of squares in a side of the grid
	//The grid corners are 
	//(dl, dr, ur, ul)-> (-max_distance_ui, -max_distance_ui), (max_distance_ui, -max_distance_ui), (max_distance_ui, max_distance_ui), (-max_distance_ui, max_distance_ui)
	
	Point furthest(*max_element(this->points.begin(), this->points.begin()+this->latsize,
	[](Point a, Point b){
		return sqrt(pow(a.first,2.)+pow(a.second,2.)) < sqrt(pow(b.first,2.)+pow(b.second,2.));
	}));
	
	double max_distance(sqrt(pow(furthest.first,2.)+pow(furthest.second, 2.)));
	int max_distance_ui(ceil(max_distance));//the square of densities will be 2*max_distance_uiX2*max_distance_ui
	max_distance_ui += side_out-max_distance_ui%side_out;
	unsigned int size_out(2*max_distance_ui/side_out+3*side_out);
	
	vector<double> neigh_nb(size_out*size_out,0.);//for the normalization with max possible number of dipoles
	vector<double> dipole_nb(size_out*size_out,0.);//number of dipoles 
	
	for(auto it(this->dip_list.begin()); it!=this->dip_list.end(); ++it){
		
		//dipole midpoint defines the position of the dipole
		Point dip_mid(make_pair(this->points[it->second].first+gsl_vector_get(var, it->second*2)+(this->points[it->first].first+gsl_vector_get(var, it->first*2)),
		 this->points[it->second].second+gsl_vector_get(var, it->second*2+1)+(this->points[it->first].second+gsl_vector_get(var, it->first*2+1))));
		dip_mid.first/=2.;
		dip_mid.second/=2.;
		
		//shift and scale the position
		Point dip_mid_shifted_scaled(make_pair((dip_mid.first+max_distance_ui)/side_out, (dip_mid.second+max_distance_ui)/side_out));
		//increment the good position in the vector
		dipole_nb[(unsigned int)floor(dip_mid_shifted_scaled.first)+(unsigned int)floor(dip_mid_shifted_scaled.second)*size_out]++;
	}
	
	//find the maximal number of dipoles
	double max_dip(*max_element(dipole_nb.begin(), dipole_nb.end()));
	cerr << "Maximal number of dipoles in a square " << max_dip << endl;
	//count the number of possible dipoles
	for(auto it(this->points.begin()); it!=this->points.begin()+this->latsize-1; ++it){
		for(auto itt(it+1); itt!=this->points.begin()+this->latsize; ++itt){
			if(sqrt(pow(it->first-itt->first,2.)+pow(it->second-itt->second,2.))<2.5*this->rest_length){
				
				Point bond_mid(make_pair((it->first+gsl_vector_get(var, (it-points.begin())*2)+itt->first+gsl_vector_get(var, (itt-points.begin())*2))/2.,
				 (it->second+gsl_vector_get(var, (it-points.begin())*2+1)+itt->second+gsl_vector_get(var, (itt-points.begin())*2+1))/2.));
				Point bond_mid_shifted_scaled(make_pair((bond_mid.first+max_distance_ui)/side_out, (bond_mid.second+max_distance_ui)/side_out));
				neigh_nb[(unsigned int)floor(bond_mid_shifted_scaled.first)+(unsigned int)floor(bond_mid_shifted_scaled.second)*size_out]++;
			}
		}
	}
	//MAX NUMBER OF POSSIBLE DIPOLES NORMALIZATION
	//for(unsigned int i(0); i<neigh_nb.size(); ++i){
		
		//if(neigh_nb[i]!=0) dipole_nb[i]/=neigh_nb[i];
	//}
	
	//MAX NUMBER OF DIPOLES IN THE FRAME NORMALIZATION
	for(auto& value: dipole_nb){
		value/=max_dip;
	}
	
	
	ofstream ofile(s, ios::out);
	//format that is used for the simulations
	//the first line gives the size of a square, the number of squares in a line, the shift and the number of dipoles
	ofile << "# "<< side_out << " " << size_out << " " << max_distance_ui << " " << this->dip_list.size() <<  endl;
	for(auto it(dipole_nb.begin()); it!=dipole_nb.end(); ++it){//one column
		ofile << *it <<endl;
	}
	//Test format to plot a text image
	//for(unsigned int j(size_out-1); j>0; --j){
		//for(unsigned int i(0); i<size_out; ++i){
			//ofile << dipole_nb[i+j*size_out] << " ";
		//}
		//ofile << endl;
	//}
	//for(unsigned int i(0); i<size_out; ++i){
		//ofile << dipole_nb[i] << " " ;
	//}
	ofile.close();
}

ostringstream& System::ref(ostringstream& name) const{
	name << this->param->get_ref() << "_link-prob"<< this->link_prob<< "_spring" << this->spring_coef << "_bend" << this->bend_coef 
	<< "_F" << this->F_magn << "_myo" << this->dip_prob<< "_ratio" << this->ratio << "_anchors-prob" << this->anchors_ratio;
	return name;
}

Couple System::get_dip(const unsigned int & i) const{
	return this->dip_list[i];
}

System::~System(){
	cerr << "System destructor entered" << endl;
}

//End of the class definitions

double compute_energy(const gsl_vector *var, void *param){
	
	System * sys=(System *)param;
	
	return sys->get_bend_E(var)+sys->get_spring_E(var)+sys->get_force_E(var);
	//return sys->get_bend_E(var)+sys->get_spring_E(var);
	//return sys->get_force_E(var);
}

void compute_gradient(const gsl_vector *var, void *param, gsl_vector * g){
	
	System *sys((System *)param);
	//cerr << "entered in gradient computing" << endl;
	for(unsigned int i(0); i<sys->get_latsize(); i++){
		if(sys->get_inpoint(i)){
			//cerr << i << endl;
			gsl_vector_set(g, 2*i, sys->get_bend_dE_x(var, i) +sys->get_spring_dE_x(var, i) + sys->get_force_dE_x(var, i));
			gsl_vector_set(g, 2*i+1, sys->get_bend_dE_y(var, i) +sys->get_spring_dE_y(var, i) + sys->get_force_dE_y(var, i));
			//gsl_vector_set(g, 2*i, sys->get_bend_dE_x(var, i) +sys->get_spring_dE_x(var, i));
			//gsl_vector_set(g, 2*i+1, sys->get_bend_dE_y(var, i) +sys->get_spring_dE_y(var, i));
			//cerr << i << " " << sys->get_points_size() << endl;
		}
		else{
			//cerr << "outpoint " << i << endl;
			gsl_vector_set(g, 2*i, 0.);
			gsl_vector_set(g, 2*i+1, 0.);
		}
	}
	//cerr << "out of first loop" << endl;
	//for(unsigned int i(sys->get_latsize()); i<sys->get_points_size(); i++){
		//if(sys->get_inpoint(i)){
			//gsl_vector_set(g, 2*i, sys->get_bend_dE_x(var, i) +sys->get_spring_dE_x(var, i));
			//gsl_vector_set(g, 2*i+1, sys->get_bend_dE_y(var, i) +sys->get_spring_dE_y(var, i));
			////cerr << i << " " << sys->get_points_size() << endl;
		//}
		//else{
			////cerr << "outpoint " << i << endl;
			//gsl_vector_set(g, 2*i, 0.);
			//gsl_vector_set(g, 2*i+1, 0.);
		//}
	//}
	//cerr << "out of first loop" << endl;
	for(unsigned int i(sys->get_latsize()); i<sys->get_points_size(); i++){
		gsl_vector_set(g, 2*i, sys->get_bend_dE_x(var, i) +sys->get_spring_dE_x(var, i));
		gsl_vector_set(g, 2*i+1, sys->get_bend_dE_y(var, i) +sys->get_spring_dE_y(var, i));
		//cerr << i << " " << sys->get_points_size() << endl;
	}
	//cerr << "done" << endl;
}

void compute_both(const gsl_vector *var, void * param, double * E, gsl_vector * g){
	
    *E=compute_energy(var, param);
    compute_gradient(var, param, g);
}
