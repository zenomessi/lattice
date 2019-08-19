//
//  system_myo_cpp.hpp
//  
//	
//  Created by Zeno Messi on 29.01.19
//

#ifndef system_myo_cpp_hpp
#define system_myo_cpp_hpp

#include <iostream>
#include "lattice_cell_cpp.hpp"
#include "system_general_cpp.hpp"
#include "variable.hpp"
#include "parseur_cpp.hpp"

class System: public System_General{
	
	protected:
		double F_magn;
		vector<Couple> dip_list;
		vector<vector<unsigned int> > dip_list_other;//new to cpp version
		double dip_prob;
		double min_dist;
		
	public:
		//constructor/destructor
		System(Parseur& par);
		System(Parseur& par, const unsigned int & i);//for the Toy model
		~System();
		
		//sizes
		unsigned int get_dip_list_size() const;
		
		//compute the energy
		double get_force_E(const gsl_vector* var) const;
		double get_force_dE_x(const gsl_vector* var, const unsigned int& index) const;
		double get_force_dE_y(const gsl_vector* var, const unsigned int& index) const;
		
		//dipoles density
		void plot_dipoles_density(const gsl_vector* var, const unsigned int& side_out, const string& s) const;
		
		ostringstream& ref(ostringstream& name) const;
		
		//plot the system
		Couple get_dip(const unsigned int &) const;
		
		//set methods
		void set_dipoles(const string& filepath);//mehtod used for the moving cells to initialize the dipoles with previous distribution
};

double compute_energy(const gsl_vector *var, void *sys);
void compute_gradient(const gsl_vector *var, void *sys, gsl_vector * g);
void compute_both(const gsl_vector *var, void * param, double * E, gsl_vector * g);

#endif /*system_myo_cpp_hpp */
