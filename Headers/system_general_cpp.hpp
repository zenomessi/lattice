//
//  system_general_cpp.hpp
//  
//
//  Created by Zeno Messi on 28.01.19
//
#ifndef system_general_cpp_hpp
#define system_general_cpp_hpp

#include <iostream>
#include "variable.hpp"
#include "parseur_cpp.hpp"
#include "lattice_cell_cpp.hpp"

using namespace std;

class System_General: public Lattice{
	
	protected:
		double spring_coef;
		double bend_coef;
		
		Parseur* param;
	
	public:
	
	//Methods
	//Constructors
		System_General(const Parseur &);
	
	//get/set
		double get_state(const gsl_vector * var, const unsigned int& i, unsigned int j) const;
		double get_var(const gsl_vector * var, const unsigned int& i, unsigned int j) const;

	//compute the energy
		double get_spring_E(const gsl_vector * var) const;
		double get_spring_dE_x(const gsl_vector* var, const unsigned int & i) const;
		double get_spring_dE_y(const gsl_vector* var, const unsigned int & i) const;
		double get_bend_E(const gsl_vector* var) const;
		double get_bend_dE_x(const gsl_vector* var, const unsigned int& i) const;
		double get_bend_dE_y(const gsl_vector* var, const unsigned int& i) const;
		double compute_bend_dE_x(const double& r_x, const double& r_y, const double& v_1_x, const double& v_1_y, const double& v_2_x, const double v_2_y) const;
		double compute_bend_dE_y(const double& r_x, const double& r_y, const double& v_1_x, const double& v_1_y, const double& v_2_x, const double v_2_y) const;
		virtual double get_force_E(const gsl_vector* var) const=0;//to be redefined in sub class
		virtual double get_force_dE_x(const gsl_vector* var, const unsigned int& i) const=0;//to be redefined in sub class
		virtual double get_force_dE_y(const gsl_vector* var, const unsigned int& i) const=0;//to be redefined in sub class
		
		virtual ostringstream& ref(ostringstream& name) const=0;//to be redefined in sub class
		virtual unsigned int get_dip_list_size() const=0;//to be redefined in sub class
		virtual Couple get_dip(const unsigned int&) const=0;//to be redefined in sub class
		virtual void plot_dipoles_density(const gsl_vector* var, const unsigned int& side_out, const string& s) const=0;//to be redefined in sub class
		virtual void set_dipoles(const string& filepath)=0;//to be redefined in sub class
	//Destructor
		virtual ~System_General();
};

#endif //system_general_cpp_hpp
