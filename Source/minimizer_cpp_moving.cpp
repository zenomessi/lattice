//this program is made to compare the output of the cpp version to the old one
//it is made to run on scitas

#include "lattice_cell_cpp.hpp"
#include "system_myo_cpp.hpp"
#include <iostream>
#include <gsl/gsl_vector.h>
#include <vector>
#include <list>
#include <iomanip>
#include <climits>
#include <gsl/gsl_multimin.h>

int main(int argc, char* argv[]){
	
	//sempiternel
	ostringstream pathname("");
	ofstream ofile;
	//creating a Parseur instance
	string fname(argv[1]);
	Parseur param(fname);
	cerr << "out of parseur constructor" << endl;
	
	//    initialization of the minimizer
	size_t iter = 0;
	int status;
	const gsl_multimin_fdfminimizer_type *T;
	gsl_multimin_fdfminimizer *s;
		
		
	//declare the system
	System_General * sys(new System(param));
	cerr << "System created" << endl;
	
	//put the reference for the files in the pathname
	pathname << "./";
	sys->ref(pathname);
	
	//indexes of each hinge
	ofile.open(pathname.str()+"_hinges.dat", ofstream::out);
	if(ofile.fail()) cerr << "Problem" << endl;
	for(unsigned int i(0); i<sys->get_hinges_size(); i++){
		ofile << sys->get_hinge_point(i, 0) << "\t" << sys->get_hinge_point(i, 1) << "\t" << sys->get_hinge_point(i, 2) << endl;
	}
	ofile.close();
	
	//just the lattice
	ofile.open(pathname.str()+"_network.dat", ofstream::out);
	cerr << "printing network.dat" << endl;
	for(unsigned int i(0); i<sys->get_hinges_size(); i++){
		ofile << sys->get_point(sys->get_hinge_point(i,0),0) << "\t" << sys->get_point(sys->get_hinge_point(i,0), 1) << endl;
		ofile << sys->get_point(sys->get_hinge_point(i,1),0) << "\t" << sys->get_point(sys->get_hinge_point(i,1), 1) << endl;
		ofile << sys->get_point(sys->get_hinge_point(i,2),0) << "\t" << sys->get_point(sys->get_hinge_point(i,2), 1) << endl;
	}
	ofile.close();
	
	//positions of outpoints
	ofile.open(pathname.str()+"_outpoints.dat", ofstream::out);
	cerr << "printing outpoints.dat" << endl;
	for(unsigned int i(0); i<sys->get_points_size(); i++){
		if(!sys->get_inpoint(i)) ofile << sys->get_point(i,0) << "\t" << sys->get_point(i,1) << endl;
	}
	ofile.close();
	
	//indexes of outpoints
	ofile.open(pathname.str()+"_boundpoints.dat", ofstream::out);
	cerr << "printing boundpoints.dat" << endl;
	for(unsigned int i(0); i<sys->get_points_size(); i++){
		if(!sys->get_inpoint(i)) ofile << i  << endl;
	}
	ofile.close();
	
	//dipoles indexes
	ofile.open(pathname.str()+"_dipoles_indexes.dat", ofstream::out);
	cerr << "printing dipoles_indexes.dat" << endl;
	for(unsigned int i(0); i<sys->get_dip_list_size(); i++){
		ofile << sys->get_dip(i).first << "\t" << sys->get_dip(i).second << endl;
	}
	ofile.close();
	
	//dipoles positions, to plot
	ofile.open(pathname.str()+"_dipoles_debut.dat", ofstream::out);
	cerr << "printing dipoles_debut.dat" << endl;
	for(unsigned int i(0); i<sys->get_dip_list_size(); i++){
		ofile << sys->get_point(sys->get_dip(i).first, 0) << "\t" << sys->get_point(sys->get_dip(i).first, 1) << endl;
		ofile << sys->get_point(sys->get_dip(i).second, 0) << "\t" << sys->get_point(sys->get_dip(i).second, 1) << endl;
	}
	ofile.close();
	
	//create the gsl variable
	gsl_vector * var=gsl_vector_calloc(sys->get_points_size()*2);
	
	//shift the nodes a little bit
	unsigned seedminimizer(std::chrono::system_clock::now().time_since_epoch().count());
	cerr << "Seedminimizer " << seedminimizer << endl;
	default_random_engine generator(seedminimizer);

	uniform_real_distribution<double> distribution(-0.015,0.015);
	
	for(unsigned int i(0); i<sys->get_points_size(); i++){
		if(sys->get_inpoint(i)){
			gsl_vector_set(var, 2*i, distribution(generator));
			gsl_vector_set(var, 2*i+1, distribution(generator));
		}
	}
	ofile.open(pathname.str()+"_state_debut.dat", ofstream::out);
	cerr << "printing state_debut.dat" << endl;
	for(unsigned int i(0); i<sys->get_hinges_size(); i++){
		ofile << sys->get_state(var,sys->get_hinge_point(i,0),0) << "\t" << sys->get_state(var,sys->get_hinge_point(i,0), 1) << endl;
		ofile << sys->get_state(var,sys->get_hinge_point(i,1),0) << "\t" << sys->get_state(var,sys->get_hinge_point(i,1), 1) << endl;
		ofile << sys->get_state(var,sys->get_hinge_point(i,2),0) << "\t" << sys->get_state(var,sys->get_hinge_point(i,2), 1) << endl;
	}
	ofile.close();
	
	ofile.open(pathname.str()+"_positions_debut.dat", ofstream::out);
	cerr << "printing positions_debut.dat" << endl;
	for(unsigned int i(0); i<sys->get_points_size(); i++){
		ofile << sys->get_state(var, i, 0) << "\t" << sys->get_state(var, i, 1) << endl;
	}
	ofile.close();
	
	cerr << "End of shifts" << endl;
	//initialization of the functions and derivatives
	gsl_multimin_function_fdf fdf_functions;
	fdf_functions.n=sys->get_points_size()*2;
	fdf_functions.f=&compute_energy;
	fdf_functions.df=&compute_gradient;
	fdf_functions.fdf=&compute_both;
	fdf_functions.params=(void *)sys;
	cerr<< "End of functions definitions "<< endl;
	
	T = gsl_multimin_fdfminimizer_vector_bfgs2;
	s = gsl_multimin_fdfminimizer_alloc (T, sys->get_points_size()*2);
	gsl_multimin_fdfminimizer_set (s, &fdf_functions, var, 0.01, .9);//line minimization value
	cerr << "FDF minimizer created" << endl;
	iter=0;
	status=0;
	
	//store the norm of the gradient in a file
	ofile.open(pathname.str()+"_gradients.dat", ofstream::out);
	
	cerr << "Début de l'optimisation" << endl;
	vector<Point> store_pos;
	
	do{
		
		
		if(iter%1000==0){
			//plot the total gradient norm and the total energy
			ofile << gsl_blas_dnrm2(s->gradient) << endl;
		}
		iter++;
		status = gsl_multimin_fdfminimizer_iterate (s);
		
		if(status) break;
		status = gsl_multimin_test_gradient (s->gradient, 5e-3);
		
	}while(status == GSL_CONTINUE && iter < 1000000);
	cerr << "event: " << gsl_strerror (status) << endl;
	cerr << "iterations :" << iter << endl;
	cerr << "gradient " << gsl_blas_dnrm2(s->gradient) << endl;
	ofile << gsl_blas_dnrm2(s->gradient) << endl;
	ofile.close();
	//minfile.close();
	
	//plot the final state
	ofile.open(pathname.str()+"_state_final.dat", ofstream::out);
	cerr << "printing state_final.dat" << endl;
	//for(unsigned int i(0); i<sys->get_hinges_size(); i++){
	for(unsigned int i(0); i<sys->get_lat_hinges_size(); i++){
		ofile << sys->get_state(s->x,sys->get_hinge_point(i,0),0) << "\t" << sys->get_state(s->x,sys->get_hinge_point(i,0), 1) << endl;
		ofile << sys->get_state(s->x,sys->get_hinge_point(i,1),0) << "\t" << sys->get_state(s->x,sys->get_hinge_point(i,1), 1) << endl;
		ofile << sys->get_state(s->x,sys->get_hinge_point(i,2),0) << "\t" << sys->get_state(s->x,sys->get_hinge_point(i,2), 1) << endl;
	}
	ofile.close();
	
	ofile.open(pathname.str()+"_positions_final.dat", ofstream::out);
	cerr << "printing positions_final.dat" << endl;
	//This is put in the file for the CGAL version
	ofile << sys->get_points_size() << endl;
	for(unsigned int i(0); i<sys->get_points_size(); i++){
		ofile << sys->get_state(s->x, i, 0) << "\t" << sys->get_state(s->x, i, 1) << endl;
	}
	ofile.close();
	
	//forces, plottable version
	ofile.open(pathname.str()+"_forces.dat", ofstream::out);
	cerr << "printing forces.dat" << endl;
	for(unsigned int i(0); i<sys->get_points_size(); i++){
		if(!sys->get_inpoint(i)){
			ofile << sys->get_state(s->x, i, 0) << "\t" << sys->get_state(s->x, i, 1) << endl;
			ofile << sys->get_state(s->x, i, 0)-sys->get_spring_dE_x(s->x,i)-sys->get_bend_dE_x(s->x,i)-sys->get_force_dE_x(s->x,i)
			 << "\t" << sys->get_state(s->x, i, 1)-sys->get_spring_dE_y(s->x,i)-sys->get_bend_dE_y(s->x,i)-sys->get_force_dE_y(s->x,i) << endl;
		}
	}
	ofile.close();
	
	//dipoles positions, to plot
	ofile.open(pathname.str()+"_dipoles_final.dat", ofstream::out);
	cerr << "printing dipoles_final.dat" << endl;
	for(unsigned int i(0); i<sys->get_dip_list_size(); i++){
		ofile << sys->get_state(s->x, sys->get_dip(i).first, 0) << "\t" << sys->get_state(s->x, sys->get_dip(i).first, 1) << endl;
		ofile << sys->get_state(s->x, sys->get_dip(i).second, 0) << "\t" << sys->get_state(s->x, sys->get_dip(i).second, 1) << endl;
	}
	ofile.close();
	
	//forces distance vectorial
	ofile.open(pathname.str()+"_forces_norms.dat", ofstream::out);
	cerr << "printing forces_norms.dat" << endl;
	for(unsigned int i(0); i<sys->get_points_size(); i++){
		if(!sys->get_inpoint(i)){
			ofile << sys->get_state(s->x, i, 0)<< "\t" << sys->get_state(s->x, i, 1) <<"\t";
			ofile << sqrt(pow(sys->get_spring_dE_x(s->x,i)+sys->get_bend_dE_x(s->x,i)+sys->get_force_dE_x(s->x,i),2.)+pow(sys->get_spring_dE_y(s->x,i)+sys->get_bend_dE_y(s->x,i)+sys->get_force_dE_y(s->x,i),2.)) <<endl;
		}
	}
	ofile.close();
	
	////forces distances norms
	//ofile.open(pathname.str()+"_forces-distances_norms.dat", ofstream::out);
	//cerr << "printing forces-distances_norms.dat" << endl;
	//for(unsigned int i(0); i<sys->get_points_size(); i++){
		//if(!sys->get_inpoint(i)){
			//ofile << sqrt(pow(sys->get_state(s->x, i, 0),2.)+pow(sys->get_state(s->x, i, 1),2.)) <<endl;
		//}
	//}
	
	string dip_density_path(pathname.str()+"_dipoles_density.dat");
	sys->plot_dipoles_density(s->x, 5, dip_density_path);
	
	cerr << "done printing the system" << endl;
	
	//free the memory
	delete sys;
	gsl_multimin_fdfminimizer_free(s);
	gsl_vector_free(var);
		
	return 0;
}
