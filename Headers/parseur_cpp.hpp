//
//  parseur.hpp
//  
//
//  Created by Zeno Messi on 20.03.18.
//

#ifndef parseur_hpp
#define parseur_hpp

#include <stdio.h>
#include"lattice_cell_cpp.hpp"
#include<iostream>
#include<fstream>
#include<map>
#include<string>
#include"variable.hpp"
#include<gsl/gsl_matrix.h>

using namespace std;

class Lattice;

class Parseur{
    
private:
    string ref;
    string cell_path;
    double* var_1d;
    gsl_matrix * var_matrix;
    map<string,Variable> liste;
    unsigned seedLAT;
    unsigned seedMYO;
public:
		Parseur(const string&);
		Parseur(const Parseur&);
		~Parseur();
		void parse(ifstream&, Lattice&);
		double get(const Variable&) const;
		double get(const unsigned int&) const;
		unsigned get_seed(const Variable&) const;
		string get_ref() const;
		string get_cell() const;
		gsl_matrix* get_mat() const;
		map<string, Variable> get_liste() const;
		void set_var1d(int variable, double valeur);
		void set_seed(const Variable&, const unsigned&);
		void set_cell(const string& s);
};



#endif /* parseur_hpp */
