//
//  lattice_cell_cpp.hpp
//  
//
//  Created by Zeno Messi on 05.12.18
//

#ifndef lattice_cell_cpp_hpp
#define lattice_cell_cpp_hpp

#include <iostream>
#include<fstream>
#include<sstream>
#include <gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include<random>
#include <vector>
#include<stdlib.h>
#include<time.h>
#include<cmath>
#include <chrono>
#include<iterator>
#include<deque>
#include <utility>
#include <tuple>
#include<string>
#include "parseur_cpp.hpp"
#include "variable.hpp"
#include "cell2.hpp"
#include<algorithm>
#include<limits>

using namespace std;

typedef std::pair<double,double> Point;
typedef std::tuple<unsigned int, unsigned int, unsigned int> Trio;
typedef std::pair<unsigned int, unsigned int> Couple;

class Parseur;

class Lattice{

	protected:
	
//         _   _       _ _         _          
//    __ _| |_| |_ _ _(_) |__ _  _| |_ ___ ___
//   / _` |  _|  _| '_| | '_ \ || |  _/ -_|_-<
//   \__,_|\__|\__|_| |_|_.__/\_,_|\__\___/__/
//                                
	//1D
	//the ones defined a priori in the param file
	double link_prob;
	double rest_length;
	double ratio;
	double anchors_ratio;
	Cell boundary;
	
	//the ones created randomly
	unsigned int insize;
	unsigned int outsize;
	unsigned int hinges_size;
	unsigned int edges_size;
	unsigned int latsize;
	unsigned int midsize;
	unsigned int lattice_hinges_size;
	double dist_cent_max;
	double dist_bound_max;
	
	//helpers
	Point cell_center;
	
	//vectors
	vector<bool> inpoints;
	vector<Point> points;
	vector<Trio> hinges;
	vector<pair<int, int> > indices;
	vector<vector<pair<unsigned int, unsigned int> > > hinges_first;//new to cpp version
	vector<vector<pair<unsigned int, unsigned int> > > hinges_mid;//new to cpp version
	vector<vector<pair<unsigned int, unsigned int> > > hinges_last;//new to cpp version
	vector<vector<unsigned int> > edges_first;//new to cpp version
	vector<vector<unsigned int> > edges_second;//new to cpp version
	
	public:
	
//              _   _            _    
//    _ __  ___| |_| |_  ___  __| |___
//   | '  \/ -_)  _| ' \/ _ \/ _` (_-<
//   |_|_|_\___|\__|_||_\___/\__,_/__/
//                                    
	//constructor
	Lattice(const Parseur&);
	
	//set the cell boundaries
	void set_bound(const string& s);
	
	bool get_inpoint(const unsigned int &) const;
	
	//get and set
	//plotting
	unsigned int get_hinge_point(const unsigned int &, const unsigned int &) const;
	double get_point(const unsigned int &, const unsigned int &) const;
	Point get_point(const unsigned int &) const;
	
	
	//sizes
	unsigned int get_hinges_size() const;
	unsigned int get_points_size() const;
	unsigned int get_latsize() const;
	unsigned int get_lat_hinges_size() const;
	
	//Destructor
	~Lattice();
};

//Depth-first algorithm implemented to determine what are the connected lattice sites
void DFS(const vector<unsigned int>& graph, const unsigned int& i, bool * discovered);
void DFS(const vector<Trio>& graph, const unsigned int& i, bool * discovered);
#endif /* lattice_cell_cpp_hpp */
