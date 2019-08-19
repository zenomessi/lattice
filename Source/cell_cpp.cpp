//New Cell class that replaces cell2

#include "../Headers/cell_cpp.hpp"
#include <cmath>
#include <iostream>

Cell::Cell(const std::string & path): A(0), cx(0), cy(0){
	
	std::ifstream ifile(path, std::ios::in);
	//check that the path is reached
	if(ifile.fail()){
		std::cerr << "Could open the file stream for creating the cell." << std::endl;
		std::cerr << "Check the path : " << path << std::endl;
	}
	else{
		//initialize the positions of the contour
		double tempx, tempy;
		ifile >> tempx >> tempy;
		while(!ifile.eof()){
			
			this->points.push_back(std::make_pair(tempx, tempy));
			ifile >> tempx >> tempy;
		}
		ifile.close();
		
		//create the characteristics
		
		//The area and the centroid in the same loop
		for(auto it(this->points.begin()); it!=this->points.end()-1; ++it){
			
			double temp(it->first*std::next(it,1)->second-std::next(it,1)->first*it->second);
			this->A+= temp;
			this->cx+=(it->first+std::next(it,1)->first)*temp;
			this->cy+=(it->second+std::next(it,1)->second)*temp;
			
			//this->A+=it->first*std::next(it,1)->second-std::next(it,1)->first*it->second;//for the area
			//this->cx+=(it->first+std::next(it,1)->first)*(it->first*std::next(it,1)->second-std::next(it,1)->first*it->second);
		}
		this->A+=points.back().first*points.front().second-points.front().first*points.back().second;
		this->A/=2.;
		this->cx+=(points.back().first+points.front().first)*(points.back().first*points.front().second-points.front().first*points.back().second);
		this->cy+=(points.back().second+points.front().second)*(points.back().first*points.front().second-points.front().first*points.back().second);
		this->cx/=6*this->A;
		this->cy/=6*this->A;
	}
}
//Get methods
unsigned int Cell::get_size() const{ return this->points.size();}

double Cell::get_area() const{ return this->A;}
double Cell::get_cx() const{ return this->cx;}
double Cell::get_cy() const{ return this->cy;}
std::pair<double, double> Cell::get_point(unsigned int i) const{
	if(i<this->points.size()) return this->points.at(i);
	else{
		std::cerr << "Cell limits exceeded, the cell size is " << this->points.size() << " . There is no index " << i << std::endl;
		throw 1;
	}
}
