//New cell_cpp class replaces cell2
//Created on 24.04.19

#ifndef cell_cpp_hpp
#define cell_cpp_hpp

#include <vector>
#include <fstream>
#include <string>
#include <utility>

class Cell{

	private:
		std::vector<std::pair<double, double> > points;
		double A;
		double cx;
		double cy;
		
	public:
		Cell(const std::string & path);
		
	//get methods
		unsigned int get_size() const;
		double get_area() const;
		double get_cx() const;
		double get_cy() const;
		std::pair<double, double> get_point(unsigned int) const;


};
#endif /* cell_cpp_hpp */
