//Header that declares the Cell class that contains the data from
// /home/zeno/Documents/data/alicia/
#ifndef cell2_hpp
#define cell2_hpp

#include <vector>
#include <fstream>
#include <string>
using std::vector;

enum Coord{X=0, Y=1};

class Cell{
	
	private:
		unsigned int size;
		double* pos;
		double area;
		double centrx;
		double centry;
	public:
	//construct/destruct
		Cell();
		Cell(int taille);
		Cell(unsigned int taille, const double tab[][2]);
		Cell(int taille, const double tab[]);
		Cell(vector< vector<double> > & tab);
		Cell(const Cell & c);//copy
		Cell(std::string path);
		~Cell();
		//access/manip
		void set(int location, int xy, double nb);
		unsigned int getsize() const;
		double get(int location, int xy) const;
		double get(int location, Coord xy) const;
		//operators
		double& operator()(int location, int xy);
		Cell& operator=(const Cell & c);
		
		//other methods
		bool incell(double posx, double posy) const;
		void setcharac();
		double getarea() const;
		//distance from an arbitrary point to the center of the cell
		double cdist(double x, double y) const;
		double cdist(unsigned int i) const;//point of the cell itself
		//position of the center of mass
		double getcx() const;
		double getcy() const;
		double centrangle(double posx, double posy) const;//gives the angular coordinate (theta of (R, theta)) of the point (posx, posy) in the centroid frame of reference
		double centrangle(int posx) const;//same for a node of the cell itself (notation like a normal C array centrangle(0) is the angle for the first node
		double getperimeter() const;//gives the perimeter of the cell
		double getdist(int i, int j) const;//gets the distance between the two nodes along the contour
		double getdist(double x, double y, unsigned int i) const;//gives the distance between the cell point i and the given (x, y) point
		unsigned int closest(double x, double y) const; //gives the index of the cell point closest to given (x, y) point
};

#endif /* cell2_hpp */
