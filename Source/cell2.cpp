//Class cell
//#include "/home/messi/Documents/Lattice/Headers/cell2.hpp"
#include "cell2.hpp"//scitas
#include <cmath>
#include <iostream>

using namespace std;

//construct/copy/destruct
Cell::Cell(): size(0), pos(new double[1]), area(0.), centrx(0.), centry(0.){
	
	std::cerr << "The cell created is empty" << std::endl;
	}
Cell::Cell(int taille): size(taille), pos(new double [2*taille]), area(0.), centrx(0.), centry(0.){
	
	std::cerr << "The cell created has no area or centre of mass" << std::endl
	<< "Use the setcharac function to create them once it is instanciated" << std::endl;
	}
Cell::Cell(unsigned int taille,const double tab[][2]): size(taille), pos(new double [2*taille]){
	
	//putting the positions in the array
	for(unsigned int i(0); i<taille; i++){
		
		pos[i]=tab[i][0];
		
	}
	for(unsigned int i(0); i<taille; i++){
		
		pos[i+taille]=tab[i][1];
	}
	//creating the other characteristics of the cell
	//old fashion
	//double temp(0);
	
	//for(unsigned int i(0); i<taille; i++){
		//temp+=pos[i]*pos[(i+1)%taille+taille]-pos[(i+1)%taille]*pos[i+taille];
	//}
	//area= temp/2.;//create the area
	
	//double temp2(0);
	//temp=0;
	
	//for(unsigned int i(0); i<taille; i++){
		//temp+=(pos[i]+pos[(i+1)%taille])*(pos[i]*pos[(i+1)%taille+taille]-pos[(i+1)%taille]*pos[i+taille]);//increment for x coordinate
		//temp2+=(pos[i+taille]+pos[(i+1)%taille+taille])*(pos[i]*pos[(i+1)%taille+taille]-pos[(i+1)%taille]*pos[i+taille]);//increment for y coordinate
	//}
	////create the cell center
	//centrx=temp/6./area;
	//centry=temp2/6./area;
	
	//Creating the characteristics with the setcharac function
	this->setcharac();
	
}
Cell::Cell(int taille, const double tab[]): size(taille), pos(new double [2*taille]){
	
	for(unsigned int i(0); i<2*taille; i++){
		pos[i]=tab[i];
	}
	
	this->setcharac();
}
Cell::Cell(vector< vector<double> > & tab): size(tab.size()), pos(new double [2*tab.size()]){
	for(unsigned int i(0); i<tab.size(); i++){
		pos[i]=tab[i][0];
		pos[i+tab.size()]=tab[i][1];
	}

	this->setcharac();
}
//copy constructor
Cell::Cell(const Cell & c): size(c.size), area(c.area), centrx(c.centrx), centry(c.centry), pos(new double[2*c.size]){
	
	
	for(unsigned int i(0); i<size*2; i++){
		
		pos[i]=(c.pos)[i];
	}
	
}

Cell::Cell(std::string path): size(0), pos(new double[1]), area(0.), centrx(0.), centry(0.){
	
	ifstream ifile(path.c_str(), std::ifstream::in);
	vector<vector< double> > vec;
	
	double temp(0.);
	unsigned int i(0);
	ifile>> temp;
	while(!ifile.eof()){
		//cerr << i << endl;
		vec.push_back(vector<double>(2));
		vec[i][0]=temp;
		ifile >> vec[i][1];
		ifile >> temp;
		i++;
	}
	ifile.close();
	*this = Cell(vec);
	
	vec.clear();
	
}

Cell::~Cell(){
	delete[] pos;
	pos =0;
}
//access/manip
void Cell::set(int location, int xy, double nb){
	
	pos[location+xy*this->size]=nb;
}
unsigned int Cell::getsize() const{
	return this->size;
}
double Cell::get(int location, int xy) const{
	
	return pos[location+xy*this->size];
}
double Cell::get(int location, Coord xy) const{
	
	return pos[location+xy*this->size];
}
//operators
double& Cell::operator()(int location, int xy){
	
	return pos[location+xy*this->size];
}
Cell& Cell::operator=(const Cell & c){
	if(&c != this){
		this->size=c.size;
		if(this->pos!=0) delete[] this->pos;
		pos = new double [2*size];
		for(unsigned int i(0); i< 2*size; i++){
			
			pos[i]=(c.pos)[i];
		}
		this->area=c.area;
		this->centrx=c.centrx;
		this->centry=c.centry;
	}
	return *this;
}
//other methods

//says if point (posx, posy) is in the cell.
bool Cell::incell(double posx, double posy) const {
	
	int crossnb(0);
	
	//loop over the vertices
	for(unsigned int i(0); i<size; i++){
		//if(posx<pos[i%size]){//checks if the vertex is to the right of the point
			////checks if the demi horizontal line crosses the last edge
			//if(pos[(i-1)%size+size]<posy && pos[(i-1)%size>posx]){//upward edge
				//if(pos[i%size+size]>posy) crossnb++;
			//}
			////downward edge
			//else if(pos[i%size+size]<posy) crossnb++;
		//}
		//if(posx<pos[i%size] && posy<pos[i%size+size] && posx<pos[(i-1)%size] && posy>pos[(i-1)%size+size]) crossnb++;
		//else if(posx<pos[i%size] && posy>pos[i%size+size] && posx<pos[(i-1)%size] && posy<pos[(i-1)%size+size]) crossnb++;
		//else if(posx<pos[i%size] && posy>pos[i%size+size] && posx>pos[(i-1)%size] && posy<pos[(i-1)%size+size]) crossnb++;
		//else if(posx>pos[i%size] && posy<pos[i%size+size] && posx<pos[(i-1)%size] && posy>pos[(i-1)%size+size]) crossnb++;
		
		//declare the points as variables
		double x1(pos[i]);
		double x2(pos[(i+1)%size]);
		double y1(pos[i+size]);
		double y2(pos[(i+1)%size+size]);
		if((y1<=posy && y2>posy) //upward edge
		or (y1>posy && y2<=posy)){//downward edge
			double vt((posy-y1)/(y2-y1));
			if(posx< x1+vt*(x2-x1)) crossnb++;
		}
		//if(x1!=x2 && y1!=y2){
			//double vx(x2-x1);
			//double vy(y2-y1);
			//double pente(vy/vx);
			//if (pente ==0.) cerr << "divide by 0" << endl;
			//double ordonn(y1-x1*pente);
			//if(posy< y2 && posy>=y1){//upward edge
				//if((posy-ordonn)/pente>posx) crossnb++;
			//}
			//else if(posy>=y2 && posy< y1){//downward edge
				//if((posy-ordonn)/pente>posx) crossnb++;
			//}
		//}
		//else if(posx<x1 && ((posy< y2 && posy>=y1) or (posy>=y2 && posy< y1))) crossnb++;
		
		//if(posx<=x2){
			
			//if(posy>=y2){
				
				//if(posy<y1){
					//if(posx<x1) crossnb++;//first case, normal downward edge
					//else {
						//double vx(x2-x1);
						//double vy(y2-y1);
						//double pente(vy/vx);
						//double ordon(y1-x1*pente);
						//if(posx*pente+ordon>posy) crossnb++;//second case, i-1 > <, i < >
					//}
				//}
				
			//}
			//else if(posx<x1 && posy>y1) crossnb++;//third case, normal upward edge
		//}
		//else if(posy<=y2 && posx<x1 && posy>y1){
			//double vx(x2-x1);
			//double vy(y2-y1);
			//double pente(vy/vx);
			//double ordon(y1-x1*pente);
			//if(posx*pente+ordon>posy) crossnb++;//fourth case, i-1 < >, i > <
		//}
	}
	return (crossnb&1);//& operator copies a bit to the result of it exists in both operands
}
//This method creates the characteristic values of the cell (area, center of mass)
void Cell::setcharac(){
	
	double temp(0.), temp2(0.);
	
	for(unsigned int i(0); i<this->size; i++){
		temp+=this->pos[i]*this->pos[(i+1)%this->size+this->size]-this->pos[(i+1)%this->size]*this->pos[i+this->size];
	}
	this->area= temp/2.;//create the area
	
	temp=0.;//reset temp to 0;
	
	//stupidly costly version
	
	for(unsigned int i(0); i<this->size; i++){
		temp+=(this->pos[i]+this->pos[(i+1)%this->size])*(this->pos[i]*this->pos[(i+1)%this->size+this->size]-this->pos[(i+1)%this->size]*this->pos[i+this->size]);//increment for x coordinate
		temp2+=(this->pos[i+this->size]+this->pos[(i+1)%this->size+this->size])*(this->pos[i]*this->pos[(i+1)%this->size+this->size]-this->pos[(i+1)%this->size]*this->pos[i+this->size]);//increment for y coordinate
	}
	//create the cell center
	this->centrx=temp/6./this->area;
	this->centry=temp2/6./this->area;
	
	////new version for computing the centroid
	//for(unsigned int i(0); i<this->size; i++){
		//temp+=pos[i];
		//temp2+=pos[i+this->size];
	//}
	//this->centrx=temp/this->size;
	//this->centry=temp2/this->size;
	
}

//gives the area of the cell
double Cell::getarea() const {
	
	//old fashion without an area attribute
	//double res(0);
	//for(int i(0); i<size; i++){
		//res+=pos[i]*pos[(i+1)%size+size]-pos[(i+1)%size]*pos[i+size];
	//}
	//return res/2.;
	
	//new fashion with an attribute
	return this->area;
}

//distance from an arbitrary point to the center of the cell

double Cell::cdist(double x, double y) const{
	
	return sqrt((x-this->centrx)*(x-this->centrx)+(y-this->centry)*(y-this->centry));
}
//special case for a node of the cell itself uses the general method.
double Cell::cdist(unsigned int i) const{
	
	return cdist(this->pos[i],this->pos[i+this->size]);
	
}

double Cell::getcx() const{
	
	//old fashion with a method that computes every time the centre of the cell
	//double res(0);
	//for(int i(0); i<size; i++){
		//res+=(pos[i]+pos[(i+1)%size])*(pos[i]*pos[(i+1)%size+size]-pos[(i+1)%size]*pos[i+size]);
	//}
	//return res/6./this->area;
	
	//New fashion with an attribute
	return this->centrx;
}

double Cell::getcy() const{
	
	//old fashion with a method that computes every time the centre of the cell
	//double res(0);
	//for(int i(0); i<size; i++){
		//res+=(pos[i+size]+pos[(i+1)%size+size])*(pos[i]*pos[(i+1)%size+size]-pos[(i+1)%size]*pos[i+size]);
	//}
	//return res/6./this->area;
	
	//New fashion
	return this->centry;
}
//atan gives a result [-PI/2;PI/2] so we map it on [-PI/2;3PI/2] 
double Cell::centrangle(double posx, double posy) const{
	
	if(posx<this->centrx){
		return atan((posy-this->centry)/(posx-this->centrx))+M_PI;
	}
	else return atan((posy-this->centry)/(posx-this->centrx));
}

double Cell::centrangle(int posx) const{
	
	if(pos[posx]<this->centrx){
		return atan((pos[posx+size]-this->centry)/(pos[posx]-this->centrx))+M_PI;
	}
	else return atan((pos[posx+size]-this->centry)/(pos[posx]-this->centrx));
}


double Cell::getperimeter() const{
	
	double res(0.);
	
	res+=sqrt((pos[0]-pos[size-1])*(pos[0]-pos[size-1])+(pos[size]-pos[2*size-1])*(pos[size]-pos[2*size-1]));
	
	for(unsigned int i(1); i<size; i++){
		
		res+=sqrt((pos[i]-pos[i-1])*(pos[i]-pos[i-1])+(pos[size+i]-pos[size+i-1])*(pos[size+i]-pos[size+i-1]));
	}
	return res;
}

double Cell::getdist(int i, int j) const{//attention i < j on contour
	
	double res(0.);
	if((j-i+size)%size*2>size){
		int temp(i);
		i=j;
		j=temp;
	}
	for(unsigned int k(1); k<=(j-i+size)%size; k++){
		res+=sqrt((pos[(i+k)%size]-pos[(i+k-1)%size])*(pos[(i+k)%size]-pos[(i+k-1)%size])+(pos[size+(i+k)%size]-pos[size+(i+k-1)%size])*(pos[size+(i+k)%size]-pos[size+(i+k-1)%size]));
	}
	
	return res;
}

double Cell::getdist(double x, double y, unsigned int i) const{
	return sqrt((this->get(i, 0)-x)*(this->get(i, 0)-x)+(this->get(i, 1)-y)*(this->get(i, 1)-y));
}

unsigned int Cell::closest(double x, double y) const{
	
	unsigned int res(0);
	double dist(this->getdist(x,y,0));
	for(unsigned int i(1); i< this->getsize(); i++){
		if(this->getdist(x,y,i)<dist){
			dist= this->getdist(x,y,i);
			res=i;
		}
	}
	//cerr << dist << endl;
	return res;
}
