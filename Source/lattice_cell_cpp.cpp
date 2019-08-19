//
//  lattice_cell_cpp.cpp
//  
//
//  Created by Zeno Messi on 05.12.18
//
//  defines the class Lattice_cell, new version of the initialization
//
//#include "/home/messi/Documents/Lattice/Headers/lattice_cell_cpp.hpp"
#include "lattice_cell_cpp.hpp"

//constructors

Lattice::Lattice(const Parseur& param): rest_length(param.get(Rest)), link_prob(param.get(Link)), ratio(param.get(Ratio)), boundary(), dist_cent_max(0.), anchors_ratio(param.get(AnchorRatio)){
	
	cerr << "Constructing the Lattice..." << endl;
	
	//First create and initialize the PRNG
	default_random_engine generator(param.get_seed(SeedLat));
	//default_random_engine generator(2593568856);
	uniform_real_distribution<double> distribution(0.,1.);
	
	//The older version used temporary vectors, now we use directly the final vectors
	//Since vectors are objects, no need to initialize them in the initializer list
	
	this->set_bound(param.get_cell());
	cerr << this->boundary.getsize() << endl;
	//compute the maximal distance to the center of the cell
	for(unsigned int i(0); i<this->boundary.getsize(); i++){
		if(this->boundary.cdist(i)>dist_cent_max) dist_cent_max=this->boundary.cdist(i);
	}
	int max_dist_cent_i((int) ceil(dist_cent_max));//integer version
	cerr << "max dist_cent " << dist_cent_max << endl;
	//Put the coordinates of the vector matrix in const instances
	const double v1x(gsl_matrix_get(param.get_mat(),0,0));
	const double v1y(gsl_matrix_get(param.get_mat(),1,0));
	const double v2x(gsl_matrix_get(param.get_mat(),0,1));
	const double v2y(gsl_matrix_get(param.get_mat(),1,1));
	
	//create the points in the network, and the indices initial vectors
	//for(int i(-3*max_dist_cent_i); i<3*max_dist_cent_i; i++){
		//for(int j(-3*max_dist_cent_i); j<3*max_dist_cent_i; j++){
	for(int i(-10*max_dist_cent_i); i<10*max_dist_cent_i; i++){
		for(int j(-10*max_dist_cent_i); j<10*max_dist_cent_i; j++){
			double x(v1x*i+v2x*j);
			double y(v1y*i+v2y*j);
			if(this->boundary.incell(x,y)){
				this->points.push_back(make_pair(x,y));
				this->indices.push_back(make_pair(i,j));
			}
		}
	}
	this->latsize=this->points.size();//we only created lattice points
	
	//this boolean vector tells if a lattice point has neighbour (false) or not (true)
  vector<bool> seuls(latsize,true);//at first nobody has a neighbour
  
  //create now the edges
  for(unsigned int i_x(0); i_x<this->latsize; i_x++){
		for(unsigned int i_y(i_x+1); i_y<this->latsize; i_y++){
			//let's put the condition hard coded instead of the bond vector
			int ixiyi(this->indices[i_y].first-this->indices[i_x].first);
			int ixiyj(this->indices[i_y].second-this->indices[i_x].second);
			if((ixiyi==0 and ixiyj==1) or //(i,j) = (0,1)
			(ixiyi==1 and ixiyj==0) or //(i,j) = (0,1)
			(ixiyi==1 and ixiyj==-1)){//(i,j) = (1,-1)
				
				if(distribution(generator)<this->link_prob){
					//position of the midpoint
					points.push_back(make_pair(points[i_x].first+(v1x*ixiyi+v2x*ixiyj)/2.,points[i_x].second+(v1y*ixiyi+v2y*ixiyj)/2.));
					//indices of the three points of the new hinge
					hinges.push_back(make_tuple(i_x, points.size()-1, i_y));
					seuls[i_x]=false;
					seuls[i_y]=false;
				}
			}
		}
	}
	//Edges created
	cerr << "Edges created " << this->hinges.size() << " hinges before DFS" << endl;
	
	//magic part with the DFS algorithm (very similar to the non-cpp version)
	vector<bool*> disco;//dynamic array that will hold all the arrays of connected points
	bool * discotot(new bool [this->latsize]);//this array will hold all the points that were already tested, when the array is full, the whole lattice has been gone through
	
	//the array is initialized with the positions of all lone points, we don't need to test them
	for(unsigned int i(0); i<this->latsize; i++){
		discotot[i]=seuls[i];
  }
  
  //the loop that goes through the whole lattice creating arrays of connected lattice points
  bool fini(false);
  while(!fini){
		fini=true;
		unsigned int debut(0);
		//this loop is to find where to start the DFS algorithm
		for(unsigned int i(0); i<this->latsize; i++){
			if(discotot[i]==false){
				fini=false;
				for(unsigned int j(0); j<this->hinges.size(); j++){
					if(get<0>(this->hinges[j])==i or get<2>(this->hinges[j])==i){
					//if(get<0>(this->hinges[j])==i){
						debut=i;
						break;
					}
				}
				break;
			}
		}
		if(!fini){
			bool* discovered(new bool[this->latsize]);
			for(unsigned int i(0); i<this->latsize; i++){
				discovered[i]=false;
			}
			DFS(this->hinges, debut, discovered);
			disco.push_back(discovered);
			//add the newly discovered points to discotot
      for(unsigned int i(0); i<latsize; i++){
				discotot[i]=(discotot[i] or discovered[i]);
			}
		}
	}
	
	//now we choose the best
	unsigned int taille_meilleur(0);
	unsigned int meilleur(0);
	for(unsigned int i(0); i<disco.size(); i++){
			
		unsigned int taille_disco(0);
		for(unsigned int j(0); j<this->latsize; j++){
			if(disco[i][j]) taille_disco++;
		}
		cerr << taille_disco << endl;
		if(taille_disco>taille_meilleur){
			taille_meilleur = taille_disco;
			meilleur=i;
		}
	}
	cerr << meilleur << " " << disco.size() << endl;
	
	//keep only the best
	vector<Trio> newhinges(this->hinges);
	this->hinges.clear();
	for(auto it(newhinges.begin()); it!=newhinges.end(); ++it){
		if(disco[meilleur][get<0>(*it)]){
			this->hinges.push_back(*it);//if a node is in the best, than its neighbours as well!!
		}
	}
	newhinges.clear();
	
	delete[] discotot;
	seuls.clear();
	for(unsigned int i(0); i<disco.size(); i++){
		delete[] disco[i];
	}
	//DFS part finished!
	
	cerr <<"Out of DFS " << this->points.size() << " Points after DFS" << endl;
	
	//now remap!!
	
	//initialize the number of midpoints and edges
	//this->midsize= this->hinges.size();
	this->midsize= this->points.size()-this->latsize;
	this->edges_size= this->hinges.size()*2;
	cerr << this->midsize <<" " << this->latsize <<" " << this->points.size() <<endl;
	//real remapping
	bool* exists(new bool[latsize+midsize]);
	int * remap(new int[latsize+midsize]);
	for(unsigned int i(0); i<latsize+midsize; i++){
		exists[i]=false;
	}
	for(auto it(this->hinges.begin()); it!=this->hinges.end(); ++it){
		exists[get<0>(*it)]=true;
		exists[get<1>(*it)]=true;
		exists[get<2>(*it)]=true;
	}
	unsigned int new_pos(0);
	for(unsigned int i(0); i<this->latsize+this->midsize; i++ ){
		if(exists[i]){
			remap[i]=new_pos;
			new_pos++;
		}
		else remap[i]=-1;
	}
	vector<Point> newpoints(new_pos);
	vector<pair<int, int> > newindices(new_pos);
	
	for(unsigned int i(0); i<this->latsize; i++){
		if(remap[i]>-1){
			newpoints[(unsigned int) remap[i]]=this->points[i];
			newindices[(unsigned int) remap[i]]=this->indices[i];
		}
	}
	for(unsigned int i(this->latsize); i<this->latsize+this->midsize; i++){
		if(remap[i]>-1){
			newpoints[(unsigned int) remap[i]]=this->points[i];
		}
	}
	this->points=newpoints;
	this->indices=newindices;
	newpoints.clear();
	newindices.clear();
	
	for(unsigned int i(0); i<this->hinges.size(); i++){
		int deplace(remap[get<0>(this->hinges[i])]);
		get<0>(this->hinges[i])=deplace;
		deplace=remap[get<1>(this->hinges[i])];
		get<1>(this->hinges[i])=deplace;
		deplace=remap[get<2>(this->hinges[i])];
		get<2>(this->hinges[i])=deplace;
	}
	delete[] exists;
	delete[] remap;
	//done remapping!
	
	cerr << "coucou" << " out of remapping" << endl;
	
	this->lattice_hinges_size=this->hinges.size();//number of hinges of format lattice point-midpoint-lattice point
	//update midsize and latsize
	this->midsize= this->lattice_hinges_size;
	this->latsize=this->points.size()-midsize;
	cerr <<" again " << this->midsize <<" " << this->latsize <<" " << this->points.size() <<endl;
	
//              _               _       _       
//   ___  _   _| |_ _ __   ___ (_)_ __ | |_ ___ 
//  / _ \| | | | __| '_ \ / _ \| | '_ \| __/ __|
// | (_) | |_| | |_| |_) | (_) | | | | | |_\__ \
//  \___/ \__,_|\__| .__/ \___/|_|_| |_|\__|___/
//                 |_|     

	this->insize=midsize;
	this->outsize=0;
	
	this->inpoints= vector<bool>(latsize+midsize,true);//all points are in the cell by construction
	
	//original version with a random distribution of outpoints
	for(unsigned int i(0); i<this->latsize; i++){
		if(distribution(generator)>this->anchors_ratio){
			this->inpoints[i]=false;
			this->outsize++;
		}
		else this->insize++;
	}
	
	
	//Now we need to add the other hinges, the ones that have a lattice point in the middle (form midpoint-lattice point-midpoint)
	//first loop is on the first lattice point (starting point) of the hinge
	for(unsigned int hinges_index(0); hinges_index<this->lattice_hinges_size; hinges_index++){
		//the vector indicating the lattice direction of the hinge
		int u_1_i(indices[get<2>(this->hinges[hinges_index])].first-indices[get<0>(this->hinges[hinges_index])].first);
		int u_1_j(indices[get<2>(this->hinges[hinges_index])].second-indices[get<0>(this->hinges[hinges_index])].second);
		
		//check if the end of the first hinge is the beginning of the second
		for(unsigned int test_index(hinges_index+1); test_index<this->lattice_hinges_size; test_index++){
			if(get<0>(this->hinges[test_index])==get<2>(this->hinges[hinges_index])){
				//same thing as u_1_a but for the second "test" hinge
				int u_2_i(indices[get<2>(this->hinges[test_index])].first-indices[get<0>(this->hinges[test_index])].first);
				int u_2_j(indices[get<2>(this->hinges[test_index])].second-indices[get<0>(this->hinges[test_index])].second);
				
				if(u_1_i==u_2_i && u_1_j==u_2_j){
					this->hinges.push_back(make_tuple(get<1>(this->hinges[hinges_index]), get<2>(this->hinges[hinges_index]), get<1>(this->hinges[test_index])));
				}
			}
		}
	}
	
	this->hinges_size = this->hinges.size();
	cerr << "Lattice hinges size " << this->lattice_hinges_size << ", hinges size " << this->hinges_size << endl;
	cerr << "coucou" << " out of putting other hinges" << endl;
	
	//we change the indices to put the midpoints as well
	for(unsigned int i(0); i<new_pos; i++){
		indices[i].first*=2;
		indices[i].second*=2;
	}
	for(unsigned int i(0); i< this->lattice_hinges_size; i++){
		//cerr << get<1>(this->hinges[i]) << endl;
		indices[get<1>(this->hinges[i])].first=indices[get<0>(this->hinges[i])].first+(indices[get<2>(this->hinges[i])].first-indices[get<0>(this->hinges[i])].first)/2;
		indices[get<1>(this->hinges[i])].second=indices[get<0>(this->hinges[i])].second+(indices[get<2>(this->hinges[i])].second-indices[get<0>(this->hinges[i])].second)/2;
	}
	cerr << "coucou" << " out of all" << endl;
	
	
  cerr << "The Lattice has been created" << endl;
  cerr << "Now we create the helping vectors" << endl;
  cerr << this->hinges_size << " " << this->hinges.size() << endl;
  this->hinges_first= vector<vector<pair<unsigned int, unsigned int> > >(this->points.size());
  this->hinges_mid=vector<vector<pair<unsigned int, unsigned int> > >(this->points.size());
  this->hinges_last=vector<vector<pair<unsigned int, unsigned int> > >(this->points.size());
  for(unsigned int i(0); i< this->hinges_size; i++){
		this->hinges_first[get<0>(this->hinges[i])].push_back(make_pair(get<1>(this->hinges[i]), get<2>(this->hinges[i])));
		this->hinges_mid[get<1>(this->hinges[i])].push_back(make_pair(get<0>(this->hinges[i]), get<2>(this->hinges[i])));
		this->hinges_last[get<2>(this->hinges[i])].push_back(make_pair(get<0>(this->hinges[i]), get<1>(this->hinges[i])));
	}
	
	this->edges_first =vector<vector<unsigned int> >(this->points.size());
	this->edges_second =vector<vector<unsigned int> >(this->points.size());
	//cerr << "hello" << endl;
	for(unsigned int i(0); i<this->lattice_hinges_size; i++){
		this->edges_first[get<0>(this->hinges[i])].push_back(get<1>(this->hinges[i]));
		this->edges_first[get<1>(this->hinges[i])].push_back(get<2>(this->hinges[i]));
		this->edges_second[get<1>(this->hinges[i])].push_back(get<0>(this->hinges[i]));
		this->edges_second[get<2>(this->hinges[i])].push_back(get<1>(this->hinges[i]));
	}
	
	cerr << "Done with the helping vectors, hope they will help..." << endl;
}


//new version with the ratio as an object member
void Lattice::set_bound(const string& s){
	Cell uncentered(s);
	cerr << "Cell file : "<< s <<" ratio :" << this->ratio << endl;
	double center[2]={uncentered.getcx(),uncentered.getcy()};
	this->cell_center.first=center[0];
	this->cell_center.second=center[1];
	//for(int i(0); i<uncentered.getsize(); i++){
		//for(int j(0); j<2; j++){
			//uncentered.set(i, j, (uncentered.get(i, j)-center[j])/this->ratio);
		//}
	//}
	uncentered.setcharac();
	this->boundary=uncentered;
	//cerr << "hello" << endl;
}

bool Lattice::get_inpoint(const unsigned int & i) const{
	return this->inpoints[i];
}

unsigned int Lattice::get_hinge_point(const unsigned int & i, const unsigned int & j) const{
	switch(j){
		case 0: return get<0>(this->hinges[i]);
		case 1: return get<1>(this->hinges[i]);
		case 2: return get<2>(this->hinges[i]);
	}
	cerr << "ERROR: in get_hinge_point out of the bonds, second argument must be [0,2], 0 returned"<< endl;
	return 0;
}
double Lattice::get_point(const unsigned int & i, const unsigned int & j) const{
	if(j==0) return this->points[i].first;
	return this->points[i].second;
}
Point Lattice::get_point(const unsigned int & i) const{
	return this->points[i];
}

unsigned int Lattice::get_hinges_size() const{
	return this->hinges_size;
}

unsigned int Lattice::get_points_size() const{
	return this->latsize+this->midsize;
}

unsigned int Lattice::get_latsize() const{
	return this->latsize;
}

unsigned int Lattice::get_lat_hinges_size() const{
	return this->lattice_hinges_size;
}
//out of the class

void DFS(const vector<unsigned int>& graph, const unsigned int& i, bool * discovered){
    
    discovered[graph[i]]=true;
    for(unsigned int j(0); j<graph.size()/3; j++){
        
        if(graph[3*j]==graph[i]){
            if(!discovered[graph[3*j+2]]) DFS(graph,3*j+2,discovered);
        }
        if(graph[3*j+2]==graph[i]){
            if(!discovered[graph[3*j]]) DFS(graph,3*j,discovered);
        }
    }
}

void DFS(const vector<Trio>& graph, const unsigned int& i, bool * discovered){
	
	discovered[i]=true;
	for(unsigned int j(0); j< graph.size(); j++){
		
		if(get<0>(graph[j])==i){
			if(!discovered[get<2>(graph[j])]) DFS(graph, get<2>(graph[j]), discovered);
		} 
		if(get<2>(graph[j])==i){
			if(!discovered[get<0>(graph[j])]) DFS(graph, get<0>(graph[j]), discovered);
		}
	}
}

Lattice::~Lattice(){
	cerr << "The Lattice instance has been correctly destroyed" << endl;
}
