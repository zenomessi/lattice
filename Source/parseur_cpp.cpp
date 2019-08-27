//
//  parseur.cpp
//  
//
//  Created by Zeno Messi on 20.03.18.
//

#include "parseur_cpp.hpp"
#include<random>
#include<stdlib.h>
#include<time.h>
#include<cmath>
#include <chrono>

Parseur::Parseur(const string& param_file_n): ref("0"), cell_path("0"), var_1d(0), var_matrix(0){
    
    const unsigned int dimension(2);
    cerr << this->liste.size() << endl;
    this->liste["Link"]=Link;
    this->liste["Rest"]=Rest;
    this->liste["Radius"]=Radius;
    this->liste["Cell"]=Zelle;
    this->liste["Spring"]=Spring;
    this->liste["Bend"]=Bend;
    this->liste["Myo"]=Myo;
    this->liste["AnchorRatio"]=AnchorRatio;
    this->liste["Force_M"]=FMagnitude;
    this->liste["Ratio"]=Ratio;
    this->liste["Matrix"]=Matrix;
    this->liste["Reference"]=Ref;
    this->liste["#"]=Comment;
    this->liste["seedLAT"]=SeedLat;
    this->liste["seedMYO"]=SeedMyo;
    
    var_1d=new double[liste.size()];
    for(unsigned int i(0); i<liste.size(); i++){
        var_1d[i]=0.;
    }
    
    this->seedLAT = std::chrono::system_clock::now().time_since_epoch().count();
    cerr << seedLAT << " seedLAT" << endl;
    ifstream in;
    in.open(param_file_n.c_str());
    
    string var_n;
    double var_val(0.);
    
    in>>var_n;
    
    //switch version
    while(!in.eof()){
		//cerr << var_n <<" " << this->liste.at(var_n) << endl;
		switch(this->liste.at(var_n)){
			
			case Comment :
				in.ignore(numeric_limits<streamsize>::max(), '\n');
				//var_1d[Comment]=1.;
				break;
			
			case Ref :
				in>>this->ref;
				var_1d[Ref]=1.;
				//cerr << this->get_ref() << endl;
				break;
				
			case Zelle :
				in>>this->cell_path;
				var_1d[Zelle]=1.;
				break;
				
			case Matrix :{
				in>> var_val;
				this->var_matrix=gsl_matrix_calloc(var_val,dimension);
                double temp_val(0.);
				for(unsigned int i(0); i<var_val ; i++){
					for(unsigned int j(0); j<dimension; j++){
                        in>>temp_val;
                        gsl_matrix_set(this->var_matrix, i,j,temp_val);
                    }
                }
                var_1d[Matrix]=1.;
                break;
            }
            default :
				in>> var_val;
				this->var_1d[this->liste.at(var_n)]=var_val;
		}
		//cerr << var_n << " " << var_val << endl;
		in>>var_n;
	}
	
    //if/else version
    //while(!in.eof()){
        //if(this->liste.at(var_n)==Comment){
            //in.ignore(numeric_limits<streamsize>::max(), '\n');
            //in>>var_n;
            ////cerr << "Ciao!" << endl;
            //continue;
        //}
        //else if(this->liste.at(var_n)==Ref){
            //in>>this->ref;
        //}
        //else{
            //in>>var_val;
            ////cerr << "Coucou" << endl;
            //cerr << var_n <<" "<< var_val << endl;
            //this->var_1d[this->liste.at(var_n)]=var_val;
            //if(this->liste.at(var_n)==Matrix){
                //this->var_matrix=gsl_matrix_calloc(var_val,dimension);
                //double temp(0.);
                    //for(unsigned int i(0); i<var_val ; i++){
                        //for(unsigned int j(0); j<dimension; j++){
                        //in>>temp;
                        //gsl_matrix_set(this->var_matrix, i,j,temp);
                    //}
                //}
            //}
        //}
        //in>>var_n;
    //}
    //in.close();
    cerr << this->get_ref() <<endl;
    //for(map<string,Variable>::iterator it(this->liste.begin()); it!=this->liste.end(); ++it){
			//cerr << it->first << " " << var_1d[it->second] << endl;
		//}
}

Parseur::~Parseur(){
    gsl_matrix_free(this->var_matrix);
    delete[] var_1d;
}

Parseur::Parseur(const Parseur& autre): ref(autre.get_ref()), cell_path(autre.get_cell()), var_1d(0), var_matrix(0), liste(autre.get_liste()), seedLAT(autre.seedLAT), seedMYO(autre.seedMYO){
	
	this->var_1d=new double[liste.size()];
	for(unsigned int i(0); i<liste.size(); i++){
		//this->var_1d[i]=autre.get(i);
		this->var_1d[i]= autre.var_1d[i];
	}
	this->var_matrix=gsl_matrix_alloc(autre.get_mat()->size1, autre.get_mat()->size2);
	gsl_matrix_memcpy(this->var_matrix, autre.get_mat());
	
}

//void Parseur::parse(ifstream& in, Lattice& sys){
//
//    string var_n;
//    double var_val(0.);
//    in>>var_n;
//    if(!in.eof() and var_n <Matrix){
//        in>>var_val;
//        //cerr << var_n <<" "<< var_val << endl;
//        //this->var_1d[this->liste.at(var_n)]=var_val;
//        sys.set(this->liste.at(var_n),var_val);
//    }
//    else{
//        in>>var_val;
//    }
//}

//get methods
double Parseur::get(const Variable& var) const{
    return this->var_1d[var];
}

double Parseur::get(const unsigned int& i) const{
	return this->var_1d[i];
}

string Parseur::get_ref() const{
    return this->ref;
}

gsl_matrix* Parseur::get_mat() const{
    return this->var_matrix;
}
string Parseur::get_cell() const{
	return this->cell_path;
}
unsigned Parseur::get_seed(const Variable& var) const{
	if(var==SeedLat) return this->seedLAT;
	if(var==SeedMyo) return this->seedMYO;
	return 0;
}

map<string, Variable> Parseur::get_liste() const{
	return this->liste;
}

//set methods
void Parseur::set_var1d(int variable, double valeur){
	this->var_1d[variable]=valeur;
}

void Parseur::set_seed(const Variable& var, const unsigned& seed){
	if(var==SeedLat) this->seedLAT=seed;
	if(var==SeedMyo) this->seedMYO=seed;
}

void Parseur::set_cell(const string& s){
	this->cell_path=s;
}
