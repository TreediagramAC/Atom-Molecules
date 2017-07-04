#include<vector>
using std::vector;
#include<sstream>
using std::ostringstream;
#include<algorithm>
using std::sort;
#include<numeric>
using std::reverse;
#include<string>
using std::string;
using std::to_string;
#include"atoms_and_molecules.h"
#include<iostream>
using std::cout;
using std::endl;

string Atom::get_id(){return(id_);}
double Atom::get_weight(){return(weight_);}
int Atom::get_bonds(){return(bonds_);}
 Molecule& operator + (const Atom& add1, const Atom& add2){
	Molecule* m=new(Molecule);
	m->add_atoms(add1);
	m->add_atoms(add2);
	return *m;
}

 Molecule& operator +(const Atom& a1, const Molecule& a2){
	 Molecule* m = new(Molecule);
	 *m = a2;
	 m->add_atoms(a1);
	 return *m;
 }
 Molecule& operator *(const Atom& a, int num){
	 Molecule* m=new(Molecule);
	 m->add_atoms(a, num);
	 return *m;
 }
 Molecule& operator *(int num, const Atom&a){
	 Molecule* m=new(Molecule);
	 m->add_atoms(a, num);
	 return *m;
 }
Hydrogenyx::Hydrogenyx() : Atom(){
   id_="Hh";
   weight_=2.4;
   bonds_=1;
  // cout << "This is an atom of " <<"Hh"<< endl;
}

Carbonyx::Carbonyx() : Atom(){
   id_="Cc";
   weight_=5.6;
   bonds_=4;
  // cout << "This is an atom of " <<"Cc"<< endl;
}

Sulphuryx::Sulphuryx() : Atom(){
   id_="Ss";
   weight_=10.8;
   bonds_=7;
//   cout << "This is an atom of " << "Ss"<< endl;
}

void Molecule::add_atoms(Atom add,int num){
   if(add.get_id()=="Hh"){
       atomH+=num;
   }
   if(add.get_id()=="Cc"){
       atomC+=num;
   }
   if(add.get_id()=="Ss"){
       atomS+=num;
   }
   //cout << "After adding atoms, the molecule has the formula " << this->get_formula() << " and weight " <<this->get_weight() << endl;
}

double Molecule::get_weight(){
   Hydrogenyx H;
   Carbonyx C;
   Sulphuryx S;
   double weight;
   weight=H.get_weight()*atomH+C.get_weight()*atomC+S.get_weight()*atomS;
   return weight;
}

string Molecule::get_formula(){
   string formula;
   Hydrogenyx H;
   Carbonyx C;
   Sulphuryx S;
   ostringstream num_str;
   if(atomC>0){
       formula+="Cc";
       if(atomC>1){
           num_str.str("");
           num_str<<atomC;
           formula+=num_str.str();
       }
   }
   if(atomH>0){
       formula+="Hh";
       if(atomH>1){
           num_str.str("");
           num_str<<atomH;
           formula+=num_str.str();
       }
   }
   if(atomS>0){
       formula+="Ss";
       if(atomS>1){
           num_str.str("");
           num_str<<atomS;
           formula+=num_str.str();
       }
   }
   return formula;
}



bool Molecule::operator == (const Molecule& comp){
   return(atomH==comp.atomH and atomC==comp.atomC and atomS==comp.atomS);
}


Molecule& Molecule::operator += (const Molecule& add){
   Hydrogenyx H;
   Carbonyx C;
   Sulphuryx S;
   add_atoms(H,add.atomH);
   add_atoms(C,add.atomC);
   add_atoms(S,add.atomS);
   return *this;
}

Molecule& Molecule::operator += (const Atom& add){
   add_atoms(add);
   return *this;
}

Molecule& Molecule::operator = (const Molecule& assign){
   atomH=assign.atomH;
   atomC=assign.atomC;
   atomS=assign.atomS;
   return *this;
}

Molecule &Molecule::operator+(const Molecule& assign)
{
	atomH += assign.atomH;
	atomC += assign.atomC;
	atomS += assign.atomS;
	return *this;
}
Molecule &Molecule::operator+(const Atom& a)
{
	this->add_atoms(a);
	return *this;
}
int Molecule::is_stable(){
  vector<long> bondlist;
  int fbond = 0;
  for (int i = 0;i<atomS;i++){
    bondlist.push_back(7);
  }
  for (int j = 0;j<atomC;j++){
    bondlist.push_back(4);
  }
  for (int k = 0;k<atomH;k++){
    bondlist.push_back(1);
  }

// This method is created by Dr. Joshua and he taught me in CSE232 helproom.
// It is used to deal with the "bond" vector to final form that we need to get remain bonds.
  bool job_finish = false;
  while (job_finish != true){
    sort(bondlist.begin(),bondlist.end());
    reverse(bondlist.begin(), bondlist.end());
    if((bondlist[0]==0) or (bondlist.empty())){
      return 0;
    }
    size_t bonds = bondlist[0];
    for( size_t i=1;i<bonds+1;i++){
      if(i>=bondlist.size()){
        job_finish = true;
        break;
      }
	  if (bondlist[i] == 0){
		  job_finish = true;
		  break;
	  }
      --bondlist[0];
      --bondlist[i];
    }
  }
  fbond = accumulate(bondlist.begin(), bondlist.end(),0);
  return fbond;
}

// Build function that is used for final evaluate_chemical_equation function
//to judge which methond the function should go.
vector<string> Molecule::splitstring(char flag, string s){
	vector<string> SetofString;
	unsigned int begin = 0;
	for (unsigned i = 0; i < s.size(); i++){
		if (((s.at(i) == flag) and (i > begin)) or (i == s.size() - 1)){
			string sub;
			if (s.size() == i + 1) sub = s.substr(begin, i - begin + 1);
			else  sub = s.substr(begin, i - begin);
			SetofString.push_back(sub);
			begin = i + 1;
		}
	}
	return SetofString;
}

//Build function to tranform molecule's formula in string from string type to a
//Molecule type. Thus it can help to do the operators in evaluate_chemical_equation
Molecule& Molecule::string_to_molecule(string str){
	Molecule* m = new(Molecule);
	unsigned int size = str.size();
	int first = 0;
	int last;
	for (unsigned int i = 2; i < size; i++){
		if (str[i] >= 'a' and str[i] <= 'z'){
			last = i - 2;
			string atom = str.substr(first, 2);
			int num = 1;
			string numstr;
			if (last - first == 1) num = 1;
			else{
				numstr = str.substr(first + 2, last - first - 1);
				num = std::stoi(numstr);
			}
			if (atom == "Hh"){
				Hydrogenyx h;
				m->add_atoms(h, num);
			}
			else if (atom == "Cc"){
				Carbonyx c;
				m->add_atoms(c, num);
			}
			else if (atom == "Ss"){
				Sulphuryx s;
				m->add_atoms(s, num);
			}
			first = last+1;
		}
		if (i == size - 1){
			last = size - 1;
			string atom = str.substr(first, 2);
			int num = 0;
			if (last - first == 1)
				num = 1;
			else{
				string numstr = str.substr(first + 2, last - first - 1);
				num = std::stoi(numstr);
			}
			if (atom == "Hh"){
				Hydrogenyx h;
				m->add_atoms(h, num);
			}
			else if (atom == "Cc"){
				Carbonyx c;
				m->add_atoms(c, num);
			}
			else if (atom == "Ss"){
				Sulphuryx s;
				m->add_atoms(s, num);
			}
			first = last+1;
		}
	}
	return *m;
}

//This is a functionthat have many condition inside itself, the first thing is to
//judge which kind it is so I use switch, can cases lead fucntion goes into different
//method that it should use to get the final output.
string Molecule::evaluate_chemical_equation(string str){
	string result = "";
	unsigned int size = str.size();
	int kind = 0;
	Molecule m = Molecule();
	Molecule m1 = Molecule();
	for (unsigned int i = 0; i < size; i++){
		if (str[i] == '-'){
			kind = 2;
			break;
		}
	}
	if (kind!= 2){
		for (unsigned int i = 0; i < size; i++){
			if (str[i] == '+'){
				kind = 1;
				break;
			}
		}
	}
	switch (kind){
	case 0:{
		result += "This is a Simple Molecule\n";
		m = string_to_molecule(str);
		result += "The molecule " + m.get_formula() + " has the weight ";
		result += to_string(m.get_weight()).substr(0,to_string(m.get_weight()).size()-5);
		if (m.is_stable() != 0)
			result += " and is unstable with " + to_string(m.is_stable()) + " loose bonds\n";
		else
			result += " and is stable!!!\n ";
		break;
	case 1:{
		vector<string> atomlists = splitstring(' ', str);
		unsigned int size = atomlists.size();
		for (unsigned int i = 0; i < size; i++){
			if (atomlists.at(i) != "+")
				m += string_to_molecule(atomlists.at(i));
		}
		result += "This is a Reaction\n";
		if (m.is_stable()!=0)
		result += "The resulting molecule is " + m.get_formula() + " and has the weight " + to_string(m.get_weight()).substr(0, to_string(m.get_weight()).size() - 5) +
			" and is unstable with " + to_string(m.is_stable()) + " loose bonds\n";
		else
			result += "The resulting molecule is " + m.get_formula() + " and has the weight " + to_string(m.get_weight()).substr(0, to_string(m.get_weight()).size() - 5) +
			" and is stable!!!\n";
		break;
	}
	case 2:{
		vector<string> atomlists = splitstring('-', str);
		string str1 = atomlists.at(0);
    string str2 ="";
		str2 = atomlists.at(1).substr(2,str2.size()-2);
		vector<string> list1 = splitstring(' ', str1);
		vector<string> list2 = splitstring(' ', str2);
		unsigned size1 = list1.size();
		unsigned size2 = list2.size();
		for (unsigned int i = 0; i < size1; i++){
			if (list1.at(i) != "+")
				m += string_to_molecule(list1.at(i));
		}
		for (unsigned int i = 0; i < size2; i++){
			if (list2.at(i) != "+")
				m1 += string_to_molecule(list2.at(i));
		}
		if (m == m1){
			if (m1.is_stable() == 0){
				result += "This is a Transformation\n";
				result += "The transformation is valid!!! and results in a stable molecule!!!\n";
				result += "The resulting molecule is " + m1.get_formula() + " and has the weight " + to_string(m1.get_weight()).substr(0, to_string(m.get_weight()).size() - 5) + "\n";
			}
			else{
				result += "This is a Transformation\n";
				result += "The transformation is valid!!! and results in an unstable molecule with "+
				to_string(m.is_stable()) + " loose bonds\n";
        result += "The resulting molecule is " + m1.get_formula() + " and has the weight " + to_string(m1.get_weight()).substr(0, to_string(m.get_weight()).size() - 5) + "\n";
			}
		}
		else{
			result+="This is a Transformation\n";
			result+=m.get_formula() +" cannot be transformed to " + m1.get_formula()+"\n";
		}
		break;
	}
	default:
		break;
	}
	}
	return result;
}
