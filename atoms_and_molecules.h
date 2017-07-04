#ifndef ATOM_AND_MOLECULES_H
#define ATOM_AND_MOLECULES_H
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

class Molecule;

class Atom{
protected:
  string id_;
  double weight_;
  int bonds_;
 public:
  string get_id();
  int get_bonds();
  void use_bonds();
  double get_weight();
  friend Molecule& operator +(const Atom&,const Atom&);
  friend Molecule& operator +(const Atom&, const Molecule &);
  friend Molecule& operator *(const Atom&, int num);
  friend Molecule& operator*(int num, const Atom&);
};

class Hydrogenyx : public Atom{
public:
 Hydrogenyx();
};

class Carbonyx : public Atom{
public:
 Carbonyx();
};

class Sulphuryx : public Atom{
public:
 Sulphuryx();
};


class Molecule{
private:
 int atomH = 0;
 int atomC = 0;
 int atomS = 0;
public:
 Molecule()=default;
 Molecule(const Molecule& cpy){atomH=cpy.atomH;atomC=cpy.atomC;atomS=cpy.atomS;
 //cout << "This molecule has the formula " << this->get_formula() << " and
 //weight " << this->get_weight() << endl;
 }
 Molecule(const Atom& add){atomC =0; atomH = 0; atomS = 0;add_atoms(add);
// cout << "This molecule has the formula " << this->get_formula() << " and
//weight " << this->get_weight() << endl;
 }
 void add_atoms(Atom ,int num = 1);
 double get_weight();
 string get_formula();
 bool operator == (const Molecule& );
 Molecule& operator += (const Molecule& );
 Molecule& operator += (const Atom& );
 Molecule& operator= (const Molecule&);
 Molecule& operator= (const Atom&);
 Molecule& operator+(const Molecule&);
 Molecule& operator+(const Atom&);
 int is_stable();

static vector<string> splitstring(char flag, string s);
static Molecule& string_to_molecule(string str);//function that use in ecq
                                                //fucntion to transform string
                                                //into Molecule form.
string evaluate_chemical_equation(string str);

//Blow is the code that I tried to use template function but failed....

 /*
  template<typename T>
  Molecule operator + (T add1, T add2)
 {
     Hydrogenyx H;
     Carbonyx C;
     Sulphuryx S;
     Molecule add3(add1);
     add3.add_atoms(H,add2.atomH);
     add3.add_atoms(C,add2.atomC);
     add3.add_atoms(S,add2.atomS);
     return add3;
 }
 */
 /*
 template<typename T>
 Molecule operator * (T mul1, T mul2)
  {
      Molecule mul3=mul1;
      add3.add_atoms(add2);
      return add3;
  }
 */
};

#endif
