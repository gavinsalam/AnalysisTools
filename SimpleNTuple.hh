#ifndef   __SIMPLENTUPLE_HH__
#define   __SIMPLENTUPLE_HH__

#include <iostream>
#include <sstream>
#include <string>


/// a very simple (not very efficient) class to write ntuples to a file
class SimpleNTupleLine {
public:
  /// start off an entry
  SimpleNTupleLine(std::ostream * ostr = NULL, const std::string & title = "") 
    : _ostr(ostr) {
    if (_ostr) _sstr << "!!" << title;
  }
  SimpleNTupleLine(std::ostream & ostr, const std::string & title = "") 
    : _ostr(&ostr) {
    if (_ostr) _sstr << "!!" << title;
  }

  /// add a member entry for this line of the ntuple
  template <class T>
  void add(const std::string & name, const T & value) {add(name.c_str(), value);}

  /// the same, straight from a C string
  template <class T>
  void add(const char * name, const T & value) {
    if (_ostr) _sstr << " " << name << "#" << value;
  }

  void add(const std::string & name, double value, int prec) {add(name.c_str(), value, prec);}
  
  void add(const char * name, double value, int prec) {
    if (_ostr) {
      int old_prec = _sstr.precision();   // store old precision
      _sstr.precision(prec);              // set it
      _sstr << " " << name << "#" << value;
      _sstr.precision(old_prec);          // reset it to the old value
    }
  }
    
  /// finish off the line once we're done
  ~SimpleNTupleLine() {
    if (_ostr) *_ostr << _sstr.str() << std::endl;
  }

private:
  std::ostream * _ostr;
  std::ostringstream _sstr;
};


#endif // __SIMPLENTUPLE_HH__
