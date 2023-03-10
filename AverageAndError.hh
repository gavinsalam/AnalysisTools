#ifndef __AVERAGEANDERROR_HH__
#define __AVERAGEANDERROR_HH__

#include<cmath>
#include<vector>
#include<limits>

/// micro class to calculate averages and errors
class AverageAndError {
public:

  /// default constructor
  AverageAndError() { _sum = 0.0; _sum2 = 0.0; _sum3 = 0.0; _sum4 = 0.0; _n=0;}
  
  /// add one entry
  inline virtual void add(double x) { _sum += x;
                              double x2 = x*x; 
                              _sum2 += x2;
                              _sum3 += x2*x;
                              _sum4 += x2*x2;
                              _n += 1;
                            }

  /// add vector with entries
  inline void add(std::vector<double> v)  { 
      for(unsigned i=0;i<v.size();i++) add(v[i]);
  }
  
  /// alternative way to add one or more entries, for consistency with 
  /// SimpleHist, AveragingHist, etc
  inline void add_entry(double x) { add(x); }
  inline void add_entry(std::vector<double> v) { add(v); }

  /// add one event with a different notation
  inline void operator+= (double x) { add(x); }

  /// return sum
  inline double sum() const { return _sum; }

  /// return sum2, second way for consistency with AveragingHist
  inline double sum2() const { return _sum2; }
  inline double sum_of_squares() const { return _sum2; }

  /// return sum3
  inline double sum3() const { return _sum3; }
  
  /// return sum4
  inline double sum4() const { return _sum4; }
  
  /// return number of events
  inline long long int n() const { return _n; }
  /// alternative way to return number of events, for consistency with
  /// SimpleHist, AveragingHist, etc
  inline long long int n_entries() const { return n(); }
  
  /// allow the user to reset the effective value of n
  inline void set_n(long long int n_in) { _n = n_in;}

  /// calculate and return average
  inline double average() const { return (_n > 0) ? _sum/_n : 0. ; }

  /// calculate and return average of squares, second way for consistency with AveragingHist
  inline double average2() const { return (_n > 0) ? _sum2/_n : 0. ; }
  inline double average_of_squares() const { return average2(); }

  /// calculate and return error
  inline double error() const { return sd()/std::sqrt(_n); }

  /// return the error on the sum (as opposed to the error on the average)
  inline double error_on_sum() { return error() * n(); }

  /// calculate and return the unbiased sample variance
  inline double variance() const { return (_n > 1) ? std::abs(_sum2 - _sum*_sum/_n)/(_n-1.0) : 0.; }

  /// calculate and return the standard deviation (i.e. sqrt of the unbiased sample variance)
  inline double sd() const { return (_n > 1) ? std::sqrt(variance()) : 0.; }

  /// calculate and return variance of the (unbiased) variance, i.e.
  /// var[S^2] = 1/n ( E[(X-E[X])^4] - (n-3)/(n-1) E[(X-E[X])^2]^2 )
  /// (see also http://www.talkstats.com/showthread.php/12302-Standard-error-of-the-sample-standard-deviation)
  inline double variance_of_variance() const {
      return  (_n > 1) ? 
               ( _sum4/_n - 4*_sum3*_sum/powd(_n,2) + ((3-_n)*_sum2*_sum2/powd(_n,2) +
                 4*(2*_n-3)*_sum2*_sum*_sum/powd(_n,3) + 2*(3-2*_n)*pow(_sum,4)/powd(_n,4) ) / (_n-1) )/_n 
	      : 0.; 
  }
  
  /// return error on variance (i.e. the square root of variance_of_variance() )
  inline double error_on_variance() const  { return (_n > 1) ? std::sqrt(variance_of_variance()) : 0.; }

  /// return error on standard deviation, given in approximate form as error of sqrt of variance
  inline double error_on_sd() const  { return (_n > 1) ? error_on_variance()/sd()/2. : 0.; }

  /// virtual dummy destructor
  virtual ~AverageAndError() {};
  
private:  
  double _sum, _sum2, _sum3, _sum4;

  // make this long long int so that we can handle subtractions safely
  // and at the same time go to large numbers of events
  long long int _n;

  /// a pow operation that converts a long long int to a double to avoid
  /// the risk of overflow 
  double powd(long long int i, int p) const {return std::pow(1.0*i,p);}

  // making this friend avoid "double-definition" compilation isses
  // (see e.g. https://en.cppreference.com/w/cpp/language/friend)
  friend std::ostream & operator<<(std::ostream & ostr, const AverageAndError & avg) {
    ostr << avg.average() << " +- " << avg.error() << " (n= " << avg.n() << " )";
    return ostr;
  }

};

/// adds functionality of maintaining info about the smallest and largest values
class AverageAndErrorWithBounds : public AverageAndError {
public:
  // default constructor will do the job here
  inline void add(double x) {
    AverageAndError::add(x);
    if (x < _min) _min = x;
    if (x > _max) _max = x;
  }

  /// return the the minimum value encountered
  double min() const {return _min;}
  /// return the the maximum value encountered
  double max() const {return _max;}
  
private:
  /// internal store of the bounds
  /// (note the _min starts as numeric_limits<double>::max() because then any thing smaller
  /// than that will immediately register as the minimum)
  double _min = std::numeric_limits<double>::max(), _max=-std::numeric_limits<double>::max();
};

#endif // __AVERAGEANDERROR_HH__
