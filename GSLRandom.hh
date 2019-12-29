#ifndef __GSLRANDOM_HH__
#define __GSLRANDOM_HH__

#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

#include <cstdio>
#include <memory>
#include <stdexcept>
#define SHARED_PTR std::shared_ptr

#include <string>

/// A (partial) C++ interface to the GSL random number generators;
///
/// - to link remember to use -lgsl -lgslcblas
///
/// - names of generators are to be found at
///     http://www.gnu.org/software/gsl/manual/html_node/Random-number-generator-algorithms.html
///   & include gsl_rng_mt19937 (default) 
///             gsl_rng_ranlxs0 / s1 / s2 (second generation ranlux)
///             gsl_rng_ranlxd1 / d2      (48 bits, double precision output)
///             gsl_rng_ranlux / ranlux289 (orig. ranlux; 389: 24 decorr bits)
///             gsl_rng_cmrg (combined multiple recursive gen [l'Ecuyer])
///      etc.... (see docs for more info)
///
/// - some enquiry functions haven't been implemented (e.g. for state, range).
///
class GSLRandom {
public:
  /// create a default GSLRandom generator
  inline GSLRandom() {
    reset(gsl_rng_alloc(gsl_rng_mt19937));
  }
  /// create a default GSLRandom generator with the specified seed
  inline GSLRandom(unsigned long int s) {
    reset(gsl_rng_alloc(gsl_rng_mt19937));
    set(s);
  }
  /// create a GSLRandom generator with the specified generator type
  inline GSLRandom(const gsl_rng_type * T) {reset(gsl_rng_alloc(T));}

  /// create a GSLRandom generator with the specified generator type
  /// and seed
  inline GSLRandom(const gsl_rng_type * T, unsigned long int s) {
    reset(gsl_rng_alloc(T));
    set(s);
  }

  /// create a GSLRandom from a gsl_rng pointer (takes over ownership)
  inline GSLRandom(gsl_rng * T) {reset(T);}

  /// creates a GSLRandom as a clone of an existing generator
  /// (including its state)
  inline GSLRandom(const GSLRandom & orig) : GSLRandom(gsl_rng_clone(orig.gsl_generator())) {}

  /// copies the state of orig into this generator; this one
  /// must be of the same type as the original.
  inline void copy_state(const GSLRandom & orig) {
    gsl_rng_memcpy(gsl_generator(), orig.gsl_generator());
  }
  
  inline void reset(gsl_rng * r) {_r.reset(r, gsl_rng_free);}

  /// set the seed
  inline void set(unsigned long int s) {gsl_rng_set(_r.get(), s);}

  /// returns in range [0,1) (includes 0, excludes 1)
  inline double uniform() const {return gsl_rng_uniform(_r.get());}

  /// returns in range (0,1) (excludes 0, excludes 1)
  inline double uniform_pos() const {return gsl_rng_uniform_pos(_r.get());}

  /// returns in range [xmin,xmax) (includes xmin, excludes xmax -- modulo rounding)
  inline double uniform(double xmin, double xmax) const {return xmin + (xmax-xmin)*gsl_rng_uniform(_r.get());}

  /// returns +-1 with equal probability
  inline double sign() const {return uniform() < 0.5 ? -1 : 1;}

  /// returns a Gaussian distributed random number with standard
  /// deviation sigma and mean zero
  inline double gaussian(double sigma) const {return gsl_ran_gaussian(_r.get(),sigma);}
  /// Gaussian with sigma=1
  inline double gaussian() const {return gsl_ran_ugaussian(_r.get());}

  /// return random x according to p(x) dx = {1 \over \mu} \exp(-x/\mu) dx
  inline double exponential(double mu) const {return gsl_ran_exponential(_r.get(),mu);}

  // returns random k according to p(k) = {\mu^k \over k!} \exp(-\mu)
  inline unsigned int poisson(double mu) const {return gsl_ran_poisson(_r.get(),mu);}

  // returns a number generated according to the gamma distribution
  inline double gamma(double a, double b) const { return gsl_ran_gamma(_r.get(),a,b);}

  /// returns integer in range [0,n-1]
  inline unsigned long int uniform_int(unsigned long int n) const {
    return n > 1 ? gsl_rng_uniform_int(_r.get(), n) : 0;}

  /// returns true with probability prob
  inline bool accept(double prob) const {return uniform() < prob;}

  /// returns true with probability local_value/max_value
  inline bool accept(double local_value, double max_value) const {return max_value*uniform() < local_value;}

  /// returns true with probability local_value/max_value
  /// Before, checks that the local value does not exceed the max value (assertion)
  inline bool check_and_accept(double local_value, double max_value) const {
    assert(local_value <= max_value);
    return accept(local_value, max_value);
  }

  inline std::string name() const {return std::string(gsl_rng_name (_r.get()));}

  inline gsl_rng * gsl_generator() {return _r.get();}
  inline const gsl_rng * gsl_generator() const {return _r.get();}

  /// writes state of the generator to the specified filename
  inline void write_state_to_file(const std::string & filename) const {
    std::FILE * file = std::fopen(filename.c_str(), "wb");
    if (file == nullptr) throw std::runtime_error("Could not open "+filename);
    int result = gsl_rng_fwrite(file, gsl_generator());
    if (result != 0) throw std::runtime_error("Problem writing gsl state to "+filename);
    result = fclose(file);
    if (result != 0) throw std::runtime_error("Problem closing "+filename);
  }

  /// reads state of the generator from the specified filename
  inline void read_state_from_file(const std::string & filename) {
    std::FILE * file = std::fopen(filename.c_str(), "rb");
    if (file == nullptr) throw std::runtime_error("Could not open "+filename+" for reading");
    int result = gsl_rng_fread(file, gsl_generator());
    if (result != 0) throw std::runtime_error("Problem reading gsl state from "+filename);
    result = fclose(file);
    if (result != 0) throw std::runtime_error("Problem closing "+filename);
  }

  
  /// no specific destructor
  //inline ~GSLRandom() {}
protected:
  SHARED_PTR<gsl_rng> _r;
};

#endif // __GSLRANDOM_HH__
