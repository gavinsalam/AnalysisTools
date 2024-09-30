#ifndef __GSLRANDOM_HH__
#define __GSLRANDOM_HH__

#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

#include <cstdio>
#include <memory>
#include <stdexcept>
#include <cassert>
// for generator state output as string
#include <sstream>
#include <iomanip>
#define SHARED_PTR std::shared_ptr

#include <string>
#include <functional>

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
  inline GSLRandom(const gsl_rng_type * T) {
    reset(gsl_rng_alloc(T));
  }

  /// create a GSLRandom generator with the specified generator type
  /// and seed
  inline GSLRandom(const gsl_rng_type * T, unsigned long int s) {
    reset(gsl_rng_alloc(T));
    _seed = s;
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
  
  /// reset the generator to a new gsl_rng object
  inline void reset(gsl_rng * r) {
    _r.reset(r, gsl_rng_free);
    _update_granularity();
  }

  /// set the seed
  inline void set(unsigned long int s) {
    _seed = s;
    gsl_rng_set(_r.get(), s);
  }

  /// return the maximum value that can be generated
  unsigned long int max() const {return gsl_rng_max(_r.get());}

  /// return half the granularity (as relevant in the double range 0,1)
  double half_granularity() const {return _half_granularity;}

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

  /// returns an integer in the range [lo, hi-1]. 
  /// NB: hi must be >= lo
  inline long int uniform_int(int lo, int hi) const {
    unsigned long int n = hi - lo;
    return lo + uniform_int(n);
  }

  /// returns true with probability prob
  inline bool accept(double prob) const {return uniform() < prob;}

  /// returns true with probability local_value/max_value
  inline bool accept(double local_value, double max_value) const {return max_value*uniform() < local_value;}

  /// returns true with probability prob. If prob is 1, accept without
  /// generating a random number
  inline bool accept_lazy(double prob) const {return ((prob==1.0) || (accept(prob))); }

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

  /// returns a string with a hexadecimal string representation of the
  /// state of the generator
  inline std::string hex_state() {
    void * state = gsl_rng_state(gsl_generator());
    std::size_t n = gsl_rng_size(gsl_generator());
    char * state_char = (char *) state;
    std::ostringstream ostr;
    for (std::size_t i = 0; i < n; i++) {
      ostr << std::setfill('0') << std::setw(2) << std::hex << (0xff & (uint8_t) state_char[i]);
    }
    return ostr.str();
  }

  /// returns a string with a C++ hash of the hexadecimal string
  /// representation of the state of the generator. 
  ///
  /// Note: the hash is whatever C++ uses for its unordered maps, and is
  /// not guaranteed to be free of collisions. It's mainly intended to
  /// be helpful in debugging, e.g. to have a quick, rough check of
  /// whether the generator is in the same state across two different
  /// processes.
  inline std::size_t hash_hex_state() {
    return std::hash<std::string>{}(hex_state());
  }
  
  unsigned long int seed() const {return _seed;}

  std::string description() const {
    std::ostringstream ostr;
    ostr << "GSLRandom with rng = " << name() << " with seed = " << seed();
    return ostr.str();
  }

  inline static unsigned long int default_seed() {return gsl_rng_default_seed;}

  /// no specific destructor
  //inline ~GSLRandom() {}
protected:
  SHARED_PTR<gsl_rng> _r;
  unsigned long int _seed = gsl_rng_default_seed;
  double _half_granularity;

  void _update_granularity() {
    _half_granularity = 0.5 / double(max());
  }
};

#endif // __GSLRANDOM_HH__
