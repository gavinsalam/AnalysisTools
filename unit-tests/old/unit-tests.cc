/// main program for unit testing.
/// NB1: done in such a way that it has an explicit main, so that
///      automated makefile generation picks it up
/// NB2: it is quite slow to compile (about 20s with -O2 and clang),
///      but need only be compiled once.
#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

int main( int argc, char* argv[] ) {
  // global setup...

  int result = Catch::Session().run( argc, argv );

  // global clean-up...

  return result;
}

