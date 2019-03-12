// program to test speed of using map v. unordered_map
// 
// 10^8 events, with 2 (unordered) map entries, whose first letters
// differ (about 15 characters each)
//
//  | unordered_map | 6.4s | 32ns/call | 
//  | map           | 5.7s | 29ns/call |
//
// 10^8 events, with 2 (unordered) map entries, whose last letters
// differ (about 15 characters each)
// 
//  | unordered_map | 6.3s |
//  | map           | 7.1s |
// 
// 10^8 events, with 2 (unordered) map entries, whose last letters
// differ (about 47 characters each)
//
//  | unordered_map | 26.7s | 134ns/call |
//  | map           | 27.8s | 139ns/call |
//
// 10^8 events, with 2 (unordered) map entries, whose first letters
// differ (about 47 characters each)
//
//  | unordered_map | 19.9s | 100ns/call |
//  | map           | 20.7s | 104ns/call  |

// 10^8 events, with 2 (unordered) map entries, whose first letters
// differ (about 47 characters each), using a reference to the entries
// in the main loop
//
//  | unordered_map | 0.3s | 1.5ns/call |
//  | map           | 0.3s | 1.5ns/call  |

#include <map>
#include <unordered_map>
#include <iostream>
#include <string>
#include "CmdLine.hh"

using namespace std;

class DefaultDouble {
public:
  DefaultDouble() : value(0.0) {}
  double value;
};

class TTT {
public:
  map<string,DefaultDouble> repo;
  decltype(repo)::mapped_type & a = repo["hello"];
};

int main(int argc, char ** argv) {
  CmdLine cmdline(argc,argv);
  unsigned int nev = cmdline.value("-nev", 1e4);
  //unordered_map<string, DefaultDouble> repo;
  map<string, DefaultDouble> repo;
  auto & a = repo["atest-string-mod-test-string-mod-test-string-mod"];
  auto & b = repo["btest-string-mod-test-string-mod-test-string-mod"];
  for (unsigned iev = 0; iev < nev; iev++) {
    a.value += 0.4;
    b.value += 0.2;
    // repo["atest-string-mod-test-string-mod-test-string-mod"].value += 0.4;
    // repo["btest-string-mod-test-string-mod-test-string-mod"].value += 0.2; 
  }
  for (const auto & a: repo) {
    cout << a.first << " " << a.second.value << endl;
  }
  // cout << repo["test-string-mod"] << endl;
  // cout << repo["test-string-mod2"] << endl;
}
