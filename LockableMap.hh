#include<map>
#include <iostream>

//----------------------------------------------------------------------
/// a map whose normal access functions can be locked, so that
/// they return an error when you try to access something that
/// does not exist.
///
/// Note that you can still modify the contents associated with a key.
template<class T, class U>
class LockableMap : public std::map<T,U> {
public:
  LockableMap() : _locked(false) {}
  
  const U & operator[](const T & key) const {
    typename std::map< T, U >::const_iterator location = this->find(key);
    if (location == this->end()) throw KeyNotPresentForConstAccess(key);
    return location->second;
  }

  /// return access to a key if possible
  U & operator[](const T & key) {
    typename std::map< T, U >::iterator location = this->find(key);
    if (location == this->end()) {
      if   (_locked) throw InvalidAccessToLockedMap(key);
      else return std::map<T,U>::operator[](key);
    } else {
      return location->second;
    }
  }

  /// lock the map 
  void lock()   {_locked = true;}
  /// unlock the map
  void unlock() {_locked = false;}

  class InvalidAccessToLockedMap {
  public:
    InvalidAccessToLockedMap() {}
    InvalidAccessToLockedMap(const T & key) {std::cerr << "invalid write access with key " << key << std::endl;}
  };
  class KeyNotPresentForConstAccess {
    KeyNotPresentForConstAccess() {}
    KeyNotPresentForConstAccess(const T & key) {std::cerr << "key not present for const access " << key << std::endl;}
  };
private:
  bool _locked;
};
