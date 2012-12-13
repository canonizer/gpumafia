#ifndef MAFIA_REF_H_
#define MAFIA_REF_H_
/** @file ref.h simple counting pointer implementation */

/** a simple counting pointer; the referenced object must have a field
    to store usage count
 */
template<class T>
class ref {
 private:
  void inc(const ref<T>& r) {
    ptr = r.ptr;
    if(ptr)
      ptr->ref_count++;
  }
	void inc(T *p) {
		ptr = p;
		if(ptr)
			ptr->ref_count++;
	}
  void dec() {
    if(ptr) {
      ptr->ref_count--;
      if(!ptr->ref_count)
	delete ptr;
    }
  }
 public:
  ref(T *p = 0) { inc(p); }
  ref(const ref<T>& r) { inc(r); }
  T* operator->() const { return ptr; }
  operator T*() const { return ptr; }
  T& operator*() const {return *ptr; }
  ~ref() { dec(); }
  ref<T> & operator=(const ref<T> &r) {
    if(this != &r) {
      dec();
      inc(r);
    }
    return *this;
  }
protected:
  T *ptr;
}; // class ref

/** a class which can be referenced using ref */
class ref_object {
 public:
  /** default constructor 
      @note reference count set to 1 on object creation
   */
  ref_object() {ref_count = 0;}
  /** copy constructor 
   @note reference count of a copy is set to 1; the object is
   typically not supposed to be copied, however
  */
  ref_object(const ref_object& o) {
    ref_count = 1;
  }
  /** just a virtual destructor, nothing to do with reference counts */
  virtual ~ref_object() {}
  // TODO: make it private, do friendship correctly
 public:
  int ref_count;
  friend class ref<ref_object>;
};

#endif
