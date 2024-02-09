/*
  This file demonstrates the Proxy design pattern.
  It contains two classes: a "SimpleSet" and a "SimpleElement".

  When the SimpleSet changes size, all of its elements are copied to new
  locations in memory and renumbered. But that's OK: since the SimpleElements
  are proxies, they always access the SimpleSet using the most up-to-date
  memory address and index.

  We use dynamic memory allocation here to make the issues obvious; you would
  almost certainly use an STL container. Because an STL container like vector 
  is rarely re-allocated, and as a result, faster for adding new elements. 
  Our proxy is also slow: it takes O(N)
  time to access an element, where N is the number of elements. You will fix
  this in your Graph. By using STL functions?

*/

#include <iostream>
#include <assert.h>

class SimpleSet {
  typedef unsigned size_type;	//create an alias, size_type, for unsigned

  struct internal_element {
    const char *text;   // The text held by an element
    size_type uid;      // The unique identification number of an element
  };

  internal_element* elements_;
  size_type size_;
  size_type next_uid_;

 public:
 
  SimpleSet()
      : elements_(), size_(0), next_uid_(0) {	
  }

  ~SimpleSet() {
    delete[] elements_;	//delete pointer
  }

  /** The proxy class */
  class SimpleElement {
   public:
   
    SimpleElement() {
    }

    const char* text() const {	// const after text(): does not change class members
      return fetch().text;		//  const char* return a pointer to a constant character
    }

    void set_text(const char *text) {
      fetch().text = text;
    }
    
   private:
    SimpleSet* set_;	 // Pointer back to the SimpleSet container
    size_type uid_;		// This element's unique identification number

    SimpleElement(const SimpleSet* set, size_t uid)
        : set_(const_cast<SimpleSet*>(set)), uid_(uid) {	// proxy: delay creation of 
        // SimpleSet*
    }

    internal_element& fetch() const {
      for (size_type i = 0; i < set_->size(); ++i)	// replace set_->size() with size_ results in error
        if (set_->elements_[i].uid == uid_)
          return set_->elements_[i];
      assert(0);	
    }
    
    friend class SimpleSet;	// Allow SimpleSet to access SimpleElement's private members
  };

  /** Return SimpleSet's size. */
  size_type size() const {
    return size_;
  }
  
  /** Return a proxy object for element @a i. */
  SimpleElement operator[](size_t i) const {	// overload []
    assert(i < size());
    return SimpleElement(this, i);
  }
  
  /** Add a new element at the end.
   * @return A proxy for the new element */
  SimpleElement push_back(const char* text) {
    internal_element* new_elements = new internal_element[size_ + 1];
    for (size_t i = 0; i < size_; ++i)
      new_elements[i] = elements_[i];
    new_elements[size_].text = text;
    new_elements[size_].uid = next_uid_;
    delete[] elements_;
    elements_ = new_elements;
    ++size_;
    ++next_uid_;
    return SimpleElement(this, next_uid_ - 1);
  }
  
  /** Remove the element at position @a i, moving later elements down. */
  void remove(size_t i) {
    assert(i < size());
    for (++i; i < size(); ++i)
      elements_[i - 1] = elements_[i];	// shift element down 1
    --size_;		// re-define size_;
  }
  
  friend class SimpleElement;	// proxy (SimpleElement) needs permission to access
  // private state
  SimpleSet(const SimpleSet&) = delete;	// explicitly delete constructor SimpleSet
  SimpleSet& operator=(const SimpleSet&) = delete;	// explicitly delete member function =
};


int main() {
  SimpleSet v;
  SimpleSet::SimpleElement e0 = v.push_back("Hello");
  SimpleSet::SimpleElement e1 = v.push_back("World");
  std::cerr << e0.text() << " " << e1.text() << std::endl;
  // prints "Hello World"

  SimpleSet::SimpleElement e0_copy = v[0];
  e0.set_text("Goodbye");
  std::cerr << e0.text() << " " << e0_copy.text() << std::endl;
  // prints "Goodbye Goodbye": since SimpleElement is a proxy, e0 and e0_copy
  // both return the most up-to-date information

  SimpleSet::SimpleElement e2 = v.push_back("Friends");
  v.remove(1);
  std::cerr << e0.text() << " " << e2.text() << std::endl;
  // prints "Goodbye Friends": SimpleElement locates its element using a
  // unique number that stays stable even after SimpleSet's internal array
  // is rearranged
}
