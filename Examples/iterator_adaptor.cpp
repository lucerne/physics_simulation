
#include <iterator>

/** @class iterator_adaptor
 * @brief A simple adaptor class for forward iterators
 * Classes may inherit and override appropriate operators
 */
template <typename Derived, typename Base, typename Value>
class iterator_adaptor {
  typedef iterator_adaptor<Derived,Base,Value> this_type;
  typedef Derived derived_type;
  typedef Base base_type;

 public:
  typedef Value value_type;
  typedef value_type* pointer;
  typedef value_type& reference;
  typedef std::forward_iterator_tag iterator_category;
  typedef std::ptrdiff_t difference_type;

  // For convenience in derived classes
  typedef this_type super_type;

  // Constructors
  explicit iterator_adaptor(const base_type& b)
      : b_(b) {
  }
  iterator_adaptor(const this_type& m)
      : b_(m.base()) {
  }

  // Inheriting classes use these to access the underlying "iterator"
  base_type& base() {
    return b_;
  }
  const base_type& base() const {
    return b_;
  }
  derived_type& derived() {
    return static_cast<derived_type&>(*this);
  }
  const derived_type& derived() const {
    return static_cast<const derived_type&>(*this);
  }

  // Overloadable forward iterator iterator operators //

  derived_type& operator++() {                 // Prefix increment (++x)
    ++base();
    return derived();
  }
  derived_type operator++(int) {               // Postfix increment (x++)
    derived_type tmp = derived();
    ++base();
    return tmp;
  }
  reference operator*() {                      // Dereference (*x)
    return *base();
  }
  const reference operator*() const {          // Const Dereference (*x)
    return *base();
  }
  pointer operator->() {                       // Ptr to Member (x->member)
    placeholder_ = operator*();
    return &placeholder_;
  }
  bool operator==(const this_type& m) const {  // Equality (x==y)
    return base() == m.base();
  }
  bool operator!=(const this_type& m) const {  // Not Equality (x!=y)
    return base() != m.base();
  }
  derived_type& operator=(const this_type& m) { // Assignment (x=y)
    base() = m.base();
    return derived();
  }
 private:
  mutable value_type placeholder_;
  base_type b_;
};


#include <vector>

// Derived class implementing a transformed return value
struct my_transform_iterator
    : public iterator_adaptor<my_transform_iterator,
                              std::vector<int>::iterator,
                              int>
{
  my_transform_iterator(std::vector<int>::iterator _m)
      : super_type(_m) {
  }
  int operator*() {
    return (*base()) * (*base());
  }
};


// Derived class implementing a strided iterator
struct my_stride_iterator
    : public iterator_adaptor<my_stride_iterator,
                              std::vector<int>::iterator,
                              int>
{
  my_stride_iterator(std::vector<int>::iterator _m)
      : super_type(_m) {
  }
  my_stride_iterator& operator++() {
    base() += 2;
    return *this;
  }
};

// Derived class implementing an iterator with an index into container
struct my_position_iterator
    : public iterator_adaptor<my_position_iterator,
                              int,
                              int>
{
  my_position_iterator(const std::vector<int>& _c, int _pos)
      : super_type(_pos), c(_c) {
  }
  int operator*() {
    return c[base()];
  }
  const std::vector<int>& c;
};


#include <algorithm>
#include <iostream>

int main() {
  std::vector<int> a(4,2);
  std::partial_sum(a.begin(), a.end(), a.begin());
  std::copy(a.begin(), a.end(), std::ostream_iterator<int>(std::cout, ", "));
  std::cout << std::endl;

  auto begin = a.begin();
  auto end = a.end();

  std::cout << "Normal sum: " << std::accumulate(begin, end, 0) << std::endl;

  auto pbegin = my_position_iterator(a, 0);
  auto pend = my_position_iterator(a, 4);

  std::cout << "By position: " << std::accumulate(pbegin, pend, 0) << std::endl;

  auto tbegin = my_transform_iterator(a.begin());
  auto tend = my_transform_iterator(a.end());

  std::cout << "Squared sum: " << std::accumulate(tbegin, tend, 0) << std::endl;

  auto sbegin = my_stride_iterator(a.begin());
  auto send = my_stride_iterator(a.end());

  std::cout << "Strided sum: " << std::accumulate(sbegin, send, 0) << std::endl;

  return 0;
}
