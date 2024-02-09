#include <iostream>

#include <map>
#include <functional>

/** Memoization functor to wrap an expensive unary function
 * Note: Requires op<(Arg,Arg)
 * Can be generalized to arbitrary number of arguments
 */
template <typename Result, typename Arg, typename Function>
struct Memoizer {
  Memoizer(Function& func)
      : func_(func) {
  }

  Result operator()(Arg arg) {
    // Use equal_range to prevent two map lookups
    auto range = cache_.equal_range(arg);

    // If the range is empty, the argument was not found
    if (range.first == range.second) {
      std::cout << "Computing" << std::endl;
      // Insert the (arg,result) pair in the correct positon (no search)
      range.first = cache_.insert(range.first, std::make_pair(arg, func_(arg)));
    } else {
      std::cout << "Found!" << std::endl;
    }

    // Return the map's value for this argument
    return range.first->second;
  }

 private:
  std::map<Arg, Result> cache_;
  Function& func_;
};

/** Construct a Memoizer functor for the function */
template <typename Result, typename Arg, typename Function>
Memoizer<Result, Arg, Function> memoize(Function func) {
  return Memoizer<Result, Arg, Function>(func);
}

/** Construct an anonymous memoizer closure for the function */
template <typename Result, typename Arg, typename Function>
std::function<Result(Arg)> memoize2(Function func_) {
  std::map<Arg, Result> cache_;
  return [=](Arg arg) mutable {
    // Use equal_range to prevent two map lookups
    auto range = cache_.equal_range(arg);

    // If the range is empty, the argument was not found
    if (range.first == range.second) {
      std::cout << "Computing" << std::endl;
      // Insert the (arg,result) pair in the correct positon (no search)
      range.first = cache_.insert(range.first, std::make_pair(arg, func_(arg)));
    } else {
      std::cout << "Found!" << std::endl;
    }

    // Return the map's value for this argument
    return range.first->second;
  };
}


/** Test the memoization.
 * Function is a function from int to int
 */
template <typename Function>
void test_int_func(Function& f) {
  auto ff = memoize2<int,int>(f);
  //auto ff = memoize2<int,int>(f);

  std::cout << ff(3) << std::endl;
  std::cout << ff(4) << std::endl;
  std::cout << ff(5) << std::endl;
  std::cout << ff(6) << std::endl;
  std::cout << ff(5) << std::endl;
  std::cout << ff(4) << std::endl;
  std::cout << ff(4) << std::endl;
  std::cout << std::endl;
}


struct MyFunctor {
  int operator()(int x) {
    return x+1;
  }
};

int addone(int x) {
  return x+1;
}

int main() {
  MyFunctor f;
  test_int_func(f);

  test_int_func(addone);

  auto my_lambda = [](int x) { return x+1; };
  test_int_func(my_lambda);
}
