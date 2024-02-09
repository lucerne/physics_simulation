#include <stddef.h>
#include <iostream>
#include <vector>
#include <assert.h>

#define USE_EXPR_TEMP 1

template <typename E>
struct MatExpr {
  inline const E& derived() const {
    return static_cast<const E&>(*this);
  }
  inline double operator()(int i, int j) const {
    return derived()(i,j);
  }
  inline size_t dim(int i) const {
    return derived().dim(i);
  }
};


template <typename E1, typename E2>
struct MatAdd : public MatExpr<MatAdd<E1, E2>> {
  MatAdd(const E1& a, const E2& b)
      : a_(a), b_(b) {
    assert(a_.dim(0) == b.dim(0) && a.dim(1) == b.dim(1));
  }
  size_t dim(int i) const {
    return a_.dim(i);
  }
  double operator()(int i, int j) const {
    return a_(i, j) + b_(i, j);
  }
private:
  const E1& a_;
  const E2& b_;
};


template <typename E1, typename E2>
struct MatMult : public MatExpr<MatMult<E1, E2>> {
  MatMult(const E1& a, const E2& b)
      : a_(a), b_(b) {
    assert(a_.dim(1) == b.dim(0));
  }
  size_t dim(int i) const {
    return (i == 0 ? a_.dim(0) : b_.dim(1));
  }
  double operator()(int i, int j) const {
    double result = 0;
    for (size_t k = 0; k < a_.dim(1); ++k)
      result += a_(i, k) * b_(k, j);
    return result;
  }
private:
  const E1& a_;
  const E2& b_;
};

template <typename E>
struct MatScale : public MatExpr<MatScale<E>> {
  MatScale(double s, const E& a)
      : s_(s), a_(a) {
  }
  size_t dim(int i) const {
    return a_.dim(i);
  }
  double operator()(int i, int j) const {
    return s_ * a_(i,j);
  }
private:
  double s_;
  const E& a_;
};


struct Matrix : public MatExpr<Matrix> {
  Matrix(size_t n, size_t m, double v = 0)
      : dim_{n, m}, v_(n * m, v) {
      }
  Matrix(const Matrix& x)
      : dim_{x.dim_[0], x.dim_[1]}, v_(x.v_) {
      }
  template <typename E>
  Matrix(const MatExpr<E>& x) {
    dim_[0] = x.dim(0);
    dim_[1] = x.dim(1);
    v_.resize(dim_[0] * dim_[1]);
    for (size_t i = 0; i < dim_[0]; ++i)
      for (size_t j = 0; j < dim_[1]; ++j)
        (*this)(i, j) = x(i, j);
  }
  size_t dim(int i) const {
    assert(i == 0 || i == 1);
    return dim_[i];
  }
  double operator()(int i, int j) const {
    return v_[i*dim_[0] + j];
  }
  double& operator()(int i, int j) {
    return v_[i*dim_[0] + j];
  }
  friend std::ostream& operator<<(std::ostream& os, const Matrix& m) {
    for (size_t i = 0; i < m.dim(0); ++i) {
      for (size_t j = 0; j < m.dim(1); ++j)
        os << m(i, j) << "\t";
      os << "\n";
    }
    return os;
  }
 private:
  size_t dim_[2];
  std::vector<double> v_;
};


#if USE_EXPR_TEMP

template <typename E1, typename E2>
MatAdd<E1, E2> operator+(const MatExpr<E1>& a, const MatExpr<E2>& b) {
  return MatAdd<E1, E2>(a.derived(), b.derived());
}

template <typename E1, typename E2>
MatMult<E1, E2> operator*(const MatExpr<E1>& a, const MatExpr<E2>& b) {
  return MatMult<E1, E2>(a.derived(), b.derived());
}

template <typename E>
MatScale<E> operator*(double s, const MatExpr<E>& a) {
  return MatScale<E>(s, a.derived());
}

#else

Matrix operator+(const Matrix& a, const Matrix& b) {
  Matrix r(a.dim(0), b.dim(1));
  for (size_t i = 0; i < r.dim(0); ++i)
    for (size_t j = 0; j < r.dim(1); ++j)
      r(i, j) = a(i, j) + b(i, j);
  return r;
}

Matrix operator*(const Matrix& a, const Matrix& b) {
  Matrix r(a.dim(0), b.dim(1));
  for (size_t i = 0; i < r.dim(0); ++i)
    for (size_t j = 0; j < r.dim(1); ++j)
      for (size_t k = 0; k < a.dim(1); ++k)
        r(i, j) += a(i, k) * b(k, j);
  return r;
}

Matrix operator*(double s, const Matrix& a) {
  Matrix r(a.dim(0), a.dim(1));
  for (size_t i = 0; i < r.dim(0); ++i)
    for (size_t j = 0; j < r.dim(1); ++j)
      r(i, j) = s * a(i,j);
  return r;
}

#endif

#include <sys/time.h>
#include <unistd.h>
/** Clock class, useful when timing code.
 */
class Clock {
 public:
  /** Construct a Clock and start timing. */
  Clock() {
    start();
  }
  /** Start the clock. */
  inline void start() {
    time_ = now();
  }
  /** Return the amount of time elapsed since the last start. */
  inline double elapsed() const {
    timeval tv = now();
    timersub(&tv, &time_, &tv);
    return tv.tv_sec + tv.tv_usec/1e6;
  }
 private:
  timeval time_;
  inline static timeval now() {
    timeval tv;
    gettimeofday(&tv, 0);
    return tv;
  }
};



int main() {
  int N = 512;
  Matrix A(N,N,1);
  Matrix B(N,N,1);
  Clock c;

  {
    c.start();
    Matrix C = 2*A + B + A + B + A + B + A + B + A + B + B + A;
    double time = c.elapsed();
    std::cout << time << "\t" << C(0,0) << std::endl;
  }
}
