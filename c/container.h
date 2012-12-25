#ifndef CONTAINER_H
#define CONTAINER_H

template<typename T>
class array2{
public:
  array2(int n, int m, const T& value=T()):n_(n), m_(m), vec_(n*m, value){};
  inline T& get(int i,int j){return vec_[i*m_+j];}
private:
  int n_;
  int m_;
  std::vector<T> vec_;
};

template<typename T>
class array3{
public:
  array3(int n, int m, int l, const T& value=T()):n_(n), m_(m), l_(l), lm_(l*m), vec_(n*m*l, value){};
  inline T& get(int i, int j, int k){return vec_[lm_*i + l_*j + k];}
private:
  int n_;
  int m_;
  int l_;
  int lm_;
  std::vector<T> vec_;
};

#endif