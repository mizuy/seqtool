#ifndef STORE_H
#define STORE_H

#include <shared_ptr>
#include <list>
#include "container.h"

/*
template<typename T>
struct sp_less {
public:
	bool operator()(shared_ptr<T>& lt, shared_ptr<T>& rt){
		return static_cast<float>(*lt) > static_cast<float>(*rt);
	}
};
*/


/*
template<typename T>
class TopHolder {
public:
	typedef std::vector<shared_ptr<T> > container;
	typedef container::iterator iterator;

	TopHolder() : max_length_(0){
	}
	void resize(int len){
		max_length_ = len;
		list_.reserve(max_length);
	}
	bool add(shared_ptr<T> value){
		if(list_.size() < max_length_){
			list_.push_back(value);
			lowest_ = min_element();
			return true;
		}
		else if(*lowest_ < *value){
			//list_.erase(lowest_);
			//list_.push_back(value);
			*lowest_ = value;
			lowest_ = min_element();
			return true;
		}
		return false;
	}
	s
	iterator begin(){
		return list_.begin();
	}
	iterator end(){
		return list_.end();
	}
	iterator min_element(){
		assert(list_.size() > 0);
		if(list_.size()==1){
			return list_.begin();
		}
		else{
			iterator lowest = list_.begin();
			for(iterator i=lowest; i!=list_.end(); i++){
				if((*i) < (*lowest))
					lowest = i;
			}
			return lowest;
		}
	}
private:
	container list_;
	iterator lowest_;
	int max_length_;
};

template<typename T>
class Store{
public:
	Store(int length, int chunk, int coverage);

	void set(int i, int j, shared_ptr<T> value);
private:
	vector<TopHolder<T> > queue_;
	int length_;
	int clength_;
	int converage_;
};
*/
#endif