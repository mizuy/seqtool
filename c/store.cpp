#include "store.h"

Store::Store(int length, int chunk, int coverage)
	: length_(length), clength_(length/chunk + (length%chunk)?1:0), converage_(converage){
	queue_.resize(clength_);
	for(auto i=queue_.begin(); i!=queue_.end(); i++){
		i->resize(coverage);
	}
}

void Store::set(int i, int j, shared_ptr<T> value){
	int ii = i/chunk_;
	int jj = j/chunk_;
	for(int k=ii+1; k<jj; k++){
		queue[k].add(value);
	}
}
