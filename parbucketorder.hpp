/*
 *  BEANDisco: parallel bucket order definition
 *  
 *  Copyright 2011 Teppo Niinimäki <teppo.niinimaki(at)helsinki.fi>
 *  
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "stacksubset.hpp"
#include "common.hpp"

#ifndef PARBUCKETORDER_HPP
#define PARBUCKETORDER_HPP


class ParBucketOrderFamily {
	int n_;
	int nChains_;
	int chainLength_;
	int bucketSize_;
	int* order_;

	int* chainSizes_;
	int* nBuckets_;
	int** bucketSizes_;
	
	size_t* nChainIdeals_;
	size_t nIdeals_;
	
	//int* chainNums_;
	//int* bucketNums_;
	//int* itemNums_;
	
	ParBucketOrderFamily(const ParBucketOrderFamily&); // disable copy constructor
	ParBucketOrderFamily& operator=(const ParBucketOrderFamily&); // disable copying
public:
	const int n;
	const int nChains;
	const int maxBucketSize;
	
	ParBucketOrderFamily(int _n, int _maxBucketSize, int _nChains) :
			n(_n), nChains(_nChains), maxBucketSize(_maxBucketSize) {
		assert(1 <= maxBucketSize && maxBucketSize < 32);
		assert(1 <= n);
		
		chainSizes_ = new int[nChains];
		nBuckets_ = new int[nChains];
		bucketSizes_ = new int*[nChains];
		for (int c = 0; c < nChains; ++c) {
			chainSizes_[c] = (n + c) / nChains;
			nBuckets_[c] = (chainSizes_[c] + maxBucketSize - 1) / maxBucketSize;
			bucketSizes_[c] = new int[nBuckets_[c]];
			for (int b = 0; b < nBuckets_[c]; ++b)
				bucketSizes_[c][b] = (chainSizes_[c] + b) / nBuckets_[c];
		}
		
		nChainIdeals_ = new size_t[nChains];
		nIdeals_ = 1;
		for (int c = 0; c < nChains; ++c) {
			nChainIdeals_[c] = 1;
			for (int b = 0; b < nBuckets(c); ++b)
				nChainIdeals_[c] += (1 << bucketSize(c, b)) - 1;
			nIdeals_ *= nChainIdeals_[c];
		}
	}
	
	~ParBucketOrderFamily() {
		delete[] chainSizes_;
		delete[] nBuckets_;
		for (int c = 0; c < nChains; ++c)
			delete[] bucketSizes_[c];
		delete[] bucketSizes_;

		delete[] nChainIdeals_;
	}
	
	int chainSize(int c) const {
		assert(0 <= c && c < nChains);
		//return (n + c) / nChains;
		return chainSizes_[c];
	}
	
	//int chainStart(int c) const {
	//	assert(0 <= c && c < nChains);
	//	return n / nChains * c + (n + c) % nChains
	//}
	
	int nBuckets(int c) const {
		//return (chainSize(c) + maxBucketSize - 1) / maxBucketSize;
		return nBuckets_[c];
	}
	
	
	int bucketSize(int c, int b) const {
		assert(0 <= b && b < nBuckets(c));
		//return (chainSize(c) + b) / nBuckets(c);
		return bucketSizes_[c][b];
	}
	
	int indexOf(int c, int b, int i) const {
		int v = 0;
		for (int cc = 0; cc < c; ++cc)
			v += chainSize(cc);
		for (int bb = 0; bb < b; ++bb)
			v += bucketSize(c, bb);
		v += i;
		return v;
	}
	
	void getPosition(int v, int& c, int& b, int& i) const {
		c = 0;
		while (v >= chainSize(c)) {
			v -= chainSize(c);
			++c;
		}
		b = 0;
		while (v >= bucketSize(c, b)) {
			v -= bucketSize(c, b);
			++b;
		}
		i = v;
	}
	
	
	size_t nIdeals() const {
		//size_t ni = 1;
		//for (int c = 0; c < nChains; ++c) {
		//	size_t nic = 1;
		//	for (int b = 0; b < nBuckets(c); ++b)
		//		nic += (1 << bucketSize(c, b)) - 1;
		//	ni *= nic;
		//}
		//return ni;
		return nIdeals_;
	}

	size_t nChainIdeals(int c) const {
		return nChainIdeals_[c];
	}
	
	//int getBucket(int v) const {
	//	return v / maxBucketSize;
	//}
	
	bool operator==(const ParBucketOrderFamily& bof) const {
		return n == bof.n && maxBucketSize == bof.maxBucketSize && nChains == bof.nChains;
	}
	
	class Order {
	private:
		const ParBucketOrderFamily& pof_;
		int* order_;
		Order(const Order&); // disable copy constructor
		//Order& operator=(const Order&); // disable copying
	public:
		Order(const ParBucketOrderFamily& pof) : pof_(pof) {
			order_ = new int[pof_.n];
			for (int i = 0; i < pof_.n; ++i)
				order_[i] = i;
		}
		~Order() {
			delete[] order_;
		}
		
		Order& operator=(const Order& o) {
			assert(pof_ == o.pof_);
			for (int i = 0; i < pof_.n; ++i)
				order_[i] = o.order_[i];
			return *this;
		}
		
//		int& operator[](size_t i) {
//			return order_[i];
//		}
		
		int operator[](size_t i) const {
			return order_[i];
		}
		
//		int* tail(size_t b) {
//			return order_ + b * pof_.maxBucketSize;
//		}
//		
//		int tailLength(size_t b) const {
//			return pof_.n - b * pof_.maxBucketSize;
//		}
		
		const int* getOrder() const {
			return order_;
		}
		
		void print() const {
			int v = 0;
			for (int c = 0; c < pof_.nChains; ++c) {
				for (int b = 0; b < pof_.nBuckets(c); ++b) {
					if (b > 0)
						printf(" ≺  ");
					for (int i = 0; i < pof_.bucketSize(c, b); ++i) {
						printf("%d ", order_[v]);
						++v;
					}
				}
				printf("\n");
			}
		}
		
		void rand() {
			for (int i = 0; i < pof_.n; ++i) {
				int j = i + randuint(pof_.n - i);
				int tmp = order_[i];
				order_[i] = order_[j];
				order_[j] = tmp;
			}
		}
		
		// vai erikseen ketju ja indeksi?
		int getIndex(int j) const {
			for (int i = 0; i < pof_.n; ++i)
				if (order_[i] == j)
					return i;
			printf("j = %d\n", j);
			assert(0);
		}
		
		void randSwap() {
			int c = randuint(pof_.nChains);
			int nb = pof_.nBuckets(c);
			int b1 = randuint(nb);
			int b2 = (b1 + randuint(nb) + 1) % nb;
			int i1 = randuint(pof_.bucketSize(c, b1));
			int i2 = randuint(pof_.bucketSize(c, b2));
			int v1 = pof_.indexOf(c, b1, i1);
			int v2 = pof_.indexOf(c, b2, i2);
			std::swap(order_[v1], order_[v2]);
		}
	};
	
	
//	class OrderEnumerator {
//	private:
//		int** index_;
//		const ParBucketOrderFamily& pof_;
//		Order po_;
//		int b_;		// current bucket
//		int i_;		// current item in current bucket
//	public:
//		OrderEnumerator(const ParBucketOrderFamily& pof) : pof_(pof), po_(pof_) {
//			// allocate index
//			index_ = new int*[pof_.nBuckets() - 1];
//			for (int b = 0; b < pof_.nBuckets() - 1; ++b)
//				index_[b] = new int[pof_.bucketSize(b) + 1];
//			
//			// init first state
//			init();
//		}
//		~OrderEnumerator() {
//			for (int b = 0; b < pof_.nBuckets() - 1; ++b)
//				index_[b] = new int[pof_.bucketSize(b)];
//			delete[] index_;
//		}
//		
//		void init() {
//			b_ = 0;
//			i_ = 0;
//			for (; b_ < pof_.nBuckets() - 1; ++b_) {
//				for (; i_ <= pof_.bucketSize(b_); ++i_)
//					index_[b_][i_] = -1;
//				i_ = 0;
//			}
//			//printf("  b_=%d   i_=%d\n", b_, i_);
//		}
//		
//		bool next() {
//			// find the next item to change
//			while (true) {
//				// no items left in this bucket?
//				if (i_ == 0) {
//					// if no more buckets left => stop
//					if (b_ == 0)
//						return false;
//					// move to the previous bucket
//					--b_;
//					i_ = pof_.bucketSize(b_);
//				}
//				// undo the previous swap, if there was such
//				if (index_[b_][i_] > index_[b_][i_-1]) {
//					std::swap(po_.tail(b_)[i_-1], po_.tail(b_+1)[index_[b_][i_]]);
//					//printf("  b_=%d   i_=%d  index_[b_][i_]=%d\n", b_, i_, index_[b_][i_]);
//					//printf("  "); po_.print(false); printf("\n");
//				}
//				// item with unused values left (possible swaps left)? => go on
//				if (++index_[b_][i_] < po_.tailLength(b_+1))
//					break;
//				// otherwise try previous item
//				--i_;
//			}
//			
//			// new swap
//			std::swap(po_.tail(b_)[i_-1], po_.tail(b_+1)[index_[b_][i_]]);
//			//printf("  b_=%d   i_=%d  index_[b_][i_]=%d\n", b_, i_, index_[b_][i_]);
//			//printf("  "); po_.print(false); printf("\n");
//			
//			// init the tail and move back to the last bucket / item
//			for (; b_ < pof_.nBuckets() - 1; ++b_) {
//				for (++i_; i_ <= pof_.bucketSize(b_); ++i_)
//					index_[b_][i_] = index_[b_][i_ - 1];
//				i_ = 0;
//			}
//			//printf("  b_=%d   i_=%d\n", b_, i_);
//			return true;
//		}
//		
//		const Order& getOrder() const {
//			return po_;
//		}
//	};
	
	
	class Ideal {
	private:
		public: // TODO: poista	publicointi
		const ParBucketOrderFamily& pof_;
		int* hatBuckets_;
		int* hatSetMasks_;
		
		Ideal(const Ideal&); // disable copy constructor
		Ideal& operator=(const Ideal&); // disable copying
	public:
		Ideal(const ParBucketOrderFamily& pof) : pof_(pof) {
			hatBuckets_ = new int[pof_.nChains];
			hatSetMasks_ = new int[pof_.nChains];
			setEmpty();
		}

		~Ideal() {
			delete[] hatBuckets_;
			delete[] hatSetMasks_;
		}
		
		void setEmpty() {
			for (int c = 0; c < pof_.nChains; ++c) {
				hatBuckets_[c] = 0;
				hatSetMasks_[c] = 0;
			}
		}
		
		void setFull() {
			for (int c = 0; c < pof_.nChains; ++c) {
				hatBuckets_[c] = pof_.nBuckets(c) - 1;
				hatSetMasks_[c] = (1 << pof_.bucketSize(c, hatBuckets_[c])) - 1;
			}
		}
		
		//void operator+= (int c, int b, int i) {
		//	int vm = 1 << i;
		//	assert(hatBuckets_[c] == b && !(hatSetMask_[c] & vm));
		//	hatSetMask_[c] += vm;
		//}
		//
		//void operator-= (int c, int b, int i) {
		//	int vm = 1 << i;
		//	assert(hatBuckets_[c] == b && (hatSetMask_[c] & vm));
		//	hatSetMask_[c] -= vm;
		//}
		
		
		
		bool isShrinkableWith(int v) {
			int c, b, i; pof_.getPosition(v, c, b, i);
			return (b == hatBuckets_[c]) && (hatSetMasks_[c] & (1 << i));
		}
		
		void expandWith(int v) {
			int c, b, i; pof_.getPosition(v, c, b, i);
			hatSetMasks_[c] += (1 << i);
		}

		void shrinkWith(int v) {
			int c, b, i; pof_.getPosition(v, c, b, i);
			hatSetMasks_[c] -= (1 << i);
		}
		
		void setSuperOf(const StackSubset& ss) {
			setEmpty();
			for (int v = 0; v < ss.size(); ++v) {
				int c, b, i;
				pof_.getPosition(ss[v], c, b, i);
				if (b > hatBuckets_[c]) {
					hatBuckets_[c] = b;
					hatSetMasks_[c] = (1 << i);
				} else if (b == hatBuckets_[c]) {
					hatSetMasks_[c] |= (1 << i);
				}
			}
		}
		
		bool next() {
			for (int c = 0; c < pof_.nChains; ++c) {
				if (++hatSetMasks_[c] < (1 << pof_.bucketSize(c, hatBuckets_[c])))
					return true;
				if (++hatBuckets_[c] < pof_.nBuckets(c)) {
					hatSetMasks_[c] = 1;
					return true;
				}
				hatBuckets_[c] = 0;
				hatSetMasks_[c] = 0;
			}
			return false;
		}

		bool prev() {
			for (int c = 0; c < pof_.nChains; ++c) {
				if (--hatSetMasks_[c] >= 0)
					return true;
				if (--hatBuckets_[c] >= 0) {
					hatSetMasks_[c] = (1 << pof_.bucketSize(c, hatBuckets_[c])) - 2;
					return true;
				}
				hatBuckets_[c] = pof_.nBuckets(c) - 1;
				hatSetMasks_[c] = (1 << pof_.bucketSize(c, hatBuckets_[c])) - 1;
			}
			return false;
		}

		//int tailSize() const {
		//	return hatBucket_ * pof_.maxBucketSize;
		//}
		//
		//const StackSubset hat() const {
		//	int bs = pof_.bucketSize(hatBucket_);
		//	StackSubset h(bs);
		//	int hsm = hatSetMask_;
		//	int i = tailSize();
		//	while (hsm) {
		//		if (hsm & 1)
		//			h.push(i);
		//		hsm >>= 1;
		//		++i;
		//	}
		//	return h;
		//}
		
		//template<typename T> friend class IdealMap;
	};
	
	
	template<typename T>
	class IdealMap {
	private:
		const ParBucketOrderFamily& pof_;
		T* data_;
		//T** buckets_;
		
		IdealMap();
	public:
		IdealMap(const ParBucketOrderFamily& pof) : pof_(pof) {
			data_ = new T[pof_.nIdeals()];
		}

		IdealMap(const IdealMap& im) : pof_(im.pof_) {
			size_t size = pof_.nIdeals();
			data_ = new T[size];
			memcpy(data_, im.data_, size * sizeof(T));
		}

		~IdealMap() {
			delete[] data_;
		}
		
		T& operator[] (const Ideal& i) {
			int j = 0;
			for (int c = 0; c < pof_.nChains; ++c) {
				int cj = 0;
				for (int b = 0; b < i.hatBuckets_[c]; ++b)
					cj += (1 << pof_.bucketSize(c, b)) - 1;
				cj += i.hatSetMasks_[c];
				j = j * pof_.nChainIdeals_[c] + cj;
			}
			return data_[j];
		}

		//T* operator[] (size_t j) {
		//	assert(0 <= j && j < pof_.nBuckets());
		//	return data_ + j * ((1 << pof_.maxBucketSize) - 1);
		//}
		
		void setAll(T value) {
			for (int i = 0; i < pof_.nIdeals(); ++i)
				data_[i] = value;
		}
		
		T getEmpty() {
			return data_[0];
		}
		
		T getFull() {
			return data_[pof_.nIdeals() - 1];
		}
		
		void fastSparseZetaTransform() {
			// for each chain
			for (int c = 0; c < pof_.nChains; ++c) {
				// for each bucket
				for (int b = 0; b < pof_.nBuckets(c); ++b) {
					int bucketSize = pof_.bucketSize(c, b);
					// for each variable in bucket
					for (int i = 0; i < bucketSize; ++i) {
						// variable mask
						int vm = 1 << i;
						// enumerate all compatible Y̌:s in the bucket
						Ideal y(pof_);
						do {
							if (y.hatBuckets_[c] != b || !(y.hatSetMasks_[c] & vm))
								continue;
							y.hatSetMasks_[c] -= vm;
							T tmp = (*this)[y];
							y.hatSetMasks_[c] += vm;
							(*this)[y] += tmp;
						} while (y.next());
					}
				}
			}
		}

		void fastSparseUpZetaTransform() {
			// for each chain
			for (int c = 0; c < pof_.nChains; ++c) {
				// for each bucket
				for (int b = pof_.nBuckets(c) - 1; b >= 0; --b) {
					int bucketSize = pof_.bucketSize(c, b);
					// for each variable in bucket
					for (int i = 0; i < bucketSize; ++i) {
						// variable mask
						int vm = 1 << i;
						// enumerate all compatible Y̌:s in the bucket
						Ideal y(pof_);
						do {
							if (y.hatBuckets_[c] != b || !(y.hatSetMasks_[c] & vm))
								continue;
							T tmp = (*this)[y];
							y.hatSetMasks_[c] -= vm;
							(*this)[y] += tmp;
							y.hatSetMasks_[c] += vm;
						} while (y.next());
					}
				}
			}
		}
		
		void sparseForwardSum(std::vector<IdealMap<T> >& alpha) {
			data_[0] = 1.0;
			Ideal y(pof_);
			while (y.next()) {
				// for each chain
				for (int c = 0; c < pof_.nChains; ++c) {
					int b = y.hatBuckets_[c];
					int bucketSize = pof_.bucketSize(c, b);
					// for each compatible variable
					for (int i = 0; i < bucketSize; ++i) {
						// variable mask
						int vm = 1 << i;
						if (y.hatSetMasks_[c] & vm) {
							int v = pof_.indexOf(c, b, i);
							y.hatSetMasks_[c] -= vm;
							T tmp = alpha[v][y] * (*this)[y];
							y.hatSetMasks_[c] += vm;
							(*this)[y] += tmp;
						}
					}
				}		
			}
		}

		void sparseBackwardSum(std::vector<IdealMap<T> >& alpha) {
			data_[pof_.nIdeals() - 1] = 1.0;
			Ideal y(pof_);
			y.setFull();
			while (y.prev()) {
				// for each chain
				for (int c = 0; c < pof_.nChains; ++c) {
					int b = y.hatBuckets_[c];
					int bucketSize = pof_.bucketSize(c, b);
					// for each compatible variable
					for (int i = 0; i < bucketSize; ++i) {
						// variable mask
						int vm = 1 << i;
						if (!(y.hatSetMasks_[c] & vm)) {
							int v = pof_.indexOf(c, b, i);
							y.hatSetMasks_[c] += vm;
							T tmp = (*this)[y];
							y.hatSetMasks_[c] -= vm;
							(*this)[y] += alpha[v][y] * tmp;
						}
					}
				}		
			}
		}
	};
};/**/





/*
class ParBucketOrder {
private:
	int n_;
	int nChains_;
	int chainLength_;
	int bucketSize_;
	int* order_;

public:
	ParBucketOrder(int n, int bucketSize, int nChains) {
		n_ = n;
		bucketSize_ = bucketSize;
		nChains_ = nChains;
		int chainCapacity = (n + nChains_ - 1) / nChains_;
		chainLength_ = (chainCapacity + bucketSize_ - 1) / bucketSize_;
		//order_ = new int[nChains_ * chainLength_ * bucketSize];
		order_ = new int[n_];
		for (int i = 0; i < n_; ++i)
			order_[i] = i;
	}
	
	~ParBucketOrder() {
		delete[] order_;
	}
	
	
	class Ideal {
	private:
		const ParBucketOrder& po_
		int* hatBuckets_;
		int* hatSetMasks_;
	public:
		Ideal(const ParBucketOrder& po) {
			po_ = po;
			hatBuckets_ = new int[po_.nChains_];
			hatSetMasks_ = new int[po_.nChains_];
			for (int i = 0; i < po_.nChains_; ++i) {
				hatBuckets_[i] = 0;
				hatSetMasks_[i] = 0;
			}
		}
		~Ideal() {
			delete[] hatBuckets_;
			delete[] hatSetMasks_;
		}
		
		bool next() {
			int i = 0;
			do {
				if (++hatSetMasks_[i] < (1 << po.bucketSize(i, hatBuckets_[i])))
					return true;
				hatSetMasks_[i] = 0;
				if (++hatBuckets_[i] < po.bucketsInChain(i))
					return true;
				hatBuckets_[i] = 0;
				++i;
			} while (i < po.nChains_);
			return false;
		}
	}
	
	template<typename T>
	class IdealMap {
	private:
		T* data_;

	public:
		IdealMap() {
			data_ = new T[1 + nBuckets * ((1 << maxBucketSize) - 1)];
		}
		~IdealMap() {
			delete[] data_;
		}
	
		T* operator[] (size_t j) {
			assert(0 <= j && j < nBuckets);
			return data_ + j * ((1 << maxBucketSize) - 1);
		}
	};
	
};
/**/

#endif
