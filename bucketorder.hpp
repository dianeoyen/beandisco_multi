/*
 *  BEANDisco: bucket order definition
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

#ifndef BUCKETORDER_HPP
#define BUCKETORDER_HPP

class BucketOrderFamily {
public:
	const int n;
	const int maxBucketSize;
	BucketOrderFamily(int _n, int _bucketSize) : n(_n), maxBucketSize(_bucketSize) {
		assert(1 <= maxBucketSize && maxBucketSize < 32);
		assert(1 <= n);
	}
	
	~BucketOrderFamily() {
	}
	
	int bucketSize(int i) const {
		assert(0 <= i && i < nBuckets());
		return min(maxBucketSize, n - i * maxBucketSize);
	}
	
	int nBuckets() const {
		return (n + maxBucketSize - 1) / maxBucketSize;
	}
	
	int nIdeals() const {
		int nb = nBuckets();
		return (nb - 1) * (1 << maxBucketSize) + (1 << bucketSize(nb - 1)) - nb + 1;
	}
	
	int getBucket(int v) const {
		return v / maxBucketSize;
	}
	
	bool operator==(const BucketOrderFamily& bof) const {
		return n == bof.n && maxBucketSize == bof.maxBucketSize;
	}
	
	class Order {
	private:
		const BucketOrderFamily& pof_;
		int* order_;
		Order(const Order&); // disable copy constructor
		//Order& operator=(const Order&); // disable copying
	public:
		Order(const BucketOrderFamily& pof) : pof_(pof) {
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
		
		//int& operator[](size_t i) {
		//	return order_[i];
		//}

		int operator[](size_t i) const {
			return order_[i];
		}
		
		int* tail(size_t b) {
			return order_ + b * pof_.maxBucketSize;
		}
		
		int tailLength(size_t b) const {
			return pof_.n - b * pof_.maxBucketSize;
		}
		
		const int* getOrder() const {
			return order_;
		}
		
		void print(bool newline = true) const {
			for (int b = 0; b < pof_.nBuckets(); ++b) {
				if (b > 0)
					printf(" ≺  ");
				for (int i = 0; i < pof_.bucketSize(b); ++i) {
					printf("%d ", order_[b * pof_.maxBucketSize + i]);
				}
			}
			if (newline)
				printf("\n");
		}
		
		void rand() {
			for (int i = 0; i < pof_.n; ++i) {
				//int j = i + std::rand() % (pof_.n - i);
				int j = i + randuint(pof_.n - i);
				int tmp = order_[i];
				order_[i] = order_[j];
				order_[j] = tmp;
			}
		}
		
		int getIndex(int j) const {
			for (int i = 0; i < pof_.n; ++i)
				if (order_[i] == j)
					return i;
			printf("j = %d\n", j);
			assert(0);
		}
		
		template <class ProbComputer, typename Real>
		bool mcmcPOStep(ProbComputer& pc, Real& p) {
			int v1 = randuint(pof_.n);
			int b1 = v1 / pof_.maxBucketSize;
			//int bs1 = min(pof_.maxBucketSize, pof_.n - b1 * pof_.maxBucketSize);
			int bs1 = pof_.bucketSize(b1);
			//int v2 = ((b1 + 1) * pof_.maxBucketSize + randuint(pof_.n - bs1)) % pof_.n;
			int v2 = (b1 * pof_.maxBucketSize + bs1 + randuint(pof_.n - bs1)) % pof_.n;
			std::swap(order_[v1], order_[v2]);
			Real pnew = pc.calcProb(pof_, *this);
			//Real pnew = pc.updateMargin(*this, v1, v2);
			if (randu() < to<double>(pnew / p)) {
				p = pnew;
				return true;
			} else {
				std::swap(order_[v1], order_[v2]);
				//pc.updateMargin(*this, v1, v2);
				return false;
			}
		}
		
		void randSwap() {
			int v1 = randuint(pof_.n);
			int b1 = v1 / pof_.maxBucketSize;
			int bs1 = pof_.bucketSize(b1);
			int v2 = (b1 * pof_.maxBucketSize + bs1 + randuint(pof_.n - bs1)) % pof_.n;
			std::swap(order_[v1], order_[v2]);
		}
		
		/*void mcmcProposeSwap(int& v1, int& v2) {
			v1 = randuint(pof_.n);
			int b1 = v1 / pof_.maxBucketSize;
			//int bs1 = min(pof_.maxBucketSize, pof_.n - b1 * pof_.maxBucketSize);
			int bs1 = pof_.bucketSize(b1);
			//int v2 = ((b1 + 1) * pof_.maxBucketSize + randuint(pof_.n - bs1)) % pof_.n;
			v2 = (b1 * pof_.maxBucketSize + bs1 + randuint(pof_.n - bs1)) % pof_.n;
			std::swap(order_[v1], order_[v2]);
		}*/
	};
	
	
	class OrderEnumerator {
	private:
		int** index_;
		const BucketOrderFamily& pof_;
		Order po_;
		int b_;		// current bucket
		int i_;		// current item in current bucket
	public:
		OrderEnumerator(const BucketOrderFamily& pof) : pof_(pof), po_(pof_) {
			// allocate index
			index_ = new int*[pof_.nBuckets() - 1];
			for (int b = 0; b < pof_.nBuckets() - 1; ++b)
				index_[b] = new int[pof_.bucketSize(b) + 1];
			
			// init first state
			init();
		}
		~OrderEnumerator() {
			for (int b = 0; b < pof_.nBuckets() - 1; ++b)
				index_[b] = new int[pof_.bucketSize(b)];
			delete[] index_;
		}
		
		void init() {
			b_ = 0;
			i_ = 0;
			for (; b_ < pof_.nBuckets() - 1; ++b_) {
				for (; i_ <= pof_.bucketSize(b_); ++i_)
					index_[b_][i_] = -1;
				i_ = 0;
			}
			//printf("  b_=%d   i_=%d\n", b_, i_);
		}
		
		bool next() {
			// find the next item to change
			while (true) {
				// no items left in this bucket?
				if (i_ == 0) {
					// if no more buckets left => stop
					if (b_ == 0)
						return false;
					// move to the previous bucket
					--b_;
					i_ = pof_.bucketSize(b_);
				}
				// undo the previous swap, if there was such
				if (index_[b_][i_] > index_[b_][i_-1]) {
					std::swap(po_.tail(b_)[i_-1], po_.tail(b_+1)[index_[b_][i_]]);
					//printf("  b_=%d   i_=%d  index_[b_][i_]=%d\n", b_, i_, index_[b_][i_]);
					//printf("  "); po_.print(false); printf("\n");
				}
				// item with unused values left (possible swaps left)? => go on
				if (++index_[b_][i_] < po_.tailLength(b_+1))
					break;
				// otherwise try previous item
				--i_;
			}
			
			// new swap
			std::swap(po_.tail(b_)[i_-1], po_.tail(b_+1)[index_[b_][i_]]);
			//printf("  b_=%d   i_=%d  index_[b_][i_]=%d\n", b_, i_, index_[b_][i_]);
			//printf("  "); po_.print(false); printf("\n");
			
			// init the tail and move back to the last bucket / item
			for (; b_ < pof_.nBuckets() - 1; ++b_) {
				for (++i_; i_ <= pof_.bucketSize(b_); ++i_)
					index_[b_][i_] = index_[b_][i_ - 1];
				i_ = 0;
			}
			//printf("  b_=%d   i_=%d\n", b_, i_);
			return true;
		}
		
		const Order& getOrder() const {
			return po_;
		}
	};
	
	
	class Ideal {
	private:
		public: // TODO: poista	
		const BucketOrderFamily& pof_;
		int hatBucket_;
		int hatSetMask_;
		
	public:
		Ideal(const BucketOrderFamily& pof) : pof_(pof) {
			setEmpty();
		}
		//Ideal(const BucketOrderFamily& pof, int v) : pof_(pof) {
		//	setFirst(v);
		//}
		~Ideal() {
		}
		
		void setEmpty() {
			hatBucket_ = 0;
			hatSetMask_ = 0;
		}
		
		//void setFirst(int v) {
		//	hatBucket_ = pof_.getBucket(v);
		//	hatSetMask_ = 0;
		//}
		
		void getPos(int v, int& b, int& i) {
			assert(0 <= v && v < pof_.n);
			b = v / pof_.maxBucketSize;
			i = v % pof_.maxBucketSize;
		}
		
		bool isExpandableWith(int v) {
			int b, i; getPos(v, b, i);
			return (b == hatBucket_) && !(hatSetMask_ & (1 << i));
		}

		bool isShrinkableWith(int v) {
			int b, i; getPos(v, b, i);
			return (b == hatBucket_) && (hatSetMask_ & (1 << i));
		}
		
		void expandWith(int v) {
			int b, i; getPos(v, b, i);
			hatSetMask_ += (1 << i);
		}

		void shrinkWith(int v) {
			int b, i; getPos(v, b, i);
			hatSetMask_ -= (1 << i);
		}
		
		void setFirstExpandableWith(int v) {
			int b, i; getPos(v, b, i);
			hatBucket_ = b;
			hatSetMask_ = 0;
		}
		
		bool nextExpandableWith(int v) {
			int b, i; getPos(v, b, i);
			int vm = (1 << i);
			assert(b == hatBucket_);
			while (++hatSetMask_ < (1 << pof_.bucketSize(hatBucket_)))
				if (!(hatSetMask_ & vm))
					return true;
			hatSetMask_ = 0;
			return false;
		}
		
		//void nextShrinkableWith(int v) {
		//	
		//}
		
		void setSuperOf(const StackSubset& ss) {
			setEmpty();
			for (int v = 0; v < ss.size(); ++v) {
				int b = ss[v] / pof_.maxBucketSize;
				int i = ss[v] % pof_.maxBucketSize;
				if (b > hatBucket_) {
					hatBucket_ = b;
					hatSetMask_ = (1 << i);
				} else if (b == hatBucket_) {
					hatSetMask_ |= (1 << i);
				}
			}
		}
		
		bool next() {
			if (++hatSetMask_ < (1 << pof_.bucketSize(hatBucket_)))
				return true;
			hatSetMask_ = 1;
			if (++hatBucket_ < pof_.nBuckets())
				return true;
			setEmpty();
			return false;
		}

		//bool next(int v) {
		//	++hatSetMask_;
		//	int vm = (1 << v);
		//	if (hatSetMask_ & vm)
		//		hatSetMask_ += v;
		//	if (hatSetMask_ < (1 << pof_.bucketSize(hatBucket_)))
		//		return true;
		//	hatSetMask_ = 1;
		//	if (++hatBucket_ < pof_.nBuckets())
		//		return true;
		//	hatSetMask_ = 0;
		//	return false;
		//}
		
		//int tailSize() const {
		//	return hatBucket_ * pof_.maxBucketSize;
		//}
		
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
		const BucketOrderFamily& pof_;
		T* data_;
		//T** buckets_;
		
		IdealMap();
	public:
		IdealMap(const BucketOrderFamily& pof) : pof_(pof) {
			data_ = new T[pof_.nIdeals()];
			//data_ = new T[(pof_.nBuckets() +  * (1 << maxBucketSize) - pof_.nBuckets() + 1];
			//buckets_ = new T*[pof_.nBuckets()];
			//for (int i = 0; i < buckets; ++i)
			//	buckets
		}

		IdealMap(const IdealMap& im) : pof_(im.pof_) {
			size_t size = pof_.nIdeals();
			data_ = new T[size];
			memcpy(data_, im.data_, size * sizeof(T));
		}

		/*void copy(const IdealMap& im) {
			assert(&pof_ == im.pof_);
			memcpy(data_, im.data_, pof_.nIdeals() * sizeof(T));
		}*/
		
		void operator=(const IdealMap& im) {
			assert(pof_ == im.pof_);
			memcpy(data_, im.data_, pof_.nIdeals() * sizeof(T));
		}
		
		void operator+=(const IdealMap& im) {
			assert(pof_ == im.pof_);
			for (int i = 0; i < pof_.nIdeals(); ++i)
				data_[i] += im.data_[i];
		}

		void operator-=(const IdealMap& im) {
			assert(pof_ == im.pof_);
			for (int i = 0; i < pof_.nIdeals(); ++i)
				data_[i] -= im.data_[i];
		}

		~IdealMap() {
			delete[] data_;
			//delete[] buckets_;
		}
		
		T& operator[] (const Ideal& i) {
			return data_[i.hatBucket_ * ((1 << pof_.maxBucketSize) - 1) + i.hatSetMask_];
		}

		T* operator[] (size_t j) {
			assert(0 <= j && j < pof_.nBuckets());
			return data_ + j * ((1 << pof_.maxBucketSize) - 1);
		}
		
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
			// for each bucket
			for (int b = 0; b < pof_.nBuckets(); ++b) {
				int bucketSize = pof_.bucketSize(b);
				// index increment from bucket number
				int bi = b * ((1 << pof_.maxBucketSize) - 1);
				// for each variable in bucket
				for (int i = 0; i < bucketSize; ++i) {
					// variable mask
					int vm = 1 << i;
					// enumerate all Y̌:s in the bucket
					for (int yHatI = 1; yHatI < (1 << bucketSize); ++yHatI)
						if (vm & yHatI)
							data_[bi + yHatI] += data_[bi + yHatI - vm];
				}
			}
		}

		void fastSparseUpZetaTransform() {
			// for each bucket
			for (int b = pof_.nBuckets() - 1; b >= 0; --b) {
				int bucketSize = pof_.bucketSize(b);
				// index increment from bucket number
				int bi = b * ((1 << pof_.maxBucketSize) - 1);
				// for each variable in bucket
				for (int i = 0; i < bucketSize; ++i) {
					// variable mask
					int vm = 1 << i;
					// enumerate all Y̌:s in the bucket
					for (int yHatI = 1; yHatI < (1 << bucketSize); ++yHatI)
						if (vm & yHatI)
							data_[bi + yHatI - vm] += data_[bi + yHatI];
				}
			}
		}
		
		void sparseForwardSum(std::vector<IdealMap<T> >& alpha) {
			data_[0] = 1.0;
			// for each bucket
			for (int b = 0; b < pof_.nBuckets(); ++b) {
				int bucketSize = pof_.bucketSize(b);
				// index increment from bucket number
				int bi = b * ((1 << pof_.maxBucketSize) - 1);
				// enumerate all Y̌:s in the bucket
				for (int yHatI = 1; yHatI < (1 << bucketSize); ++yHatI) {
					data_[bi + yHatI] = 0.0;
					// for each direct subideal
					for (int u = 0; u < bucketSize; ++u) {
						// variable mask
						int vMask = (1 << u);
						if (yHatI & vMask) {
							int sHatI = yHatI - vMask;
							// variable index
							int v = b * pof_.maxBucketSize + u;
							data_[bi + yHatI] += alpha[v].data_[bi + sHatI] * data_[bi + sHatI];
						}
					}
				}
			}
		}

		void sparseBackwardSum(std::vector<IdealMap<T> >& alpha) {
			data_[pof_.nIdeals() - 1] = 1.0;
			// for each bucket
			for (int b = pof_.nBuckets() - 1; b >= 0; --b) {
				int bucketSize = pof_.bucketSize(b);
				// index increment from bucket number
				int bi = b * ((1 << pof_.maxBucketSize) - 1);
				// enumerate all Y̌:s in the bucket
				for (int yHatI = (1 << bucketSize) - 2; yHatI >= 0; --yHatI) {
					data_[bi + yHatI] = 0.0;
					// for each direct subideal
					for (int u = 0; u < bucketSize; ++u) {
						// variable mask
						int vMask = (1 << u);
						if (!(yHatI & vMask)) {
							int tHatI = yHatI + vMask;
							// variable index
							int v = b * pof_.maxBucketSize + u;
							data_[bi + yHatI] += alpha[v].data_[bi + yHatI] * data_[bi + tHatI];
						}
					}
				}
			}
		}
	};
};


#endif


