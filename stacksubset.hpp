/*
 *  BEANDisco: stacksubset class
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


#include <cassert>
#include <iostream>


#ifndef STACKSUBSET_HPP
#define STACKSUBSET_HPP


struct StackSubset {
private:
	size_t size_;
	int* items_;
	
	size_t maxSize_;
	
	StackSubset() {
	}
	
	StackSubset& operator=(const StackSubset&); // disable assignment
public:
	StackSubset(int maxSize) {
		maxSize_ = maxSize;
		items_ = new int[maxSize_];
		size_ = 0;
	}
	
	StackSubset(const StackSubset& ss) {
		maxSize_ = ss.maxSize_;
		items_ = new int[maxSize_];
		size_ = 0;
		copyToEnd(ss);
	}
	
	~StackSubset() {
		delete[] items_;
	}
	
	size_t size() const {
		return size_;
	}
	
	int& operator[] (size_t i) {
		assert(0 <= i && i < size_);
		return items_[i];
	}

	int operator[] (size_t i) const {
		//if (!(0 <= i && i < size_))
		//	print_trace(stderr, __FILE__, __LINE__);
		assert(0 <= i && i < size_);
		return items_[i];
	}
	
	void clear() {
		size_ = 0;
	}
	
	void push(int v) {
		assert(size_ < maxSize_);
		items_[size_++] = v;
	}

	int pop() {
		assert(size_ > 0);
		return items_[--size_];
	}
	
	void copyToEnd(const StackSubset& b) {
		for (int i = 0; i < b.size(); ++i)
			push(b.items_[i]);
	}
	
	void fromMask(int mask, int start = 0) {
		int v = start;
		size_ = 0;
		while (mask) {
			assert(size_ < maxSize_);
			if (mask & 1)
				items_[size_++] = v;
			++v;
			mask >>= 1;
		}
	}
	
	bool next(int valuesFirst, int valuesEnd, int maxSize) {
		assert(size_ <= maxSize);
		//assert(maxSize <= valuesEnd - valuesFirst);
		//assert(valuesFirst < valuesEnd);
		if (valuesEnd <= valuesFirst)
			return false;
		if (maxSize == 0)
			return false;
		if (size_ == 0 || (size_ < maxSize && items_[size_ - 1] > valuesFirst)) {
			push(valuesFirst);
			return true;
		} else {
			while (true) {
				++items_[size_ - 1];
				if (size_ == 1) {
					if (items_[size_ - 1] >= valuesEnd) {
						--size_;
						return false;
					} else
						return true;
				} else if (items_[size_ - 1] < items_[size_ - 2])
					return true;
				--size_;
			}
		}
	}
	
	bool contains(int x) const {
		for (int i = 0; i < size_; ++i)
			if (x == items_[i])
				return true;
		return false;
	}
	
	friend std::ostream& operator<<(std::ostream& os, const StackSubset& ss);
};/**/


std::ostream& operator<<(std::ostream& os, const StackSubset& ss) {
	if (ss.size() == 0)
		os << "∅";
	else {
		os << "{" << ss.items_[0];
		for (int i = 1; i < ss.size(); ++i)
			os << "," << ss.items_[i];
		os << "}";
	}
	//for (int i = s.length; i < width; ++i)
	//	os << "  ";
	return os;
}


void translateSubset(const int* order, const StackSubset& a, StackSubset& b) {
	b.clear();
	for (int i = 0; i < a.size(); ++i)
		b.push(order[a[i]]);
}


#endif

