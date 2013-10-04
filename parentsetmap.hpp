/*
 *  BEANDisco: parentsetmap class
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

#include "common.hpp"
#include "stacksubset.hpp"

#ifndef PARENTSETMAP_HPP
#define PARENTSETMAP_HPP

unsigned long numKSubsets(unsigned long n, unsigned long k) {
	unsigned long sum = 1;
	unsigned long perLevel = 1;
	for (unsigned long i = 1; i <= k; ++i) {
		perLevel *= (n-i+1);
		perLevel /= i;
		sum += perLevel;
	}
	return sum;
}


class SubsetDirectory {
public:
	const int nElements;
	const int maxSize;

private:
	size_t* lattice_;

	void buildRecursive(int itemsLeft, size_t x, bool* inSet, int lastItem, size_t& nextFree) {
		//for (int i = 0; i < nElements; ++i)
		//	printf("%d", inSet[i]);
		//printf(" => %d \n", x);
	
		if (!itemsLeft)
			return;
		for (int i = 0; i < nElements; ++i) {
			if (inSet[i])
				continue;
			size_t y;
			if (i < lastItem) {
				y = neighbour(x, lastItem);
				y = neighbour(y, i);
				y = neighbour(y, lastItem);
				lattice_[x * nElements + i] = y;
				lattice_[y * nElements + i] = x;
			} else {
				y = nextFree;
				++nextFree;
				++lastItem;
				lattice_[x * nElements + i] = y;
				lattice_[y * nElements + i] = x;
				inSet[i] = true;
				buildRecursive(itemsLeft - 1, y, inSet, lastItem, nextFree);
				inSet[i] = false;
			}
		}
	}
	
	size_t neighbour(size_t index, int element) const {
		return lattice_[index * nElements + element];
	}
	
	SubsetDirectory(const SubsetDirectory&); // disable copying
	SubsetDirectory& operator=(const SubsetDirectory&); // disable copying
public:
	SubsetDirectory(int n, int k) : nElements(n), maxSize(k) {
		unsigned long nSubsets = numKSubsets(nElements, maxSize);
		//printf("numKSubsets = %d\n", nss);
		//printf("nSubsets * nElements = %d\n", nSubsets * nElements);
		lattice_ = new size_t[nSubsets * nElements];
		
		bool* inSet = new bool[nElements];
		for (int i = 0; i < nElements; ++i)
			inSet[i] = false;
		size_t nextFree = 1;
		buildRecursive(maxSize, 0, inSet, -1, nextFree);
		delete[] inSet;
		//printf("%d\n", nextFree);
	}
	
	~SubsetDirectory() {
		delete[] lattice_;
	}
	
	size_t getIndex(const StackSubset& ss) const {
		assert(ss.size() <= maxSize);
		size_t index = 0;
		for (int i = 0; i < ss.size(); ++i)
			index = neighbour(index, ss[i]);
		return index;
	}
};




template <class T>
class ParentsetMap {
public:
	const int nNodes;
	const int maxParents;
private:
	//const SubsetDirectory& directory_;
	const SubsetDirectory directory_;
	T** data_;
	
	ParentsetMap(const ParentsetMap&); // disable copying
	ParentsetMap& operator=(const ParentsetMap&); // disable copying

public:
	/*ParentsetMap(const SubsetDirectory& directory) :
			directory_(directory), nNodes(directory.nElements), maxParents(directory.maxSize) {
		int nss = numKSubsets(nNodes, directory_.maxSize);
		data_ = new T*[nNodes];
		for (int i = 0; i < nNodes; ++i) {
			data_[i] = new T[nss];
			//memset(data_[i], 0, nss * sizeof(T));
		}
	}/**/

	ParentsetMap(int n, int k) :
			nNodes(n), maxParents(k), directory_(n, k) {
		size_t nSubsets = numKSubsets(nNodes, maxParents);
		data_ = new T*[nNodes];
		for (int i = 0; i < nNodes; ++i) {
			data_[i] = new T[nSubsets];
			//memset(data_[i], 0, nss * sizeof(T));
		}
	}/**/
	
	~ParentsetMap() {
		for (int i = 0; i < nNodes; ++i)
			delete[] data_[i];
		delete[] data_;
	}
	
	T* operator[] (size_t i) {
		assert(0 <= i && i < nNodes);
		return data_[i];
	}
	
	T& operator() (size_t i, const StackSubset& subset) {
		assert(0 <= i && i < nNodes);
		size_t j = directory_.getIndex(subset);
		return data_[i][j];
	}

	T operator() (size_t i, const StackSubset& subset) const {
		assert(0 <= i && i < nNodes);
		size_t j = directory_.getIndex(subset);
		return data_[i][j];
	}

	T& operator() (size_t i, size_t j) {
		assert(0 <= i && i < nNodes);
		return data_[i][j];
	}
	
	size_t getParentsetIndex(const StackSubset& subset) const {
		return directory_.getIndex(subset);
	}
};


#endif



