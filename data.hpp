/*
 *  BEANDisco: data class for reading and writing data
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

#include <cstddef>
#include <exception>
#include <boost/format.hpp>
#include <iostream>
#include <sstream>
#include <string>

#include "common.hpp"

#ifndef DATA_HPP
#define DATA_HPP

typedef unsigned char Datum;


/*class DataReadException : public Exception {
public:
	DataReadException(const string& msg) : Exception(msg) {
	}
};*/


struct Data {
private:
	Data(const Data&); // disable copying
	Data& operator=(const Data&); // disable copying
	
	void computeArities() {
		assert(arities == NULL);
		arities = (int*)malloc(sizeof(int) * nVariables);
		for (int v = 0; v < nVariables; ++v) {
			int arity = 0;
			for (int i = 0; i < nSamples; ++i) {
				if ((*this)(v,i) >= arity)
					arity = (*this)(v,i) + 1;
			}
			arities[v] = arity;
		}
	}
public:
	int nVariables;
	int nSamples;
	Datum* data;
	int* arities;
	
	Data() {
		nVariables = 0;
		nSamples = 0;
		data = NULL;
		arities = NULL;
	}
	
	void clear() {
		nVariables = 0;
		nSamples = 0;
		if (data)
			free(data);
		data = NULL;
		if (arities)
			free(arities);
		arities = NULL;
	}
	
	~Data() {
		clear();
	}
	
	Datum& operator()(int v, int i) {
		//return data[v * nSamples + i];
		return data[i * nVariables + v];
	}
	
	Datum operator()(int v, int i) const {
		//return data[v * nSamples + i];
		return data[i * nVariables + v];
	}
	
        // May need to adjust arity if another task found higher arity
        void adjustArity(int node, int arity) {
	  assert(arities[node] < arity);
	  arities[node] = arity;
	}

        // Modified so that nVars OR nSamps can be specified
	void read(std::istream& file, int nVars, int nSamps) {
	  if (nVars == 0)
	    nVariables = 1;
	  else
	    nVariables = nVars;
	  if (nSamps == 0)
	    nSamples = 1;
	  else
	    nSamples = nSamps;
	  data = (Datum*)malloc(sizeof(Datum) * nVariables * nSamples);

	  int i = 0; // row counter
	  std::string row;

	  // read first row (and get the number of variables, if not set already)
	  if (nVars == 0) {
	    getline(file, row);
	    std::istringstream rowStream(row);
	    int v = 0;
	    rowStream >> std::ws;
	    while(!rowStream.eof()) {
	      if (v >= nVariables) {
		nVariables *= 2;
		data = (Datum*)realloc(data, sizeof(Datum) * nVariables * nSamples);
	      }
	      int tmp;
	      rowStream >> tmp;
	      if (rowStream.fail())
		throw Exception("Could not read value on row 1 column %d.") % (v+1);
	      data[v] = (Datum) tmp;
	      rowStream >> std::ws;
	      ++v;
	    }
	    file >> std::ws;
	    nVariables = v;
	    data = (Datum*)realloc(data, sizeof(Datum) * nVariables * nSamples);
	    i++;
	  }

	  // load the data (starting value of i depends on whether 1st already read)
	  if (nSamps == 0) { // allocate dynamically if nSamples unknown
	    while(!file.eof()) {
	      if (i >= nSamples) {
		nSamples *= 2;
		data = (Datum*)realloc(data, sizeof(Datum) * nVariables *nSamples);
	      }
	      getline(file, row);
	      std::istringstream rowStream(row);
	      for (int v = 0; v < nVariables; ++v) {
		int tmp;
		rowStream >> tmp;
		if (rowStream.fail())
		  throw Exception("Could not read %dth value on row %d") % (v+1) % (i+1);
		data[i * nVariables + v] = (Datum) tmp;
	      }
	      file >> std::ws;
	      ++i;
	    }
	    nSamples = i;
	    data = (Datum*)realloc(data, sizeof(Datum) * nVariables * nSamples);
	  
	  // nSamples known
	  } else { 
	    for (; i < nSamples; ++i) {
	      if (file.eof())
		throw Exception("Not enough rows (%d while %d expected).") % i % nSamples;
	      getline(file, row);
	      std::istringstream rowStream(row);
	      for (int v = 0; v < nVariables; ++v) {
		int tmp;
		//file >> tmp;
		rowStream >> tmp;
		if (rowStream.fail())
		  throw Exception("Could not read value on row %d column %d") % (i+1) % (v+1);
		data[i * nVariables + v] = (Datum) tmp;
	      }
	    }
	  }	
	  computeArities();
	}

	void read(std::istream& file) {
		nVariables = 1;
		nSamples = 1;
		data = (Datum*)malloc(sizeof(Datum) * nVariables * nSamples);
		
		std::string row;
		
		// read the first row (and get the number of variables)
		getline(file, row);
		std::istringstream rowStream(row);
		int v = 0;
		rowStream >> std::ws;
		while(!rowStream.eof()) {
			if (v >= nVariables) {
				nVariables *= 2;
				data = (Datum*)realloc(data, sizeof(Datum) * nVariables * nSamples);
			}
			int tmp;
			rowStream >> tmp;
			if (rowStream.fail())
				throw Exception("Could not read value on row 1 column %d.") % (v+1);
			data[v] = (Datum) tmp;
			rowStream >> std::ws;
			++v;
		}
		file >> std::ws;
		nVariables = v;
		data = (Datum*)realloc(data, sizeof(Datum) * nVariables * nSamples);
		
		// load the data rest of the data
		int i = 1;
		while(!file.eof()) {
			if (i >= nSamples) {
				nSamples *= 2;
				data = (Datum*)realloc(data, sizeof(Datum) * nVariables * nSamples);
			}
			getline(file, row);
			std::istringstream rowStream(row);
			for (int v = 0; v < nVariables; ++v) {
				int tmp;
				rowStream >> tmp;
				if (rowStream.fail())
					throw Exception("Could not read %dth value on row %d") % (v+1) % (i+1);
				data[i * nVariables + v] = (Datum) tmp;
			}
			file >> std::ws;
			++i;
		}
		nSamples = i;
		data = (Datum*)realloc(data, sizeof(Datum) * nVariables * nSamples);
		
		computeArities();
	}

};

#endif

