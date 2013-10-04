/*
 *  BEANDisco: general stuff
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

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/format.hpp>

#ifndef COMMON_HPP
#define COMMON_HPP

// min and max utility function
template <typename T>
T min(T x, T y) {
	return x < y ? x : y;
}

template <typename T>
T max(T x, T y) {
	return x > y ? x : y;
}


// random number generators
typedef boost::mt19937 Rng;
Rng rng;

double randu() {
	boost::uniform_01<double> dist;
	return dist(rng);
}
//boost::uniform_01<Rng&> randu(rng);

int randuint(int ceiling) {
	boost::uniform_int<> dist(0, ceiling - 1);
	return dist(rng);
	//boost::variate_generator<Rng&, boost::uniform_int<> > gen(rng, dist);
	//return gen();
}

/*double randu() {
	return rand() / (double) RAND_MAX;
}

int randuint(int ceiling) {
	return rand() % ceiling;
}/**/


using boost::format;
using std::string;

class Exception {
private:
	format msg_;
	
public:
	Exception(const string& msg) : msg_(msg) {
	}
	
	Exception(const format& msg) : msg_(msg) {
	}
	
	template <typename T>
	Exception operator%(const T& x) {
		return Exception(msg_ % x);
	}
	
	const char* what() {
		return str(msg_).c_str();
	}
};



#endif


