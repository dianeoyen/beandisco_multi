/*
 *  BEANDisco: class for handling number as logarithms
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

//#include <cmath>
//#include <cfloat>
#include <limits>

//#define NDEBUG
#include <cassert>
//#undef NDEBUG


#ifndef LOGNUM_HPP
#define LOGNUM_HPP


//template<typename T>
//struct __identity { typedef T type; };


template <class T>
struct Lognum {
private:
	T logx_;
	
public:
	
	Lognum() {
		logx_ = -std::numeric_limits<T>::infinity();
	}
	
	template <class U>
	Lognum(U x) {
		logx_ = log(x);
	}

	Lognum operator=(Lognum x) {
		logx_ = x.logx_;
		return *this;
	}
	
	template <class U>
	Lognum operator=(U x) {
		logx_ = log(x);
		return *this;
	}
	
	void setLog(T e) {
		logx_ = e;
	}
	
	T getLog() {
		return logx_;
	}
	
	T getVal() {
		return exp(logx_);
	}
	
	template <typename F>
	operator F() {
		return (F) getVal();
	}
	
	//template <typename FL>
	//operator Lognum<FL>() {
	//	Lognum<FL> tmp;
	//	tmp.setLog(getLog());
	//	return tmp;
	//}
	
	Lognum operator+(Lognum b) {
		Lognum res;
		if (logx_ == -std::numeric_limits<T>::infinity())
			res = b;
		else if (b.logx_ == -std::numeric_limits<T>::infinity())
			res = *this;
		else
			res.logx_ = (logx_ > b.logx_) ?
					logx_ + log1p(exp(b.logx_ - logx_)) :
					b.logx_ + log1p(exp(logx_ - b.logx_));
					//logx_ + log(1.0 + exp(b.logx_ - logx_)) :
					//b.logx_ + log(1.0 + exp(logx_ - b.logx_));
		return res;
	}
	
	Lognum operator+=(Lognum b) {
		//logx_ += log(1 + exp(b.logx_ - logx_));
		*this = *this + b;
		return *this;
	}
	
	Lognum operator-(Lognum b) {
		//assert(logx_ >= b.logx_);
		Lognum res;
		if (logx_ == -std::numeric_limits<T>::infinity())
			res.logx_ = -std::numeric_limits<T>::infinity();
		else if (logx_ < b.logx_)
			res.logx_ = -std::numeric_limits<T>::infinity();
		else
			res.logx_ = logx_ + log1p(-exp(b.logx_ - logx_));
			//res.logx_ = logx_ + log(1.0 - exp(b.logx_ - logx_));
		return res;
	}
	
	Lognum operator-=(Lognum b) {
		*this = *this - b;
		return *this;
	}
	
	Lognum operator*(Lognum b) {
		Lognum res;
		res.logx_ = logx_ + b.logx_;
		return res;
	}
	
	Lognum operator*=(Lognum b) {
		logx_ += b.logx_;
		return *this;
	}
	
	Lognum operator/(Lognum b) {
		Lognum res;
		res.logx_ = logx_ - b.logx_;
		return res;
	}
	
	Lognum operator/=(Lognum b) {
		logx_ -= b.logx;
		return *this;
	}


	bool operator<(Lognum b) {
		return logx_ < b.logx_;
	}

	bool operator>(Lognum b) {
		return logx_ > b.logx_;
	}

	bool operator<=(Lognum b) {
		return logx_ <= b.logx_;
	}

	bool operator>=(Lognum b) {
		return logx_ >= b.logx_;
	}
};


template <typename T>
T log(Lognum<T> x){
	return x.getLog();
}


template <typename T, typename F>
T to(F x){
	return T(x);
}

//template <typename T, typename FL>
//T to(Lognum<FL> x){
//	return to<T>(x.getVal());
//}


//template <typename T>
//T to(Lognum<T> x){
//	return x.getVal();
//}

//template <typename T, typename FL>
//T to(Lognum<FL> x){
//	return T(x.getVal());
//}

//template <typename TL, typename FL>
//Lognum<TL> to<Lognum<TL> >(Lognum<FL> x){
//	return Lognum<TL>(x);
//}



//template <class U>
//U to(double x){
//	return U(x);
//}
//
//template <class U>
//U to(Lognum<double> x){
//	return U(x);
//}
//
//template <>
//double to<double>(Lognum<double> x){
//	return x.getVal();
//}

#endif


