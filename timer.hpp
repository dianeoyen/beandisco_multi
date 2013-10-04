/*
 *  BEANDisco: timer class
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

#include <sys/times.h>
//#include <unistd.h>

#ifndef TIMER_HPP
#define TIMER_HPP

class Timer {
private:
	//double startTime_;
	clock_t startTime_;

public:
	void start() {
		//startTime_ = clock() / (double)CLOCKS_PER_SEC;
		struct tms timeValues;
		times(&timeValues);
		startTime_ = timeValues.tms_utime;
	}
	
	double elapsed() {
		//double currTime = clock() / (double)CLOCKS_PER_SEC;
		//return currTime - startTime_;
		struct tms timeValues;
		times(&timeValues);
		return (timeValues.tms_utime - startTime_) / (double)sysconf(_SC_CLK_TCK);
	}
};

#endif
