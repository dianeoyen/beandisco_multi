/*
 *  BEANDisco: functions for computing scores
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


#include <cmath>

#include "data.hpp"
#include "stacksubset.hpp"

#ifndef SCORES_HPP
#define SCORES_HPP

double* logGammas = NULL;
void computeLogGammas(int n) {
	logGammas = new double[n + 1];
	for (int i = 0; i <= n; ++i)
		logGammas[i] = lgamma(i);
}

void freeLogGammas() {
	assert(logGammas != NULL);
	delete[] logGammas;
}

double BDELogScore(int nValues, int nParentValues, double* counts, double pseudocount) {
	double score = 0;
	for (int pv = 0; pv < nParentValues; ++pv) {
		double cumCount = 0;
		double cumPCount = 0;
		for (int v = 0; v < nValues; ++v) {
			double c = counts[pv * nValues + v];
			double pc = pseudocount;
			//score += lgamma(c + pseudocount) - lgamma(pseudocount);
			score += lgamma(c + pc) - lgamma(pc);
			cumCount += c;
			cumPCount += pc;
		}
		//score += lgamma(nValues * pseudocount) - lgamma(cumCount + nValues * pseudocount);
		score += lgamma(cumPCount) - lgamma(cumCount + cumPCount);
	}
	return score;
}/**/

/*
  Computes the score for the given node and parent set for all tasks. Returns an
  array of nTasks scores.
  Just do standard single-task counts, will need them later.
 */
void computeScore(double* scores, const Data* data, const StackSubset& parents, 
		  int node, int nTasks) {
        // need vector of Data for all tasks
	int nParentValues = 1;
	for (int i = 0; i < parents.size(); ++i)
		nParentValues *= data[0].arities[parents[i]];
	int nNodeValues = data[0].arities[node];
	// need vector of counts, i.e. 2d array and weighted counts
	double** counts = new double*[nTasks];
	int nValues = nParentValues * nNodeValues;
	for (int t = 0; t < nTasks; t++) {
	  counts[t] = new double[nValues];
	  for (int i = 0; i < nValues; ++i) {
		counts[t][i] = 0;
	  }
	}
	//memset(counts, 0, nParentValues * nNodeValues * sizeof(int));
	
	for (int t = 0; t < nTasks; t++) {
	  for (int j = 0; j < data[t].nSamples; ++j) {
		int index = 0;
		for (int i = 0; i < parents.size(); ++i)
			index = index * data[t].arities[parents[i]] + data[t](parents[i], j);
		index = index * data[t].arities[node] + data[t](node, j);
		++counts[t][index];
	  }
	  scores[t] = BDELogScore(nNodeValues, nParentValues, counts[t], 1.0);
	}
	//for (int t = 0; t < nTasks; t++)
	//{
	//  for (int i = 0; i < nParentValues * nNodeValues; ++i)
	//    printf("%d ", counts[t][i]);
	//  printf("\n");
	//}
	    
	// cleanup
	for (int t = 0; t < nTasks; t++) {
	  delete[] counts[t];
	}
	delete[] counts;
}/**/

#endif

