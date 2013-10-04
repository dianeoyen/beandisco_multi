/*
 *  BEANDisco: main program
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

#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_result.h>

#include "common.hpp"
#include "logger.hpp"
#include "lognum.hpp"
#include "timer.hpp"
#include "data.hpp"
#include "stacksubset.hpp"
#include "scores.hpp"
#include "parentsetmap.hpp"
#include "bucketorder.hpp"
#include "parbucketorder.hpp"

//#define NDEBUG
#include <cassert>


#define BEAND_VERSION_STRING  "1.0.1"


// type definitions
typedef Lognum<double> Real;


// create a logger
Logger logger;


struct Arc {
	int tail;
	int head;
	
	void setFirst() {
		head = 0;
		tail = 1;
	}
	
	bool next(int nNodes) {
		++tail;
		if (head == tail)
			++tail;
		if (tail >= nNodes) {
			tail = 0;
			++head;
			if (head >= nNodes) {
				head = 0;
				return false;
			}
		}
		return true;
	}

	bool holds(int v, const StackSubset& pa) {
		return (head != v || pa.contains(tail));
	}
};

std::ostream& operator<<(std::ostream& os, const Arc& arc) {
	return os << arc.tail << " -> " << arc.head;
}

const Arc NullArc = { -1, -1 };



/**
 * A templated map data structure with Arc as index type.
 */
template <class T>
class ArcMap {
private:
	int nNodes_;
	T* data_;
	
	ArcMap(const ArcMap&); // disable copy constructor
	ArcMap& operator=(const ArcMap&); // disable copying
	
public:
	ArcMap(int nNodes) {
		nNodes_ = nNodes;
		data_ = new T[nNodes * nNodes];
	}
	
	~ArcMap() {
		delete[] data_;
	}
	
	void setAll(T value) {
		for (int i = 0; i < nNodes_ * nNodes_; ++i)
			data_[i] = value;
	}
	
	T& operator[] (Arc arc) {
		return data_[arc.head + arc.tail * nNodes_];
	}

	T operator[] (Arc arc) const {
		return data_[arc.head + arc.tail * nNodes_];
	}
};



/**
 * Computes single-task K2 scores for all node-parentset pairs for all tasks
 */
void computeScores(const Data* data, ParentsetMap<Real>** scores, 
		   int nTasks) {
  //computeLogGammas(2 * data[0].nSamples);
  StackSubset parents(scores[0]->maxParents);
  for (int node = 0; node < scores[0]->nNodes; ++node) {
    parents.clear();
    do {
      if (parents.contains(node))
	continue;
      double *logScores = new double[nTasks];
      computeScore(logScores, data, parents, node, 
		   nTasks);
      for (int t = 0; t < nTasks; t++)
	{
	  Lognum<double> tmp;
	  tmp.setLog(logScores[t]);
	  (*scores[t])(node, parents) = to<Real>(tmp);
	}
    } while (parents.next(0, scores[0]->nNodes, scores[0]->maxParents));
  }
  freeLogGammas();
}

int calcDelta(StackSubset p1, StackSubset p2) {
  int delta = 0;
  for (int i = 0; i < p1.size(); i++) {
    if (!p2.contains(p1[i])) 
      delta++;
  }
  return delta;
}

Real* preCalcPDeltas(int maxN) {
  Real* pDeltas = new Real[maxN];
  if (maxN == 1) {
    pDeltas[0] = 1;
    return pDeltas;
  }
  double logZ = maxN * log(4);
  double loghyp;
  for (int i = 0; i < maxN; i++) {
    loghyp = log(gsl_sf_hyperg_2F1(1.0, maxN - 1.0, i + 2.0, 0.25));
    pDeltas[i].setLog(loghyp - logZ - log(i+1));
  }
  return pDeltas;
}


/**
 * Computes multitask scores for all node-parentset pairs for all tasks with parentset prior
 **/
void computeMultiScores(ParentsetMap<Real>** multiScores,  ParentsetMap<Real>** scores,
			int nTasks, int maxParentSets) {
  Real* pDeltas = preCalcPDeltas(scores[0]->nNodes);
  StackSubset parents(scores[0]->maxParents);
  StackSubset pa(scores[0]->maxParents);
  size_t** topPa = new size_t*[nTasks];
  int numParentSets = 0;

  for (int node = 0; node < scores[0]->nNodes; ++node) {
    // Sort top family scores for each node in each task
    if (maxParentSets > 0) {
      for (int task = 0; task < nTasks; ++task) {
	topPa[task] = new size_t[maxParentSets];
	numParentSets = 0;
	do {
	  if (pa.contains(node))
	    continue;
	  Real scoreTemp = (*scores[task])(node, pa);

	  // topPa full, see if this score will be added
	  if ((numParentSets == maxParentSets) && 
	      ((*scores[task])(node, topPa[task][numParentSets-1])))
	    continue;
	  // otherwise, it will be added, just need to find where
	  if (numParentSets == 0) {
	    topPa[task][0] = scores[task]->getParentsetIndex(pa);
	    numParentSets++;
	    continue;
	  }
	  if (numParentSets < maxParentSets) {
	    numParentSets++;
	  }
	  int i = numParentSets - 1;
	  while ((i>0) && (scoreTemp > (*scores[task])(node, topPa[task][i-1]))) {
	    topPa[task][i] = topPa[task][i-1];
	    i--;
	  }
	  topPa[task][i] = scores[task]->getParentsetIndex(pa);

	}  while(pa.next(0, scores[0]->nNodes, scores[0]->maxParents));
      }
    }

    // calculate scores
    parents.clear();
    do {
      if (parents.contains(node))
	continue;
      // Initialize prior, parent sets to loop through
      Real* prior = new Real[nTasks];
      for (int task = 0; task < nTasks; ++task) {
	prior[task] = 0;
      }

      // approximate with just top scoring parent sets
      if (maxParentSets > 0) {
	for (int i = 0; i < numParentSets; i++) {
	  for (int task = 0; task < nTasks; ++task) {
	    int delta = calcDelta(parents, topPa[task][i]);
	    prior[task] += (*scores[task])(node, topPa[task][i]) * pDeltas[delta];
	  }
	}
      }

      else { // do it old way
	pa.clear();      
	// loop through parent sets (to examine pairs of parent sets)
	do {
	  if (pa.contains(node))
	    continue;
	  int delta = calcDelta(parents, pa);
	  for (int task = 0; task < nTasks; ++task) {
	    prior[task] += (*scores[task])(node, pa) * pDeltas[delta];
	  }
	} while (pa.next(0, scores[0]->nNodes, scores[0]->maxParents));
      }

      // apply prior from all tasks to each task
      for (int task = 0; task < nTasks; ++task) {
	(*multiScores[task])(node, parents) = (*scores[task])(node, parents);
	for (int j = 0; j < nTasks; ++j) {
	  if (j==task)
	    continue;
	  //if (prior[j] == 0)
	  //  prior[j] = 1;
	  (*multiScores[task])(node, parents) *= prior[j];
	}
      }
    } while (parents.next(0, scores[0]->nNodes, scores[0]->maxParents));
  }

  if (maxParentSets > 0) {
    for (int t; t < nTasks; t++)
      delete[] topPa[t];
  }
  delete[] topPa;
  delete[] pDeltas;
}

template <class T>
T binom(int n, int k) {
	return round(exp(lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1)));
}

template <>
Lognum<double> binom<Lognum<double> >(int n, int k) {
	Lognum<double> res;
	res.setLog(lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1));
	return res;
}


Real* preCalcInvBinoms(int n) {
	Real* invBinoms = new Real[n];
	for (int k = 0; k < n; ++k)
		invBinoms[k] = Real(1.0) / binom<Real>(n, k);
	return invBinoms;
}

/**
 * Divides scores by the number of parent sets with same size.
 */
void weightScores(ParentsetMap<Real>& scores) {
	Real* invBinoms = preCalcInvBinoms(scores.nNodes - 1);
	StackSubset pa(scores.maxParents);
	for (int v = 0; v < scores.nNodes; ++v) {
		//pa.length = 0;
		do {
			scores(v, pa) *= invBinoms[pa.size()];
		} while (pa.next(0, scores.nNodes, scores.maxParents));
	}
	delete[] invBinoms;
}



/*
template <class POF>
void translateParentsetMap(
		const POF& pof,
		const typename POF::Order& po,
		ParentsetMap<Real>& values,
		ParentsetMap<Real>& transValues
		) {
	for (int v = 0; v < values.nNodes; ++v) {
		int vt = po.getOrder()[v];
		StackSubset x(values.maxParents);
		StackSubset xt(values.maxParents);
		do {
			translateSubset(po.getOrder(), x, xt);
			transValues(v, xi) = values(vt, xt);
		} while (x.next(0, values.nNodes, values.maxParents));
	}
}/**/


/*template <class POF>
void calcTailSums(
		const POF& pof,
		const typename POF::Order& po,
		ParentsetMap<Real>& scores,
		int node,
		Arc feature,
		typename POF::template IdealMap<Real>& tailSums
		) {
	typename POF::Ideal y(pof);
	do {
		const StackSubset yHat = y.hat();
		Real sum = 0;
		if (yHat.size() <= scores.maxParents) {
			int maxTailSubsetSize = min((size_t)y.tailSize(), scores.maxParents - yHat.size());
			StackSubset tailPart(maxTailSubsetSize);
			StackSubset x(max((size_t)scores.maxParents, yHat.size()));
			StackSubset xTrans(scores.maxParents);
			//tailPart.clear();
			do {
				x.clear();
				x.copyToEnd(yHat);
				x.copyToEnd(tailPart);
				translateSubset(po.getOrder(), x, xTrans);
				if (feature.holds(node, xTrans))
					sum = sum + scores(node, xTrans);
				//if (feature.holds(x))
				//	sum = sum + scores(node, x);
			} while (tailPart.next(0, y.tailSize(), maxTailSubsetSize));
		}
		tailSums[y] = sum;
	} while (y.next());
}/**/


/**
 * Computes tail sums of scores in given partial order.
 */
template <class POF>
void calcTailSums(
		const POF& pof,
		const typename POF::Order& po,
		ParentsetMap<Real>& scores,
		int node,
		Arc feature,
		typename POF::template IdealMap<Real>& tailSums
		) {
	tailSums.setAll(0.0);
	StackSubset x(scores.maxParents);
	StackSubset xt(scores.maxParents); // translated x
	typename POF::Ideal xId(pof);
	do {
		translateSubset(po.getOrder(), x, xt);
		//size_t xti = scores.getParentsetIndex(xt);
		xId.setSuperOf(x);
		if (feature.holds(node, xt))
			tailSums[xId] += scores(node, xt);
	} while (x.next(0, scores.nNodes, scores.maxParents));
}/**/


/*template <class POF>
void updateTailSums(
		const POF& pof,
		const typename POF::Order& po,
		ParentsetMap<Real>& scores,
		int node,
		typename POF::template IdealMap<Real>& tailSums,
		int node1,
		int node2
		) {
	if (node == node1 || node == node2) {
		calcTailSums(pof, po, scores, po[node], NullArc, tailSums);
		return;
	}
	StackSubset x(scores.maxParents);
	StackSubset x1t(scores.maxParents); // translated x
	StackSubset x2t(scores.maxParents); // translated x
	typename POF::Ideal x1Id(pof);
	typename POF::Ideal x2Id(pof);
	typename POF::template IdealMap<Real> additions(pof);
	typename POF::template IdealMap<Real> subtractions(pof);
	additions.setAll(0.0);
	subtractions.setAll(0.0);
	int nodet = po[node];
	do {
		if (x.contains(node1) || x.contains(node2))
			continue;
		x.push(node1);
		translateSubset(po.getOrder(), x, x1t);
		x1Id.setSuperOf(x);
		x.pop();
		x.push(node2);
		translateSubset(po.getOrder(), x, x2t);
		x2Id.setSuperOf(x);
		x.pop();
		//tailSums[x1Id] = tailSums[x1Id] + scores(node, x1t) - scores(node, x2t);
		//tailSums[x2Id] = tailSums[x2Id] + scores(node, x2t) - scores(node, x1t);
		additions[x1Id] += scores(nodet, x1t);
		subtractions[x1Id] += scores(nodet, x2t);
		additions[x2Id] += scores(nodet, x2t);
		subtractions[x2Id] += scores(nodet, x1t);
	} while (x.next(0, scores.nNodes, scores.maxParents - 1));

	{
		printf("      sub: ");
		typename POF::Ideal id(pof);
		do {
			printf("%8.6g ", to<double>(log(subtractions[id])));
		} while (id.next());
		printf("\n");
	}
	tailSums -= subtractions;
	{
		printf("        => ");
		typename POF::Ideal id(pof);
		do {
			printf("%8.6g ", to<double>(log(tailSums[id])));
		} while (id.next());
		printf("\n");
	}
	{
		printf("      add: ");
		typename POF::Ideal id(pof);
		do {
			printf("%8.6g ", to<double>(log(additions[id])));
		} while (id.next());
		printf("\n");
	}
	tailSums += additions;
	//{
	//	printf("        => ");
	//	typename POF::Ideal id(pof);
	//	do {
	//		printf("%8.6g ", to<double>(log(tailSums[id])));
	//	} while (id.next());
	//	printf("\n");
	//}
}/**/



/**
 * Computes alphas from scores.
 */
template <class POF>
void calcAlphas(
		const POF& pof,
		const typename POF::Order& po,
		ParentsetMap<Real>& scores,
		Arc arc,
		std::vector<typename POF::template IdealMap<Real> >& alphas
		) {
	for (int v = 0; v < pof.n; ++v) {
//printf("  v = %d\n", v);
		calcTailSums(pof, po, scores, po[v], arc, alphas[v]);
		//calcTailSums(pof, po, scores, v, arc, alphas[v]);

//{
//	printf("    ");
//	typename POF::Ideal i(pof);
//	do {
//		printf("%g  ", to<double>(log(alphas[v][i])));
//	} while (i.next());
//	printf("\n");
//}
		
		alphas[v].fastSparseZetaTransform();

//{
//	printf("    ");
//	typename POF::Ideal i(pof);
//	do {
//		printf("%g  ", to<double>(log(alphas[v][i])));
//	} while (i.next());
//	printf("\n");
//}
	}
}



/*
template <class POF>
void fastSparseZetaTransform(
		const POF& pof,
		typename POF::template IdealMap<Real>& im
	) {
	// for each variable
	for (int v = 0; v < pof.n; ++v) {
		// iterate over all ideals y
		typename POF::Ideal y(pof);
		//y.setFirstExpandableWith(v);
		//do {
		//	y.shrinkWith(v);
		//	Real tmp = im[y];
		//	y.expandWith(v);
		//	im[y] += tmp;
		//} while(y.nextExpandableWith(v));
		do {
			if (y.isShrinkableWith(v)) {
				y.shrinkWith(v);
				Real tmp = im[y];
				y.expandWith(v);
				im[y] += tmp;
			}
		} while(y.next());
	}
}/**/



/*template <class POF>
Real basicForwardBackwardSum(
		const POF& pof,
		typename POF::template IdealMap<Real>& fp,
		typename POF::template IdealMap<Real>& bp,
		std::vector<typename POF::template IdealMap<Real> >& alpha,
		int v) {
	int b = v / pof.maxBucketSize;
	int i = v % pof.maxBucketSize;
	int bucketSize = pof.bucketSize(b);
	// index increment from bucket number
	int bi = b * ((1 << pof.maxBucketSize) - 1);
	// variable mask
	int vm = 1 << i;
	Real sum = 0.0;
	for (int yHatI = 0; yHatI < (1 << bucketSize); ++yHatI)
		if (!(vm & yHatI))
			//sum += alpha[v][bi + yHatI] * fp[bi + yHatI] * bp[bi + yHatI + vm];
			sum += alpha[v][b][yHatI] * fp[b][yHatI] * bp[b][yHatI + vm];
	return sum;
}/**/

/*template <class POF>
void calcGammas(
		const POF& pof,
		typename POF::template IdealMap<Real>& fp,
		typename POF::template IdealMap<Real>& bp,
		std::vector<typename POF::template IdealMap<Real> >& gamma
		) {
	int v = 0;
	// for each bucket
	for (int b = 0; b < pof.nBuckets(); ++b) {
		int bucketSize = pof.bucketSize(b);
		// index increment from bucket number
		int bi = b * ((1 << pof.maxBucketSize) - 1);
		// for each variable in bucket
		for (int i = 0; i < bucketSize; ++i) {
			// variable mask
			int vm = 1 << i;
			// enumerate all Y̌:s in the bucket
			for (int yHatI = 0; yHatI < (1 << bucketSize) - 1; ++yHatI)
				if (!(vm & yHatI))
					gamma[v][b][yHatI] = fp[b][yHatI] * bp[b][yHatI + vm];
			gamma[v].fastSparseUpZetaTransform();
			++v;
		}
	}
}/**/


/**
 * Computes gammas from forward and backward sums.
 */
template <class POF>
void calcGammas(
		const POF& pof,
		typename POF::template IdealMap<Real>& fp,
		typename POF::template IdealMap<Real>& bp,
		std::vector<typename POF::template IdealMap<Real> >& gamma
		) {
	// for each variable
	for (int v = 0; v < pof.n; ++v) {
		// iterate over all ideals y
		typename POF::Ideal y(pof);
		//y.setFirstExpandableWith(v);
		//do {
		//	y.expandWith(v);
		//	Real tmp = bp[y];
		//	y.shrinkWith(v);
		//	gamma[v][y] = fp[y] * tmp;
		//} while(y.nextExpandableWith(v));
		do {
			if (y.isShrinkableWith(v)) {
				Real tmp = bp[y];
				y.shrinkWith(v);
				gamma[v][y] = fp[y] * tmp;
				y.expandWith(v);
			}
		} while(y.next());
		gamma[v].fastSparseUpZetaTransform();
	}
}/**/


/**
 * Computes the final (unnormalized) probabilities for each arc from gammas and local scores.
 */
template <class POF>
void addParentsetSums(
		const POF& pof,
		ParentsetMap<Real>& scores,
		std::vector<typename POF::template IdealMap<Real> >& gammas,
		const typename POF::Order& po,
		ArcMap<Real>& sums
		) {
	StackSubset pa(scores.maxParents);
	StackSubset pat(scores.maxParents);
	typename POF::Ideal paId(pof);
	do {
		translateSubset(po.getOrder(), pa, pat);
		//size_t pai = scores.getParentsetIndex(pa);
		size_t pati = scores.getParentsetIndex(pat);
		paId.setSuperOf(pa);
		
		Arc arc;
		for (int i = 0; i < pat.size(); ++i) {
			arc.tail = pat[i];
			for (int headt = 0; headt < pof.n; ++headt) {
				arc.head = po.getOrder()[headt];
				if (arc.head == arc.tail)
					continue;
				sums[arc] += scores(arc.head, pati) * gammas[headt][paId];
			}
		}
	} while (pa.next(0, scores.nNodes, scores.maxParents));
	
}/**/


/**
 * Computes the (unnormalized) probability of given arc in given partial order.
 */
template <class POF>
Real calcUnnormProb(
		const POF& pof,
		ParentsetMap<Real>& scores,
		const typename POF::Order& po,
		Arc arc
	) {
//po.print(); printf("\n");
	
	// compute alphas
	std::vector<typename POF::template IdealMap<Real> >
			alphas(pof.n, typename POF::template IdealMap<Real>(pof));
	calcAlphas(pof, po, scores, arc, alphas);
	
	// compute the probability
	typename POF::template IdealMap<Real> fp(pof);
	fp.sparseForwardSum(alphas);

//{		
//	printf("sparse forward sum =>\n  ");
//	typename POF::Ideal i(pof);
//	do {
//		printf("%g  ", to<double>(log(fp[i])));
//	} while (i.next());
//	printf("\n");
//}

	Real p = fp.getFull();
	//typename POF::template IdealMap<Real> fp(pof);
	//bp.sparseBackwardSum(alphas);
	//Real p = bp.getEmpty();
	//Real p = basicForwardBackwardSum(pof, fp, bp, alphas, po.getIndex(arc.tail));
	
	return p;
}


/**
 * Computes the (unnormalized) probabilities of all arc simultaneously in given partial order.
 */
template <class POF>
void calcUnnormArcProbs(
		const POF& pof,
		ParentsetMap<Real>& scores,
		const typename POF::Order& po,
		ArcMap<Real>& probs
	) {
	// compute alphas for null feature
	std::vector<typename POF::template IdealMap<Real> >
			nullAlphas(pof.n, typename POF::template IdealMap<Real>(pof));
	calcAlphas(pof, po, scores, NullArc, nullAlphas);
	
	// compute forward and backward functions
	typename POF::template IdealMap<Real> fp(pof);
	fp.sparseForwardSum(nullAlphas);
	typename POF::template IdealMap<Real> bp(pof);
	bp.sparseBackwardSum(nullAlphas);
	
	// compute gammas
	std::vector<typename POF::template IdealMap<Real> >
			gammas(pof.n, typename POF::template IdealMap<Real>(pof));
	calcGammas(pof, fp, bp, gammas);
	
	// compute all arc probs at once
	probs.setAll(0.0);
	addParentsetSums(pof, scores, gammas, po, probs);
}



/*
class MarginProbComputer {
private:
	ParentsetMap<Real>& scores_;

	MarginProbComputer(const MarginProbComputer&); // disable copying
	MarginProbComputer& operator=(const MarginProbComputer&); // disable copying
public:
	MarginProbComputer(ParentsetMap<Real>& scores) : scores_(scores) {}
	
	template <class POF>
	Real calcProb(const POF& pof, const typename POF::Order& po) {
		return calcUnnormProb(pof, scores_, po, NullArc);
	}
};/**/


/*template <class POF>
class ProbComputer {
private:
	const POF& pof_;
	ParentsetMap<Real>& scores_;
	std::vector<typename POF::template IdealMap<Real> > marginTailSums_;
	std::vector<typename POF::template IdealMap<Real> > marginAlphas_;
	//const typename POF::Order& po;
	
	ProbComputer(const ProbComputer&); // disable copying
	ProbComputer& operator=(const ProbComputer&); // disable copying
public:
	ProbComputer(ParentsetMap<Real>& scores, const POF& pof) :
		scores_(scores),
		marginTailSums_(pof.n, typename POF::template IdealMap<Real>(pof)),
		marginAlphas_(pof.n, typename POF::template IdealMap<Real>(pof)),
		pof_(pof)
		{}
	
	Real initMargin(const typename POF::Order& po) {
		for (int v = 0; v < pof_.n; ++v) {
			calcTailSums(pof_, po, scores_, po[v], NullArc, marginTailSums_[v]);
			//calcTailSums(pof, po, scores, v, arc, alphas[v]);
			marginAlphas_[v] = marginTailSums_[v];
			marginAlphas_[v].fastSparseZetaTransform();
		}
		
		typename POF::template IdealMap<Real> fp(pof_);
		fp.sparseForwardSum(marginAlphas_);
		Real p = fp.getFull();
		return p;
	}
	
	Real updateMargin(const typename POF::Order& po, int node1, int node2) {
		std::vector<typename POF::template IdealMap<Real> > marginTailSums2(marginTailSums_);
		std::vector<typename POF::template IdealMap<Real> > marginAlphas2(marginAlphas_);
		printf("\nSWAP %d and %d\n", node1, node2);
		for (int v = 0; v < pof_.n; ++v) {
			printf(" v = %d\n", v);
			{
				printf("## origin: ");
				typename POF::Ideal id(pof_);
				do {
					printf("%8.6g ", to<double>(log(marginTailSums_[v][id])));
				} while (id.next());
				printf("\n");
			}
			//updateTailSums(pof_, po, scores_, v, marginTailSums_[v], node1, node2);
			updateTailSums(pof_, po, scores_, v, marginTailSums2[v], node1, node2);
			{
				printf("## update: ");
				typename POF::Ideal id(pof_);
				do {
					//printf("%8.6g ", to<double>(log(marginTailSums_[v][id])));
					printf("%8.6g ", to<double>(log(marginTailSums2[v][id])));
				} while (id.next());
				printf("\n");
			}
			calcTailSums(pof_, po, scores_, po[v], NullArc, marginTailSums_[v]);
			{
				printf("## recalc: ");
				typename POF::Ideal id(pof_);
				do {
					printf("%8.6g ", to<double>(log(marginTailSums_[v][id])));
				} while (id.next());
				printf("\n");
			}
			printf("\n");
			marginAlphas_[v] = marginTailSums_[v];
			marginAlphas_[v].fastSparseZetaTransform();

			marginAlphas2[v] = marginTailSums2[v];
			marginAlphas2[v].fastSparseZetaTransform();
		}
		
		{
			typename POF::template IdealMap<Real> fp(pof_);
			fp.sparseForwardSum(marginAlphas2);
			Real p = fp.getFull();
			printf("log(p) = %g\n", to<double>(log(p)));
		}
		
		typename POF::template IdealMap<Real> fp(pof_);
		fp.sparseForwardSum(marginAlphas_);
		Real p = fp.getFull();
		printf("log(p) = %g\n", to<double>(log(p)));
		return p;
		//return initMargin(po);
	}
	
	void calcUnnormArcProbs(
			const typename POF::Order& po,
			ArcMap<Real>& probs
		) {
		
		// compute forward and backward functions
		typename POF::template IdealMap<Real> fp(pof_);
		fp.sparseForwardSum(marginAlphas_);
		typename POF::template IdealMap<Real> bp(pof_);
		bp.sparseBackwardSum(marginAlphas_);

		// compute gammas
		std::vector<typename POF::template IdealMap<Real> >
				gammas(pof_.n, typename POF::template IdealMap<Real>(pof_));
		calcGammas(pof_, fp, bp, gammas);
		
		// compute all arc probs at once
		probs.setAll(0.0);
		addParentsetSums(pof_, scores_, gammas, po, probs);
	}
};/**/



template <class POF>
class ExactArcProbComputer {
private:
	ParentsetMap<Real>& scores_;
	const POF& pof_;
	
public:
	ExactArcProbComputer(ParentsetMap<Real>& scores, const POF& pof) :
		scores_(scores), pof_(pof)
	{
	}
	
	~ExactArcProbComputer() {
	}
	
	double calcProb(Arc arc) {
		//ParentsetMap<Real> transScores;
		
		Real cumMarginalLikelihood = 0;
		Real cumArcLikelihood = 0;
		typename POF::OrderEnumerator poe(pof_);
		do {
			Real lhPO = calcUnnormProb(pof_, scores_, poe.getOrder(), NullArc);
			Real lhFPO = calcUnnormProb(pof_, scores_, poe.getOrder(), arc);
			//translateParentsetMap(pof_, subsetDirectory_, poe.getOrder(), scores_, transScores);
			//Real lhPO = calcUnnormProb(pof_, subsetDirectory_, transScores, poe.getOrder(), NullArc);
			//Arc arct;
			//arct.head = poe.getOrder().getIndex(arc.head);
			//arct.tail = poe.getOrder().getIndex(arc.tail);
			//Real lhFPO = calcUnnormProb(pof_, subsetDirectory_, transScores, poe.getOrder(), arct);
			cumMarginalLikelihood += lhPO;
			cumArcLikelihood += lhFPO;
		} while(poe.next());
		
		return to<double>(cumArcLikelihood / cumMarginalLikelihood);
	}
	
	void printAllProbs(std::ostream& resStream) {
		Arc arc; arc.setFirst();
		do {
			double p = calcProb(arc);
			resStream << arc << "   " << p << std::endl;
		} while (arc.next(pof_.n));
	}

	void printArcProbs(std::ostream& resStream) {
		Real cumMarginalProb = 0;
		ArcMap<Real> cumArcProbs(pof_.n);
		ArcMap<Real> probs(pof_.n);
		cumArcProbs.setAll(0.0);
		typename POF::OrderEnumerator poe(pof_);
		do {
			Real marginalProb = calcUnnormProb(pof_, scores_, poe.getOrder(), NullArc);
			calcUnnormArcProbs(pof_, scores_, poe.getOrder(), probs);
			cumMarginalProb += marginalProb;
			Arc arc; arc.setFirst();
			do {
				cumArcProbs[arc] += probs[arc];
			} while (arc.next(pof_.n));
		} while(poe.next());
		
		Arc arc; arc.setFirst();
		do {
			double p = to<double>(cumArcProbs[arc] / cumMarginalProb);
			resStream << arc << "   " << p << std::endl;
		} while (arc.next(pof_.n));
	}
};



template <class POF>
class MCMCArcProbComputer {
private:
	ParentsetMap<Real>& scores_;
	const POF& pof_;
	
	Real marginUnnormProb_;
	typename POF::Order po_;
	
	//Real* burnInProbs_;
	
	//MarginProbComputer marginProbComputer;
	//ProbComputer<POF> probComputer;

	int nAccepts_;
	int nSteps_;
	
	std::ostream* marginStream_;
	
public:
	MCMCArcProbComputer(ParentsetMap<Real>& scores, const POF& pof) :
		scores_(scores), pof_(pof), po_(pof), marginStream_(NULL)
		//, marginProbComputer(scores)
		//, probComputer(scores, pof)
	{
		//printf("Create starting state (random permutation)...\n");
		logger.println(1, "  Initialize starting state (random permutation)...");
		po_.rand();
		
		logger.println(1, "  Compute initial probability...");
		marginUnnormProb_ = calcUnnormProb(pof_, scores_, po_, NullArc);
		//marginUnnormProb_ = probComputer.initMargin(po_);
		
		//burnInProbs_ = NULL;
		
		resetStats();
	}
	
	~MCMCArcProbComputer() {
	}

	void resetStats() {
		nAccepts_ = 0;
		nSteps_ = 0;
	}
	
	double getAcceptRatio() {
		return nAccepts_ / (double) nSteps_;
	}
	
	void setMarginStream(std::ostream& targetStream) {
		marginStream_ = &targetStream;
		//marginStream_rdbuf(targetStream.rdbuf());
	}
	
	//void setBurnInProbArray(Real* arr) {
	//	burnInProbs_ = arr;
	//}
	
	void mcmcStep(int nSwaps = 1) {
		typename POF::Order poNew(pof_);
		poNew = po_;
		for (int i = 0; i < nSwaps; ++i)
			poNew.randSwap();
		Real pnew = calcUnnormProb(pof_, scores_, poNew, NullArc);
		//Real pnew = pc.updateMargin(*this, v1, v2);
		++nSteps_;
		if (randu() < to<double>(pnew / marginUnnormProb_)) {
			po_ = poNew;
			marginUnnormProb_ = pnew;
			++nAccepts_;
			//return true;
		} else {
			//pc.updateMargin(*this, v1, v2);
			//return false;
		}
		if (marginStream_)
			(*marginStream_) << to<double>(log(marginUnnormProb_)) << std::endl;
	}
	
	
	void temperedIdleRun(int nSteps, int nSwaps = 1) {
		for (int i = 1; i <= nSteps; ++i) {
			int ns = 1 + nSwaps * (nSteps - i) / nSteps;
			//printf("%d swaps\n", ns);
			mcmcStep(ns);
		}
	}
	
	void idleRun(int nSteps) {
		for (int i = 0; i < nSteps; ++i) {
			//po_.mcmcPOStep(scores_, marginUnnormProb_);
			//po_.mcmcPOStep(marginProbComputer, marginUnnormProb_);
			//nAccepts = mcmcStep();
			//po_.mcmcPOStep(probComputer, marginUnnormProb_);
			//++nSteps_;
			mcmcStep();
			//if (burnInProbs_)
			//	burnInProbs_[i] = marginUnnormProb_;
//printf("marginUnnormProb_ = %g\n", to<double>(log(marginUnnormProb_)));
		}
	}
	
	double calcProb(int nSamples, int nStepsPerSample, Arc arc, double* probs = NULL) {
		//ParentsetMap<Real> transScores;
		
		double psum = 0.0;
		for (int i = 0; i < nSamples; ++i) {
			for (int j = 0; j < nStepsPerSample; ++j) {
				//nAccepts += po_.mcmcPOStep(scores_, marginUnnormProb_);
				//nAccepts += po_.mcmcPOStep(marginProbComputer, marginUnnormProb_);
				//nAccepts += mcmcStep();
				//nAccepts += po_.mcmcPOStep(probComputer, marginUnnormProb_);
				//++nSteps_;
				mcmcStep();
			}
			Real arcUnnormProb = calcUnnormProb(pof_, scores_, po_, arc);
			//Arc arct;
			//arct.head = po_.getIndex(arc.head);
			//arct.tail = po_.getIndex(arc.tail);
			//Real arcUnnormProb = calcUnnormProb(pof_, transScores, po_, arct);
			double pi = to<double>(arcUnnormProb / marginUnnormProb_);
			if (probs)
				probs[i] = pi;
			psum += pi;
			//if (isnan(psum)) {
			//	printf("%g  %g\n", to<double>(arcUnnormProb), to<double>(marginUnnormProb_));
			//	exit(1);
			//}
		}
		double p = psum / nSamples;
		
		return p;
	}
	
	void printAllProbs(std::ostream& resStream, int nSamples, int nStepsPerSample,
			bool printSamples = false) {
		double* samples = NULL;
		if (printSamples)
			samples = new double[nSamples];
		Arc arc; arc.setFirst();
		do {
			double p = calcProb(nSamples, nStepsPerSample, arc, samples);
			resStream << arc << "   " << p;
			if (printSamples)
				for (int i = 0; i < nSamples; ++i)
					resStream << "  " << samples[i];
			resStream << std::endl;
		} while (arc.next(pof_.n));
		if (printSamples)
			delete[] samples;
	}
	
	void printArcProbs(std::ostream& resStream, int nSamples, int nStepsPerSample,
			bool printSamples = false) {
		//ParentsetMap<Real> transScores;
		
		ArcMap<double> cumArcProbs(pof_.n);
		ArcMap<Real> arcUnnormProbs(pof_.n);
		
		cumArcProbs.setAll(0.0);
		for (int i = 0; i < nSamples; ++i) {
			for (int j = 0; j < nStepsPerSample; ++j) {
				//nAccepts += po_.mcmcPOStep(scores_, marginUnnormProb_);
				//nAccepts += po_.mcmcPOStep(marginProbComputer, marginUnnormProb_);
				//nAccepts += mcmcStep();
				//nAccepts += po_.mcmcPOStep(probComputer, marginUnnormProb_);
				//++nSteps;
				mcmcStep();
			}
			
			//translateParentsetMap(pof_, subsetDirectory_, po_, scores_, transScores);
			
			calcUnnormArcProbs(pof_, scores_, po_, arcUnnormProbs);
			//probComputer.calcUnnormArcProbs(po_, arcUnnormProbs);
			//Arc arct;
			//arct.head = po_.getIndex(arc.head);
			//arct.tail = po_.getIndex(arc.tail);
			//Real arcUnnormProb = calcUnnormProb(pof_, subsetDirectory_, transScores, po_, arct);
			
			Arc arc; arc.setFirst();
			do {
				double p = to<double>(arcUnnormProbs[arc] / marginUnnormProb_);
				if (printSamples)
					resStream << "  " << p;
				cumArcProbs[arc] += p;
			} while (arc.next(pof_.n));
			
			if (printSamples) {
				resStream << std::endl;
			}
		}
		
		if (!printSamples) {
			Arc arc; arc.setFirst();
			do {
				double p = cumArcProbs[arc] / nSamples;
				resStream << arc << "   " << p << std::endl;
			} while (arc.next(pof_.n));
		}
	}
};



/**
 * Writes local scores to file.
 */
void writeScores(std::ostream& file, const ParentsetMap<Real>& scores) {
	file << scores.nNodes << std::endl;
	file << scores.maxParents << std::endl;
	file.precision(16);
	StackSubset parents(scores.maxParents);
	for (int node = 0; node < scores.nNodes; ++node) {
		parents.clear();
		do {
			if (parents.contains(node))
				continue;
			file << scores(node, parents).getLog() << " ";
		} while (parents.next(0, scores.nNodes, scores.maxParents));
		file << std::endl;
	}
}


/**
 * Reads local scores from file.
 */
ParentsetMap<Real>* readScores(std::istream& file) {
	int nNodes, maxParents;
	file >> nNodes;
	file >> maxParents;
	ParentsetMap<Real>* scores = new ParentsetMap<Real>(nNodes, maxParents);
	StackSubset parents(maxParents);
	for (int node = 0; node < nNodes; ++node) {
		parents.clear();
		do {
			if (parents.contains(node))
				continue;
			double tmp;
			file >> tmp;
			if (file.fail())
				throw Exception("File corrupted; could not read all scores.");
			Lognum<double> tmp2;
			tmp2.setLog(tmp);
			(*scores)(node, parents) = to<Real>(tmp2);
		} while (parents.next(0, nNodes, maxParents));
	}
	file >> std::ws;
	if (!file.eof())
		throw Exception("File corrupted; contains more data than expected.");
	return scores;
}





using namespace std;


/**
 * Compute and print exact probabilities.
 */
template <class POF>
void printExact(std::ostream& resStream, std::ostream& logStream,
		const POF& pof, ParentsetMap<Real>& scores) {
	logger.println(1, "Initialize...");
	ExactArcProbComputer<POF> eapc(scores, pof);
	logger.printfln(1, "Actual computation...");
	//eapc.printAllProbs(resStream);
	eapc.printArcProbs(resStream);
}


/**
 * Compute and print MCMC approximated probabilities.
 */
template <class POF>
void printMCMC(std::ostream& resStream, std::ostream& marginStream, std::ostream& logStream,
		const POF& pof, ParentsetMap<Real>& scores, int nBurnInSteps, int nSamples,
		int nStepsPerSample, int nBurnOutSteps, bool printSamples) {
	logger.println(1, "Initialize...");
	MCMCArcProbComputer<POF> mapc(scores, pof);
	mapc.setMarginStream(marginStream);
	
	// burn-in
	if (nBurnInSteps > 0) {
		Timer burninTimer; burninTimer.start();
		logger.printfln(1, "Burn-in (%d steps)...", nBurnInSteps);
		//mapc.idleRun(nBurnInSteps);
		mapc.temperedIdleRun(nBurnInSteps, 1);
		double burninTime = burninTimer.elapsed();
		logStream << "burnin_time = " << burninTime << endl;
		logger.printfln(1, "  Elapsed %.2f s.", burninTime);
	}
	
	// actual sampling
	if (nSamples > 0) {
		logger.printfln(1, "Actual computation (%d samples x %d steps)...",
				nSamples, nStepsPerSample);
		mapc.resetStats();
		Timer samplingTimer; samplingTimer.start();
		//mapc.printAllProbs(resStream, nSamples, nStepsPerSample, printSamples);
		mapc.printArcProbs(resStream, nSamples, nStepsPerSample, printSamples);
		double samplingTime = samplingTimer.elapsed();
		logStream << "sampling_time = " << samplingTime << endl;
		logger.printfln(1, "  Elapsed %.2f s.", samplingTime);
		double acceptRatio = mapc.getAcceptRatio();
		logger.printfln(1, "  Acceptance ratio was %.3g.", acceptRatio);
		logStream << "sampling_acceptance_ratio = " << acceptRatio << endl;
	}
	
	// burn-out
	if (nBurnOutSteps > 0) {
		Timer burnoutTimer; burnoutTimer.start();
		logger.printfln(1, "Burn-out (%d steps)...", nBurnOutSteps);
		mapc.idleRun(nBurnOutSteps);
		double burnoutTime = burnoutTimer.elapsed();
		logStream << "burnout_time = " << burnoutTime << endl;
		logger.printfln(1, "  Elapsed %.2f s.", burnoutTime);
	}
}



/*
 * Main program.
 */

#include <string>
#include <iomanip>
#include <boost/program_options.hpp>

namespace opts = boost::program_options;

int main(int argc, char** argv) {

	string inFilename;
	string outFilename;
	string scoreFilename;
	string logFilename;
	//string weightFilename;
	
	int nTasks;
	//double** weights;

	int nVariables;
	int maxIndegree;
	int maxParentSets;
	int nDataSamples;

	unsigned int rngSeed;
	
	int verbosity;
	
	opts::options_description desc("Options");
	desc.add_options()
		("help,h",          "produce help message")
		("verbose,v",       opts::value<int>(&verbosity)->default_value(0)->implicit_value(1),
		                    "set verbosity level")
		("quiet,q",         "use quiet mode, does not print anything unnecessary")
		("num-rows,r",      opts::value<int>(&nDataSamples)->default_value(0),
		                    "set number of data rows (samples)")
		("num-tasks,t",      opts::value<int>(&nTasks)->default_value(1),
		                    "set number of tasks")
		("num-variables,n", opts::value<int>(&nVariables)->default_value(0),
		                    "set number of variables")
		("max-indegree,m",  opts::value<int>(&maxIndegree)->default_value(0),
		                    "set maximum indegree")
      	        ("approx,a",        opts::value<int>(&maxParentSets)->default_value(0),
	  	                    "only use top <a> scoring parent sets for transfer approximation")
		("seed",            opts::value<unsigned int>(&rngSeed)->default_value(time(0)),
		                    "set seed for random number generator")
		("input-file,i",    opts::value<string>(&inFilename)->default_value(""),
		                    "set input file for data")
		("score-file",      opts::value<string>(&scoreFilename)->default_value(""),
		                    "set score file for reading in STL scores")
	  //("weight-file",     opts::value<string>(&weightFilename)->default_value(""),
	  //	                    "read task-relatedness weights from file")
		("output-file,o",   opts::value<string>(&outFilename)->default_value("-"),
		                    "set output file for MTL scores")
		("log-file",        opts::value<string>(&logFilename)->default_value(""),
		                    "set log file to write statistics about computation")
		;

	opts::positional_options_description pdesc;
	pdesc.add("input-file", 1);
	pdesc.add("output-file", 1);
	
	opts::variables_map vm;
	
	try {
		opts::store(opts::command_line_parser(argc, argv).options(desc).positional(pdesc).run(), vm);
		opts::notify(vm);
	} catch (opts::error& err) {
		logger.println(-1, "Error: ", err.what());
		logger.println(-1, "Aborting.");
		return 1;
	}
	
	if (vm.count("help")) {
		logger.println(-1, "BEANDiscoMulti - Bayesian Exact and Approximate Network Discovery");
		logger.println(-1, "Version " BEAND_VERSION_STRING);
		logger.println(-1);
		logger.println(-1, "Usage:");
		logger.printfln(-1, "  %s [options] [infile [outfile]]", argv[0]);
		logger.println(-1);
		logger.println(-1, desc);
		return 1;
	}
	
	if (vm.count("quiet"))
		logger.setVerbosity(-1);
	else
		logger.setVerbosity(verbosity);
	
	
	logger.println(2, "Parameters:");
	logger.printfln(2, "  data file = %s", inFilename.c_str());
	logger.printfln(2, "  scoreFilename = %s", scoreFilename.c_str());
	logger.printfln(2, "  outFilename = %s", outFilename.c_str());
	logger.printfln(2, "  number_tasks = %d", nTasks);	
	logger.printfln(2, "  nVariables = %d", nVariables);
	logger.printfln(2, "  number_rows = %d", nDataSamples);
	logger.printfln(2, "  maxIndegree = %d", maxIndegree);
	logger.printfln(2, "  approx = %d", maxParentSets);
	
        // if weight file given, read it in
	/*
	if (!weightFilename.empty()) {
	  ifstream wFile;
	  istream weightFile(0);
	  string row;
	  wFile.open(weightFilename.c_str());
	  if (!wFile) {
	    logger.printfln(-1, "Error: Could not open file '%s' for reading.", weightFilename.c_str());
	    return 1;
	  }
	  weightFile.rdbuf(wFile.rdbuf());
	  for (int t = 0; t < nTasks; t++) {
	    if (weightFile.eof()) {
	      logger.printfln(-1, "Error: Not enough rows (%d while %d expected).", t, nTasks);
	      return 1;
	    }
	    getline(weightFile, row);
	    std::istringstream rowStream(row);
	    for (int j = 0; j < nTasks; j++) {
	      double tmp;
	      rowStream >> tmp;
	      if (rowStream.fail()) {
		logger.printfln(-1, "Error: Could not read value on row %d column %d", (t+1), (j+1));
		return 1;
	      }
	      weights[t][j] = tmp;
	    }
	  }
	} else {
	  // weights are uniform (standard MTL)
	  weights = new double*[nTasks];
	  for (int t = 0; t < nTasks; t++)
	    weights[t] = new double[nTasks];
	  for (int t = 0; t < nTasks; t++)
	    for (int j = 0; j < nTasks; j++)
	      weights[t][j] = 1;
	      }*/
	
	// open log stream for statistics
	ofstream logFile;
	ostream logStream(0);
	if (!logFilename.empty()) {
		if (logFilename == "-") {
			logStream.rdbuf(cout.rdbuf());
		} else {
			logFile.open(logFilename.c_str());
			if (!logFile) {
				logger.printfln(-1, "Error: couldn't open file '%s' for writing.", logFilename.c_str());
				return 1;
			}
			logStream.rdbuf(logFile.rdbuf());
		}
	}
	logStream.setf(ios::fixed);
	logStream.precision(2);
	
	logStream << "data_file = " << inFilename << endl;
	logStream << "score_file = " << scoreFilename << endl;
	logStream << "output_file = " << outFilename << endl;
	//logStream << "weight_file = " << weightFilename << endl;
	logStream << "variables = " << nVariables << endl;
	logStream << "rows = " << nDataSamples << endl;
	logStream << "maximum_indegree = " << maxIndegree << endl;
	logStream << "maximum_parent_sets (approx) = " << maxParentSets << endl;
	logStream << "seed = " << rngSeed << endl;
	
	// initialize rng
	rng.seed(rngSeed);
	
	//nBuckets = (nVariables + maxBucketSize - 1) / maxBucketSize;
	
	// Initialize timers
	Timer timer;
	timer.start();
	Timer scoreTimer; 
	double scoreTime;
	
	// map for local scores
	ParentsetMap<Real>** tscores;
	tscores = new ParentsetMap<Real>*[nTasks];
	ParentsetMap<Real>** scores;
	scores = new ParentsetMap<Real>*[nTasks];
	
	// if data file given, read the data and compute the scores
	if (!inFilename.empty()) {
		if (maxIndegree <= 0) {
			logger.println(-1, "Error: The maximum in-degree not given.");
			logger.println(-1, "Aborting.");
			return 1;
		}
		
		logger.println(1, "Reading data...");
		Data *data = new Data[nTasks];
		int *maxarities = NULL;
		for (int t = 0; t < nTasks; t++)
		{
		  istream inStream(0);
		  ifstream inFile;
		  string taskfile;
		  char taskno[10] = "";
		  if (inFilename == "-") {
		    inStream.rdbuf(cin.rdbuf());
		  } else {
		    snprintf(taskno, 10, ".task%d", t+1);
		    taskfile = inFilename + taskno;
		    inFile.open(taskfile.c_str());
		    if (!inFile) {
		      logger.printfln(-1, "Error: Could not open file '%s' for reading.", inFilename.c_str());
		      return 1;
		    }
		    inStream.rdbuf(inFile.rdbuf());
		  }
		  try {
		    if (nDataSamples > 0 || nVariables > 0)
		      data[t].read(inFile, nVariables, nDataSamples);
		    else
		      data[t].read(inStream);
		  } catch (Exception& e) {
		    logger.printfln(-1, "Error: While reading data file '%s': %s", inFilename.c_str(), e.what());
		    return 1;
		  }
		  if (inFile.is_open())
		    inFile.close();
		  // track max arity for each variable
		  if (maxarities == NULL) {
		    //memset(maxarities, 0, data[t].nVariables*sizeof(int));
		    maxarities = new int[data[t].nVariables];
		    for (int v = 0; v < data[t].nVariables; v++) {
		      maxarities[v] = 0;
		    }
		  }
		  for (int v = 0; v < data[t].nVariables; v++) {
		    if (data[t].arities[v] > maxarities[v]) {
		      maxarities[v] = data[t].arities[v];
		    }
		  }
		}
		nVariables = data[0].nVariables;
		for (int v = 0; v < nVariables; v++) {
		  for (int t = 0; t < nTasks; t++) {
		    if (data[t].arities[v] < maxarities[v]) {
		      data[t].adjustArity(v, maxarities[v]);
		    }
		  }
		}
		
		logger.println(1, "Computing STL scores...");
		scoreTimer.start();
		logger.println(2, "  Allocating a score map...");
		for (int t = 0; t < nTasks; t++)
		  tscores[t] = new ParentsetMap<Real>(nVariables, maxIndegree);
		logger.println(2, "  Computing the actual scores...");
		computeScores(data, tscores, nTasks);
		scoreTime = scoreTimer.elapsed();
		logStream << "score_stl_computation_time = " << scoreTime << endl;
		logger.printfln(1, "  Elapsed %.2f s.", scoreTime);

		// For backwards compatibility, if input and score file given,
		// write output MTL scores to "score file", rather than
		// reading in STL scores from "score file"
		if ((!scoreFilename.empty()) && outFilename=="-") {
		  logger.println(1, "Assuming score_file is for output (backwards compatibility");
		  outFilename = scoreFilename;
		  scoreFilename = "";
		}
	} // end reading data to compute STL scores

	else if (!scoreFilename.empty()) { //read STL scores from file
	  logger.println(1, "Reading STL scores...");
	  for (int t = 0; t < nTasks; t++) {
	    istream inStream(0);
	    ifstream inFile;
	    string taskfile;
	    char taskno[10] = "";
	    snprintf(taskno, 10, ".task%d", t+1);
	    taskfile = scoreFilename + taskno;
	    ifstream scoreFile(taskfile.c_str());
	    if (!scoreFile) {
	      logger.printfln(-1, "Error: Could not open file '%s' for reading.",
			      taskfile.c_str());
	      return 1;
	    }
	    try {
	      tscores[t] = readScores(scoreFile);
	    } catch (Exception& e) {
	      logger.printfln(-1, "Error: While reading score file '%s': %s",
			      taskfile.c_str(), e.what());
	      return 1;
	    }
	    scoreFile.close();
	  } // done reading task score files
	  nVariables = tscores[0]->nNodes;
	  maxIndegree = tscores[0]->maxParents;
	} else { // complain if neither data nor scores was given
	  logger.printfln(-1, "Error: either data or STL score file should be provided.");
	  logger.printfln(-1, "Type '%s --help' to see help message.", argv[0]);
	  return 1;
	}

	logger.println(1, "  Computing multitask scores...");
	scoreTimer.start();
	for (int t = 0; t < nTasks; t++)
	  scores[t] = new ParentsetMap<Real>(nVariables, maxIndegree);
	computeMultiScores(scores, tscores, nTasks, maxParentSets);
	scoreTime = scoreTimer.elapsed();
	logStream << "multi_score_computation_time = " << scoreTime << endl;
	logger.printfln(1, "  Elapsed %.2f s.", scoreTime);
		
	// Write MTL scores to file
	if (!outFilename.empty()) {
	  string taskfile;
	  char taskno[10] = "";
	  logger.println(1, "Writing MTL scores...");
	  if (outFilename == "-") {
	    for (int t = 0; t < nTasks; t++)
	      writeScores(cout, *scores[t]);
	  } else {
	    for (int t = 0; t < nTasks; t++) {
	      snprintf(taskno, 10, ".task%d", t+1);
	      taskfile = outFilename + taskno;
	      ofstream outFile(taskfile.c_str());
	      if (!outFile) {
		logger.printfln(-1, "Error: Could not open file '%s' for writing.",
				outFilename.c_str());
		return 1;
	      }
	      writeScores(outFile, *scores[t]);
	      outFile.close();
	    }
	  }
	} // done writing MTL scores
	
	// print timing
	double elapsedTime = timer.elapsed();	
	logger.printfln(0, "MTL Elapsed %g s total.", elapsedTime);
	logStream << "time = " << elapsedTime << endl;

	// release resources
	for (int t = 0; t < nTasks; t++)
	  delete scores[t];
	delete scores;
	
	if (logFile.is_open())
		logFile.close();
	
	return 0;
}


