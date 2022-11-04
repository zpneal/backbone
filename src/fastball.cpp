#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
Rcpp::List fastball_cpp(Rcpp::List inputList, int numSwaps) {

  //get number of rows
  int numRows = inputList.length();

  //convert input list into a 2D std::vector
  std::vector<std::vector<int> > oneLocs (numRows);
  for(int i = 0; i < numRows; i++) {
    oneLocs[i] = Rcpp::as<std::vector<int> > (inputList[i]);
  }

  //conduct row swaps a total of (numSwaps) times
  for (int i = 1; i <= numSwaps; i++) {

    //get two random row numbers
    int r1Index = R::runif(0,1) * numRows;
    int r2Index = r1Index;
    while (r2Index == r1Index) {
      r2Index = R::runif(0,1) * numRows;
    }

    //create references to the two rows being mixed
    std::vector<int> & r1 = oneLocs[r1Index];
    std::vector<int> & r2 = oneLocs[r2Index];
    if (r1.size() == 0 || r2.size() == 0) {
      continue;
    }

    //generate iterators for first pass through rows
    std::vector<int>::iterator first1 = r1.begin();
    std::vector<int>::iterator last1 = r1.end();
    std::vector<int>::iterator first2 = r2.begin();
    std::vector<int>::iterator last2 = r2.end();
    size_t intersectionLength = 0;

    //find the length of the intersection
    while (first1!=last1 && first2!=last2)
    {
      if (*first1<*first2) ++first1;
      else if (*first2<*first1) ++first2;
      else {
        intersectionLength += 1;
        ++first1; ++first2;
      }
    }

    //calculate length of symmetric difference
    size_t r1SymDiffSize = r1.size() - intersectionLength;
    size_t r2SymDiffSize = r2.size() - intersectionLength;
    size_t symDiffSize = r1SymDiffSize + r2SymDiffSize;

    if (symDiffSize == 0) {
      continue;
    }

    //create vector of zeros and ones
    //represents which row elements of the symmetric difference are placed in
    std::vector<int> swapLocations (symDiffSize);
    std::fill(swapLocations.begin(), swapLocations.begin() + r1SymDiffSize, 0);
    std::fill(swapLocations.begin() + r1SymDiffSize, swapLocations.end(), 1);

    //shuffle swapLocations using Fisher-Yates shuffle
    for (size_t i = 0; i < swapLocations.size() - 1; i++) {
      size_t j = i + R::runif(0,1) * (swapLocations.size() - i);
      std::swap(swapLocations[i],swapLocations[j]);
    }

    //create vectors to store output of curveball swaps
    std::vector<std::vector<int> > curveballRows (2);
    curveballRows[0].reserve(r1.size());
    curveballRows[1].reserve(r2.size());

    //generate iterators for sweep through r1 and r2
    first1 = r1.begin();
    last1 = r1.end();
    first2 = r2.begin();
    last2 = r2.end();
    std::vector<int>::iterator swapIterator = swapLocations.begin();

    //compare elements in r1 and r2 until end of a vector is reached
    while (first1!=last1 && first2!=last2)
    {
      //element in row1 is less than row2
      if (*first1<*first2) {
        //use swapLocations to add element to n1 or n2
        curveballRows[*swapIterator].push_back(*first1);
        //increment iterators
        ++swapIterator;
        ++first1;
      }
      //element in row1 is greater than row2
      else if (*first2<*first1) {
        //use swapLocations to add element to n1 or n2
        curveballRows[*swapIterator].push_back(*first2);
        //increment iterators
        ++swapIterator;
        ++first2;
      }
      //element in row1 is equal to row2
      else {
        //add element to poth arrays and increment both iterators
        curveballRows[0].push_back(*first1);
        curveballRows[1].push_back(*first2);
        ++first1; ++first2;
      }
    }

    //pass through remainder of r1
    while (first1 != last1) {
      curveballRows[*swapIterator].push_back(*first1);
      ++swapIterator;
      ++first1;
    }

    //pass through remainder of r2
    while (first2 != last2) {
      curveballRows[*swapIterator].push_back(*first2);
      ++swapIterator;
      ++first2;
    }

    //clear the rows for r1 and r2 in oneLocs
    r1.clear();
    r2.clear();

    //create references to the shuffled rows
    std::vector<int> & newV1 = curveballRows[0];
    std::vector<int> & newV2 = curveballRows[1];

    //insert the data for the shuffled rows back into oneLocs
    r1.insert(r1.end(), newV1.begin(), newV1.end());
    r2.insert(r2.end(), newV2.begin(), newV2.end());
  }

  //Return randomized adjacency list
  Rcpp::List randomizedList (numRows);

  for (int i = 0; i < numRows; i++) {
    Rcpp::IntegerVector temp = Rcpp::wrap(oneLocs[i]);
    randomizedList[i] = temp;
  }
  return randomizedList;
}
