#include "beaconfilter.h"

#include <assert.h>
#include <math.h>

#include <iostream>
#include <vector>
#include <time.h>

using beaconfilter::BeaconFilter;

int main(int argc, char **argv) {
  clock_t startTime,endTime;
  size_t total_items = 50000;
  // Create a beacon filter where each item is of type size_t and
  // use 12 bits for each item:
  //    BeaconFilter<size_t, 12> filter(total_items);
  // To enable semi-sorting, define the storage of beacon filter to be
  // PackedTable, accepting keys of size_t type and making 13 bits
  // for each key:
  //   BeaconFilter<size_t, 13, BeaconFilter::PackedTable> filter(total_items);
  double lf = 0.0;
  //double zero = 0.0;
  //double fpr = 0.0;
  for(int i=0; i<1000; ++i)
  {
  BeaconFilter<size_t, 12> filter(total_items);

  // Insert items to this beacon filter
  size_t num_inserted = 0;
  
  for (size_t i = 0; i < total_items; i++, num_inserted++) {
    //if (filter.Add(i) != beaconfilter::Ok || filter.GetLoadFactor() > 0.1) {
    if (filter.Add(i) != beaconfilter::Ok) {
      break;
    }
    //else
    //  std::cout << i << std::endl;
  }

  // Check if previously inserted items are in the filter, expected
  // true for all items
  startTime = clock();
  
  for (size_t i = 0; i < num_inserted; i++) {
    //assert(filter.Contain(i) == beaconfilter::Ok);
    if (filter.Contain(i) != beaconfilter::Ok)
    {
      //sum++;
    }
  }
  endTime = clock();
  //std::cout << "Totle Time : " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << std::endl;

  // Check non-existing items, a few false positives expected
  double total_queries = 0.0;
  size_t false_queries = 0;
  for (size_t i = total_items; i < 2 * total_items; i++) {
    if (filter.Contain(i) == beaconfilter::Ok) {
      false_queries++;
    }
    total_queries++;
  }
  //if(false_queries == 0) zero++;
  // Output the measured false positive rate
  //std::cout << "false positive rate is "
  //          << 100.0 * false_queries / total_queries << "%\n";
  //std::cout << (double)false_queries / total_queries<< std::endl;
  //fpr+=((double)false_queries / total_queries);
  std::cout << "lf:" <<filter.GetLoadFactor() << std::endl;
  lf += filter.GetLoadFactor();
}
  std::cout << lf/1000 << std::endl;
  return 0;
}
