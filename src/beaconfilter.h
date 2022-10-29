#ifndef BEACON_FILTER_BEACON_FILTER_H_
#define BEACON_FILTER_BEACON_FILTER_H_

#include <assert.h>
#include <algorithm>
#include <vector>
#include <bitset>

#include "debug.h"
#include "hashutil.h"
#include "packedtable.h"
#include "printutil.h"
#include "singletable.h"

namespace beaconfilter {
// status returned by a BEACON filter operation
enum Status {
  Ok = 0,
  NotFound = 1,
  NotEnoughSpace = 2,
  NotSupported = 3,
};

// maximum number of BEACON kicks before claiming failure
const size_t kMaxBEACONCount = 500;

// A BEACON filter class exposes a Bloomier filter interface,
// providing methods of Add, Delete, Contain. It takes three
// template parameters:
//   ItemType:  the type of item you want to insert
//   bits_per_item: how many bits each item is hashed into
//   TableType: the storage of table, SingleTable by default, and
// PackedTable to enable semi-sorting
template <typename ItemType, size_t bits_per_item,
          template <size_t> class TableType = SingleTable,
          typename HashFamily = TwoIndependentMultiplyShift>
class BeaconFilter {


  // Number of items stored
  size_t num_items_;

  typedef struct {
    size_t index;
    uint32_t tag;
    uint32_t flag;
    bool used;
  } VictimCache;

  VictimCache victim_;

  HashFamily hasher_;

  inline size_t IndexHash(uint64_t fingerprint) const {
    // table_->num_buckets is always a power of two, so modulo can be replaced
    // with
    // bitwise-and:
    //return hv & (table_->NumBuckets() - 1);
    return fingerprint & (table_->NumBuckets() - 1);  // tag_left
  }

  inline uint32_t TagHash(uint32_t fingerprint) const {
    uint32_t tag_right;
    tag_right = fingerprint & ((1ULL << bits_per_item) - 1);
    tag_right += (tag_right == 0);
    return tag_right;
  }

  inline void GenerateIndexTagHash(const ItemType& item, size_t* index,
                                   uint32_t* tag) const {
    const uint64_t hash = hasher_(item);
    uint64_t tag_h = hash & ((1ULL << (long long unsigned int)(log2(table_->NumBuckets())+bits_per_item))-1);
    *index = IndexHash(tag_h >> bits_per_item);
    *tag = TagHash(tag_h);
  }

  inline size_t AltIndex(const uint32_t index, const uint32_t tag) const {
    // NOTE(binfan): originally we use:
    // index ^ HashUtil::BobHash((const void*) (&tag), 4)) & table_->INDEXMASK;
    // now doing a quick-n-dirty way:
    // 0x5bd1e995 is the hash constant from MurmurHash2
    return IndexHash((uint32_t)(index ^ (tag * 0x5bd1e995)));
  }


  Status AddImpl(const size_t i, const uint32_t tag);

  // load factor is the fraction of occupancy
  double LoadFactor() const { return 1.0 * Size() / table_->SizeInTags(); }

  double BitsPerItem() const { return 8.0 * table_->SizeInBytes() / Size(); }

 public:

  // Storage of items
  TableType<bits_per_item> *table_;
  

  explicit BeaconFilter(const size_t max_num_keys) : num_items_(0), victim_(), hasher_() {
    size_t assoc = 4;
    size_t num_buckets = upperpower2(std::max<uint64_t>(1, max_num_keys / assoc));
    double frac = (double)max_num_keys / num_buckets / assoc;
    std::cout << "frac:" << frac << std::endl;
    if (frac > 0.96) {
      num_buckets <<= 1;
    }
    victim_.used = false;
    table_ = new TableType<bits_per_item>(num_buckets);
  }

  ~BeaconFilter() { delete table_; }


  // Add an item to the filter.
  Status Add(const ItemType &item);

  // Report if the item is inserted, with false positive rate.
  Status Contain(const ItemType &item) const;

  // Delete an key from the filter
  Status Delete(const ItemType &item);

  /* methods for providing stats  */
  // summary infomation
  std::string Info() const;

  // number of current inserted items;
  size_t Size() const { return num_items_; }

  // size of the filter in bytes.
  size_t SizeInBytes() const { return table_->SizeInBytes(); }

  double GetLoadFactor() const { return LoadFactor(); }
};

template <typename ItemType, size_t bits_per_item,
          template <size_t> class TableType, typename HashFamily>
Status BeaconFilter<ItemType, bits_per_item, TableType, HashFamily>::Add(
    const ItemType &item) {
  size_t i;
  uint32_t tag;

  if (victim_.used) {
    return NotEnoughSpace;
  }

  GenerateIndexTagHash(item, &i, &tag);
  return AddImpl(i, tag);
}

template <typename ItemType, size_t bits_per_item,
          template <size_t> class TableType, typename HashFamily>
Status BeaconFilter<ItemType, bits_per_item, TableType, HashFamily>::AddImpl(
    const size_t i, const uint32_t tag) {
  size_t curindex = i;
  uint32_t curtag = tag;
  uint32_t oldtag;
  uint32_t curflag = 0;
  uint32_t oldflag = 0;

  for (uint32_t count = 0; count < kMaxBEACONCount; count++) {
    bool kickout = count > 0;
    oldtag = 0;
    oldflag = 0;
    if (table_->InsertTagToBucket(curindex, curtag, curflag, kickout, oldtag, oldflag, table_->mat)) {
      num_items_++;
      //std::cout << num_items_ << std::endl;
      return Ok;
    }
    if (kickout) {
      curtag = oldtag;
    }
    curflag = (uint32_t)1 - oldflag;
    curindex = AltIndex(curindex, curtag);
  }

  victim_.index = curindex;
  victim_.tag = curtag;
  victim_.flag = curflag;
  victim_.used = true;
  return Ok;
}

template <typename ItemType, size_t bits_per_item,
          template <size_t> class TableType, typename HashFamily>
Status BeaconFilter<ItemType, bits_per_item, TableType, HashFamily>::Contain(
    const ItemType &key) const {
  bool found = false;
  size_t i1, i2;
  uint32_t tag;

  GenerateIndexTagHash(key, &i1, &tag);
  i2 = AltIndex(i1, tag);

  assert(i1 == AltIndex(i2, tag));

  found = victim_.used && (tag == victim_.tag) && (0 == victim_.flag) &&
          (i1 == victim_.index || i2 == victim_.index);

  if (found || table_->FindTagInBuckets(i1, i2, tag)) {
    return Ok;
  } else {
    return NotFound;
  }
}

template <typename ItemType, size_t bits_per_item,
          template <size_t> class TableType, typename HashFamily>
Status BeaconFilter<ItemType, bits_per_item, TableType, HashFamily>::Delete(
    const ItemType &key) {
  size_t i1, i2;
  uint32_t tag;

  GenerateIndexTagHash(key, &i1, &tag);
  i2 = AltIndex(i1, tag);

  if (table_->DeleteTagFromBucket(i1, tag)) {
    num_items_--;
    goto TryEliminateVictim;
  } else if (table_->DeleteTagFromBucket(i2, tag)) {
    num_items_--;
    goto TryEliminateVictim;
  } else if (victim_.used && tag == victim_.tag &&
             (i1 == victim_.index || i2 == victim_.index)) {
    // num_items_--;
    victim_.used = false;
    return Ok;
  } else {
    return NotFound;
  }
TryEliminateVictim:
  if (victim_.used) {
    victim_.used = false;
    size_t i = victim_.index;
    uint16_t tag = victim_.tag;
    AddImpl(i, tag);
  }
  return Ok;
}

template <typename ItemType, size_t bits_per_item,
          template <size_t> class TableType, typename HashFamily>
std::string BeaconFilter<ItemType, bits_per_item, TableType, HashFamily>::Info() const {
  std::stringstream ss;
  ss << "BeaconFilter Status:\n"
     << "\t\t" << table_->Info() << "\n"
     << "\t\tKeys stored: " << Size() << "\n"
     << "\t\tLoad factor: " << LoadFactor() << "\n"
     << "\t\tHashtable size: " << (table_->SizeInBytes() >> 10) << " KB\n";
  if (Size() > 0) {
    ss << "\t\tbit/key:   " << BitsPerItem() << "\n";
  } else {
    ss << "\t\tbit/key:   N/A\n";
  }
  return ss.str();
}
}  // namespace BeaconFilter
#endif  // BEACON_FILTER_BEACON_FILTER_H_
