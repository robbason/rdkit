//
//  Copyright (C) 2014 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#pragma once
#include <list>
#include <algorithm>
#include <cmath>
#include "SubstructMatchCustom.h"

namespace RDKit {

namespace FMCS {
class RDKIT_FMCS_EXPORT RingMatchTableSet {
  class RingMatchTable {
    FMCS::MatchTable MatchMatrix;
    std::map<const INT_VECT*, unsigned> RingIndex;

   public:
    inline void clear() {
      MatchMatrix.clear();
      RingIndex.clear();
    }
    inline void resize(unsigned s1, unsigned s2) {
      MatchMatrix.resize(s1, s2);
      for (size_t i = 0; i < s1; i++)
        for (size_t j = 0; j < s2; j++) MatchMatrix.set(i, j, false);
    }
    inline void makeRingIndex(const ROMol* mol2) {
      unsigned i = 0;
      // for each TARGET ring
      const RingInfo::VECT_INT_VECT& rings2 = mol2->getRingInfo()->bondRings();
      for (RingInfo::VECT_INT_VECT::const_iterator r2 = rings2.begin();
           r2 != rings2.end(); r2++)
        RingIndex[&*r2] = i++;
    }
    inline bool isEqual(unsigned i, const INT_VECT* r2) const {
      return MatchMatrix.at(i, getRingIndex(r2));
    }
    inline void setMatch(unsigned i, const INT_VECT* r2) {
      MatchMatrix.set(i, getRingIndex(r2), true);
    }

   private:
    inline unsigned getRingIndex(const INT_VECT* r2) const {
      std::map<const INT_VECT*, unsigned>::const_iterator j =
          RingIndex.find(r2);
      if (RingIndex.end() == j) throw - 1;
      return j->second;
    }
  };

 private:
  std::vector<std::vector<size_t>>* QueryBondRingsIndeces{nullptr};
  std::map<const ROMol*, std::vector<std::vector<size_t>>>
      TargetBondRingsIndecesSet;  // by target molecules

  std::map<const ROMol*, RingMatchTable> MatchMatrixSet;  // by target molecules
  std::map<const INT_VECT*, unsigned> QueryRingIndex;

 public:
  RingMatchTableSet()  {}

  inline void clear() {
    if (QueryBondRingsIndeces) QueryBondRingsIndeces->clear();
    TargetBondRingsIndecesSet.clear();
    MatchMatrixSet.clear();
    QueryRingIndex.clear();
  }

  inline bool isQueryBondInRing(unsigned bi) const {
    return (*QueryBondRingsIndeces)[bi].empty();
  }
  inline const std::vector<size_t>& getQueryBondRings(unsigned bi) const {
    return (*QueryBondRingsIndeces)[bi];
  }

  inline bool isTargetBondInRing(const ROMol* target, unsigned bi) const {
    std::map<const ROMol*, std::vector<std::vector<size_t>>>::const_iterator i =
        TargetBondRingsIndecesSet.find(target);
    if (TargetBondRingsIndecesSet.end() == i) throw - 1;  // never
    return i->second[bi].empty();
  }
  inline const std::vector<size_t>& getTargetBondRings(const ROMol* target,
                                                       unsigned bi) const {
    std::map<const ROMol*, std::vector<std::vector<size_t>>>::const_iterator i =
        TargetBondRingsIndecesSet.find(target);
    if (TargetBondRingsIndecesSet.end() == i) throw - 1;  // never
    return i->second[bi];
  }

  inline bool isEqual(const INT_VECT* r1, const INT_VECT* r2,
                      const ROMol* mol2) const {
    const RingMatchTable& m = getTargetMatchMatrix(mol2);
    unsigned i = getQueryRingIndex(r1);
    return m.isEqual(i, r2);
  }

  void init(const ROMol* query);

  inline void addTargetBondRingsIndeces(const ROMol* mol2) {
    std::vector<std::vector<size_t>>& m = TargetBondRingsIndecesSet[mol2];
    m.resize(mol2->getNumBonds());

    size_t ri = 0;
    const RingInfo::VECT_INT_VECT& rings = mol2->getRingInfo()->bondRings();
    for (RingInfo::VECT_INT_VECT::const_iterator r = rings.begin();
         r != rings.end(); r++, ri++)
      for (INT_VECT::const_iterator bi = r->begin(); bi != r->end();
           bi++)  // all bonds in the ring
        m[*bi].push_back(ri);
  }

  void computeRingMatchTable(
      const ROMol* query, const ROMol* targetMolecule,
      const MCSParameters& parameters,
      MCSCompareFunctionsData& compareFunctionsData);

 private:
  void makeRingGraph(FMCS::Graph& g, const INT_VECT& ring,
    const ROMol* mol) const;

  inline unsigned getQueryRingIndex(const INT_VECT* r1) const {
    std::map<const INT_VECT*, unsigned>::const_iterator i =
        QueryRingIndex.find(r1);
    if (QueryRingIndex.end() == i) throw - 1;  // never
    return i->second;
  }
  inline const RingMatchTable& getTargetMatchMatrix(const ROMol* mol2) const {
    std::map<const ROMol*, RingMatchTable>::const_iterator mi =
        MatchMatrixSet.find(mol2);
    if (MatchMatrixSet.end() == mi) throw - 1;  // never
    return mi->second;
  }

  inline RingMatchTable& addTargetMatchMatrix(const ROMol* mol2, unsigned s1,
                                              unsigned s2) {
    RingMatchTable& m = MatchMatrixSet[mol2];
    m.resize(s1, s2);
    m.makeRingIndex(mol2);
    return m;
  }
};
}  // namespace FMCS
}  // namespace RDKit
