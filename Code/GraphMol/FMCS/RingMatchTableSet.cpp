#include "SubstructMatchCustom.h"
#include "RingMatchTableSet.h"

namespace RDKit {
namespace FMCS {

void RingMatchTableSet::init(const ROMol* query){
    MatchMatrixSet.clear();
    // fill out QueryRingIndex
    unsigned i = 0;
    const RingInfo::VECT_INT_VECT& rings = query->getRingInfo()->bondRings();
    for (RingInfo::VECT_INT_VECT::const_iterator r = rings.begin();
         r != rings.end(); r++)
      QueryRingIndex[&*r] = i++;
    TargetBondRingsIndecesSet.clear();
    QueryBondRingsIndeces = &TargetBondRingsIndecesSet[query];
    QueryBondRingsIndeces->resize(query->getNumBonds());

    size_t ri = 0;
    for (RingInfo::VECT_INT_VECT::const_iterator r = rings.begin();
         r != rings.end(); r++, ri++)
      for (INT_VECT::const_iterator bi = r->begin(); bi != r->end();
           bi++)  // all bonds in the ring
        (*QueryBondRingsIndeces)[*bi].push_back(ri);
}

void RingMatchTableSet::computeRingMatchTable(
      const ROMol* query, const ROMol* targetMolecule,
      const MCSParameters& parameters,
      MCSCompareFunctionsData& compareFunctionsData) {// call it for all targets
    const RingInfo::VECT_INT_VECT& rings1 = query->getRingInfo()->bondRings();
    const RingInfo::VECT_INT_VECT& rings2 =
        targetMolecule->getRingInfo()->bondRings();
    RingMatchTable& m =
        addTargetMatchMatrix(targetMolecule, rings1.size(), rings2.size());
    unsigned i = 0;
    // for each query ring
    for (RingInfo::VECT_INT_VECT::const_iterator r1 = rings1.begin();
         r1 != rings1.end(); r1++, i++) {
      FMCS::Graph graph1;
      makeRingGraph(graph1, *r1,
                    query);  // for each query ring bond ADD all atoms and bonds

      // for each TARGET ring
      for (RingInfo::VECT_INT_VECT::const_iterator r2 = rings2.begin();
           r2 != rings2.end(); r2++) {
        if (r1->size() != r2->size())  // rings are different
          continue;
        FMCS::Graph graph2;
        makeRingGraph(
            graph2, *r2,
            targetMolecule);  // for each TAG ring bond ADD all atoms and bonds

        // check ring substruct match
        MCSBondCompareParameters bp = parameters.BondCompareParameters;
        bp.RingMatchesRingOnly = false;
        bp.CompleteRingsOnly = false;
        bool match =
#ifdef NEVER_xxx_PRECOMPUTED_TABLES_MATCH  // not computed yet, because
                                           // MatchTable computation uses this
                                           // ring info table
            FMCS::SubstructMatchCustomTable(graph2, graph1, tag->AtomMatchTable,
                                            tag->BondMatchTable);
#else  // noticeable slowly:
            FMCS::SubstructMatchCustom(
                graph2, *targetMolecule, graph1, *query, parameters.AtomTyper,
                parameters.BondTyper, nullptr, parameters.AtomCompareParameters,
                bp, compareFunctionsData);
#endif
        if (match) m.setMatch(i, &*r2);
      }
    }
  }

void RingMatchTableSet::makeRingGraph(FMCS::Graph& g, const INT_VECT& ring,
                     const ROMol* mol) const {  // ADD all atoms and bonds
    std::map<const Atom*, unsigned> atomMap;

    for (size_t i = 0; i < ring.size(); i++) {
      const Bond* bond = mol->getBondWithIdx(ring[i]);
      const Atom* atom1 = bond->getBeginAtom();
      const Atom* atom2 = bond->getEndAtom();
      unsigned j1 = NotSet;
      unsigned j2 = NotSet;
      std::map<const Atom*, unsigned>::const_iterator ai;
      ai = atomMap.find(atom1);
      if (atomMap.end() != ai) j1 = ai->second;
      ai = atomMap.find(atom2);
      if (atomMap.end() != ai) j2 = ai->second;
      if (NotSet == j1) {
        j1 = g.m_vertices.size();
        atomMap[atom1] = j1;
        g.addAtom(atom1->getIdx());
      }
      if (NotSet == j2) {
        j2 = g.m_vertices.size();
        atomMap[atom2] = j2;
        g.addAtom(atom2->getIdx());
      }
      g.addBond(ring[i], j1, j2);
    }
  }

}
}
