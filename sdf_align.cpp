#include <iostream>
#include <sstream>
#include <stdexcept>

#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/Fingerprints/FingerprintUtil.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/RWMol.h>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/flyweight.hpp>
#include <boost/flyweight/key_value.hpp>
#include <boost/flyweight/no_tracking.hpp>
#include <RDGeneral/BoostEndInclude.h>


#include "shape_functions.hpp"

using namespace std;

#define ERRORTHROW(msg_stream)                                                 \
  do {                                                                         \
    ostringstream os;                                                          \
    os << msg_stream;                                                          \
    throw runtime_error(os.str());                                             \
  } while (0)

#define DEBUG_MSG(msg_stream) cout << msg_stream << '\n'
// #define DEBUG_MSG(msg_stream)

using namespace RDKit;

// Bondi radii
//  can find more of these in Table 12 of this publication:
//   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3658832/
const double radius_carbon = 1.70;
const double radius_nitrogen = 1.55;
const double radius_oxygen = 1.52;
const double radius_fluorine = 1.47;
const double radius_silicon = 2.10;
const double radius_phosphorous = 1.80;
const double radius_sulfur = 1.80;
const double radius_chlorine = 1.75;
const double radius_bromine = 1.85;
const double radius_iodine = 1.98;
const double radius_color =
    1.08265; // same radius for all feature/color "atoms"

typedef boost::flyweight<boost::flyweights::key_value<std::string, RDKit::MorganFingerprints::ss_matcher>,
                         boost::flyweights::no_tracking>
    pattern_flyweight;
// Definitions for feature points adapted from:
// Gobbi and Poppinger, Biotech. Bioeng. _61_ 47-54 (1998)
std::vector<std::string> smartsPatterns={
    "[$([N;!H0;v3,v4&+1]),\
$([O,S;H1;+0]),\
n&H1&+0]",                                                  // Donor
    "[$([O,S;H1;v2;!$(*-*=[O,N,P,S])]),\
$([O,S;H0;v2]),\
$([O,S;-]),\
$([N;v3;!$(N-*=[O,N,P,S])]),\
n&H0&+0,\
$([o,s;+0;!$([o,s]:n);!$([o,s]:c:n)])]",                    // Acceptor
//    "[a]",                                                  // Aromatic
//    "[F,Cl,Br,I]",                                          // Halogen
    "[#7;+,\
$([N;H2&+0][$([C,a]);!$([C,a](=O))]),\
$([N;H1&+0]([$([C,a]);!$([C,a](=O))])[$([C,a]);!$([C,a](=O))]),\
$([N;H0&+0]([C;!$(C(=O))])([C;!$(C(=O))])[C;!$(C(=O))])]",  // Basic
    "[$([C,S](=[O,S,P])-[O;H1,-1])]"                        // Acidic
};
std::vector<const ROMol *> *getPh4Patterns() {
  static std::unique_ptr<std::vector<const ROMol *>> patterns;
  if (!patterns) {
    patterns.reset(new std::vector<const ROMol *>());
    for (auto smarts : RDKit::MorganFingerprints::defaultFeatureSmarts) {
      const ROMol *matcher =
          pattern_flyweight(smarts)
              .get()
              .getMatcher();
      CHECK_INVARIANT(matcher, "bad smarts");
      patterns->push_back(matcher);
    }
  }

  return patterns.get();
}

void PrepareConformer(
    ROMol &mol, vector<float> &coord, vector<double> &alpha_vector,
    vector<unsigned int> &atom_type_vector,
    vector<unsigned int> &volumeAtomIndexVector,
    map<unsigned int, vector<unsigned int>> &colorAtomType2IndexVectorMap,
    double &sov, double &sof) {
  coord.clear();
  alpha_vector.clear();
  atom_type_vector.clear();
  volumeAtomIndexVector.clear();
  colorAtomType2IndexVectorMap.clear();
  sov = 0.0;
  sof = 0.0;

  // unpack features (PubChem-specific property from SDF)
  // NOTE: this unpacking assumes that RWMol-atom-index = SDF-atom-number - 1
  //   e.g. RWMol uses [0..N-1] and SDF uses [1..N], with atoms in the same
  //   order

  vector<pair<vector<unsigned int>, unsigned int>> feature_idx_type;

  if (mol.hasProp("PUBCHEM_PHARMACOPHORE_FEATURES")) {

    // regular atoms have type 0; feature "atoms" (features represented by a
    // single point+radius) must have type > 0
    static const map<string, unsigned int> atomTypes = {
        {"acceptor", 1}, {"anion", 2},      {"cation", 3},
        {"donor", 4},    {"hydrophobe", 5}, {"rings", 6},
    };

    string features;
    mol.getProp("PUBCHEM_PHARMACOPHORE_FEATURES", features);
    istringstream iss(features);
    string line;
    unsigned int n = 0;
    while (getline(iss, line)) {

      if (n == 0) {
        feature_idx_type.resize(stoul(line));
      }

      else {
        unsigned int f = n - 1;
        if (f >= feature_idx_type.size())
          ERRORTHROW("Too many features");

        istringstream iss2(line);
        string token;
        unsigned int t = 0;
        while (getline(iss2, token, ' ')) {
          if (t == 0) {
            feature_idx_type[f].first.resize(stoul(token));
          } else if (t <= feature_idx_type[f].first.size()) {
            feature_idx_type[f].first[t - 1] = stoul(token) - 1;
          } else {
            map<string, unsigned int>::const_iterator type =
                atomTypes.find(token);
            if (type == atomTypes.end())
              ERRORTHROW("Invalid feature type");
            feature_idx_type[f].second = type->second;
          }
          ++t;
        }
        if (t != (feature_idx_type[f].first.size() + 2))
          ERRORTHROW("Wrong number of tokens in feature");
      }

      ++n;
    }
    if (n != (feature_idx_type.size() + 1))
      ERRORTHROW("Wrong number of features");

    DEBUG_MSG("# features: " << feature_idx_type.size());
  } else {
#if 0
    const auto patts = getPh4Patterns();
    feature_idx_type.clear();

    for(auto p=0u; p<patts->size();++p ){
        auto matches = SubstructMatch(mol,*patts->at(p));
        for(auto match : matches){
            std::vector<unsigned int> ats;
            for(const auto &pr : match){
                ats.push_back(pr.second);
            }
            feature_idx_type.emplace_back(ats,p+1);
        }
    }
  #endif
  }

  // unpack atoms

  if (mol.getNumConformers() != 1)
    ERRORTHROW("Must have exactly one conformer");
  Conformer &conformer = mol.getConformer();
  if (!conformer.is3D())
    ERRORTHROW("Conformer must be 3D");

  unsigned int nAtoms = mol.getNumAtoms();
  // DEBUG_MSG("num atoms: " << nAtoms);

  unsigned int nHeavyAtoms = mol.getNumHeavyAtoms();
  unsigned int i;
  DEBUG_MSG("num heavy atoms: " << nHeavyAtoms);

  unsigned int nAlignmentAtoms = nHeavyAtoms + feature_idx_type.size();

  vector<double> rad_vector(nAlignmentAtoms);
  atom_type_vector.resize(nAlignmentAtoms);

  unsigned int array_idx = 0;
  RDGeom::Point3D ave;
  for (i = 0; i < nAtoms; ++i) {
    unsigned int Z = mol.getAtomWithIdx(i)->getAtomicNum();
    if (Z > 1) {

      atom_type_vector[array_idx] = 0;

      const RDGeom::Point3D &pos = conformer.getAtomPos(i);
      ave += pos;

      switch (Z) {
      case 6:
        rad_vector[array_idx] = radius_carbon;
        break;
      case 7:
        rad_vector[array_idx] = radius_nitrogen;
        break;
      case 8:
        rad_vector[array_idx] = radius_oxygen;
        break;
      case 9:
        rad_vector[array_idx] = radius_fluorine;
        break;
      case 14:
        rad_vector[array_idx] = radius_silicon;
        break;
      case 15:
        rad_vector[array_idx] = radius_phosphorous;
        break;
      case 16:
        rad_vector[array_idx] = radius_sulfur;
        break;
      case 17:
        rad_vector[array_idx] = radius_chlorine;
        break;
      case 35:
        rad_vector[array_idx] = radius_bromine;
        break;
      case 53:
        rad_vector[array_idx] = radius_iodine;
        break;
      default:
        // FIX: add a lookup from the periodic table here.
        //   the results won't be 100% consistent, but it's better than nothing
        ERRORTHROW("Can't use molecules with element Z=" << Z);
      }

      ++array_idx;
    }
  }

  // translate steric center to origin
  ave /= nHeavyAtoms;
  DEBUG_MSG("steric center: (" << ave << ")");

  coord.resize(nAlignmentAtoms * 3);

  array_idx = 0;
  for (i = 0; i < nAtoms; ++i) {

    // translate all atoms
    RDGeom::Point3D &pos = conformer.getAtomPos(i);
    pos -= ave;

    // but use only non-H for alignment
    if (mol.getAtomWithIdx(i)->getAtomicNum() > 1) {
      coord[array_idx * 3] = pos.x;
      coord[(array_idx * 3) + 1] = pos.y;
      coord[(array_idx * 3) + 2] = pos.z;
      ++array_idx;
    }
  }

  // get feature coordinates - simply the average of coords of all atoms in the
  // feature
  for (i = 0; i < feature_idx_type.size(); ++i) {

    RDGeom::Point3D floc;
    for (unsigned int j = 0; j < feature_idx_type[i].first.size(); ++j) {
      unsigned int idx = feature_idx_type[i].first[j];
      if (idx >= nAtoms || mol.getAtomWithIdx(idx)->getAtomicNum() <= 1)
        ERRORTHROW("Invalid feature atom index");
      const RDGeom::Point3D &pos = conformer.getAtomPos(idx);
      floc += pos;
    }
    floc /= feature_idx_type[i].first.size();
    DEBUG_MSG("feature type " << feature_idx_type[i].second << " (" << floc
                              << ")");

    array_idx = nHeavyAtoms + i;
    coord[array_idx * 3] = floc.x;
    coord[(array_idx * 3) + 1] = floc.y;
    coord[(array_idx * 3) + 2] = floc.z;
    rad_vector[array_idx] = radius_color;
    atom_type_vector[array_idx] = feature_idx_type[i].second;
  }

  Align3D::setAlpha(&(rad_vector[0]), rad_vector.size(), alpha_vector);

  // regular atom self overlap
  Align3D::getVolumeAtomIndexVector(
      &(atom_type_vector[0]), atom_type_vector.size(), volumeAtomIndexVector);
  sov = Align3D::ComputeShapeOverlap(&(coord[0]), alpha_vector,
                                     volumeAtomIndexVector, &(coord[0]),
                                     alpha_vector, volumeAtomIndexVector);
  DEBUG_MSG("sov: " << sov);

  // feature self overlap
  if (feature_idx_type.size() > 0) {
    Align3D::getColorAtomType2IndexVectorMap(&(atom_type_vector[0]),
                                             atom_type_vector.size(),
                                             colorAtomType2IndexVectorMap);
    sof = Align3D::ComputeFeatureOverlap(
        &(coord[0]), alpha_vector, colorAtomType2IndexVectorMap, &(coord[0]),
        alpha_vector, colorAtomType2IndexVectorMap);
    DEBUG_MSG("sof: " << sof);
  }
}

void Neighbor(
    // inputs
    ROMol &ref, ROMol &fit, SDWriter &writer, bool write_ref, double opt_param,
    unsigned int max_preiters, unsigned int max_postiters,
    // outputs
    float *matrix, double &nbr_st, double &nbr_ct) {
  vector<float> ref_coord;
  vector<double> alpha_ref_vector;
  vector<unsigned int> ref_volumeAtomIndexVector;
  map<unsigned int, vector<unsigned int>> ref_colorAtomType2IndexVectorMap;
  double ref_sov;
  double ref_sof;
  vector<unsigned int> ref_atom_type_vector;

  vector<float> fit_coord;
  vector<double> alpha_fit_vector;
  vector<unsigned int> fit_volumeAtomIndexVector;
  map<unsigned int, vector<unsigned int>> fit_colorAtomType2IndexVectorMap;
  double fit_sov;
  double fit_sof;
  vector<unsigned int> fit_atom_type_vector;

  Align3D::setUseCutOff(true);

  DEBUG_MSG("Reference details:");
  PrepareConformer(ref, ref_coord, alpha_ref_vector, ref_atom_type_vector,
                   ref_volumeAtomIndexVector, ref_colorAtomType2IndexVectorMap,
                   ref_sov, ref_sof);

  DEBUG_MSG("Fit details:");
  PrepareConformer(fit, fit_coord, alpha_fit_vector, fit_atom_type_vector,
                   fit_volumeAtomIndexVector, fit_colorAtomType2IndexVectorMap,
                   fit_sov, fit_sof);

  // write out ref and fit, both translated to origin but not yet aligned
  if (write_ref) {
    writer.write(ref);
  }

  set<unsigned int> jointColorAtomTypeSet;
  Align3D::getJointColorTypeSet(
      ref_atom_type_vector.data(), ref_atom_type_vector.size(),
      fit_atom_type_vector.data(), fit_atom_type_vector.size(),
      jointColorAtomTypeSet);
  Align3D::restrictColorAtomType2IndexVectorMap(
      ref_colorAtomType2IndexVectorMap, jointColorAtomTypeSet);
  Align3D::restrictColorAtomType2IndexVectorMap(
      fit_colorAtomType2IndexVectorMap, jointColorAtomTypeSet);

  DEBUG_MSG("Running alignment...");
  Align3D::Neighbor_Conformers(
      ref_coord.data(), alpha_ref_vector, ref_volumeAtomIndexVector,
      ref_colorAtomType2IndexVectorMap, ref_sov, ref_sof, fit_coord.data(),
      alpha_fit_vector, fit_volumeAtomIndexVector,
      fit_colorAtomType2IndexVectorMap, fit_sov, fit_sof,
      !jointColorAtomTypeSet.empty(), max_preiters, max_postiters,
      opt_param, matrix, nbr_st, nbr_ct);

  DEBUG_MSG("Done!");
  DEBUG_MSG("nbr_st: " << nbr_st);
  DEBUG_MSG("nbr_ct: " << nbr_ct);

  // transform fit coords
  Conformer &fit_conformer = fit.getConformer();
  vector<float> orig(fit.getNumAtoms() * 3), transformed(fit.getNumAtoms() * 3);
  unsigned int i;
  for (i = 0; i < fit.getNumAtoms(); ++i) {
    const RDGeom::Point3D &pos = fit_conformer.getAtomPos(i);
    orig[i * 3] = pos.x;
    orig[(i * 3) + 1] = pos.y;
    orig[(i * 3) + 2] = pos.z;
  }

  Align3D::VApplyRotTransMatrix(transformed.data(), orig.data(),
                                fit.getNumAtoms(), matrix);

  for (i = 0; i < fit.getNumAtoms(); ++i) {
    RDGeom::Point3D &pos = fit_conformer.getAtomPos(i);
    pos.x = transformed[i * 3];
    pos.y = transformed[(i * 3) + 1];
    pos.z = transformed[(i * 3) + 2];
  }
  fit.setProp("shape_align_shape_tanimoto", nbr_st);
  fit.setProp("shape_align_color_tanimoto", nbr_ct);

  writer.write(fit);
}

int main(int argc, char **argv) {
  int status = -1;

  try {

    if (argc != 2)
      ERRORTHROW("Usage: sdf_align <ref_conformer.sdf> <fit_conformer.sdf> "
                 "opt_param max_preiters max_postiters");

    SDMolSupplier suppl(argv[1]);
    SDWriter writer("sdf_align.out.sdf");
    unique_ptr<ROMol> ref{suppl[0]};
    if (!ref.get()) {
      ERRORTHROW("Failed to read ref conformer");
    }

    auto opt_param = 0.5;
    auto max_preiters = 3u;
    auto max_postiters = 16u;

    for (auto i = 1u; i < suppl.length(); ++i) {
      std::unique_ptr<ROMol> fit{suppl[i]};
      if (!fit) {
        continue;
      }
      vector<float> matrix(12, 0.0);
      double nbr_st = 0.0;
      double nbr_ct = 0.0;
      Neighbor(*ref, *fit, writer, i == 1, opt_param, max_preiters,
               max_postiters, &(matrix[0]), nbr_st, nbr_ct);
    }

    status = 0;

  } catch (std::exception &e) {
    cerr << "Caught std::exception: " << e.what() << '\n';
  } catch (...) {
    cerr << "Caught unknown exception\n";
  }

  return status;
}
