/*******************************************************************************

Copyright 2024 by Greg Landrum and the pubchem_shape contributors

This file is part of pubchem_shape

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

***********************************************************************/

#include <boost/python.hpp>
#include <vector>

#include "../rdkit_shape.hpp"
#include "../shape_functions.hpp"
#include <GraphMol/RDKitBase.h>

namespace python = boost::python;

namespace helpers {
python::tuple alignMol(const RDKit::ROMol &ref, RDKit::ROMol &probe,
                       int refConfId, int probeConfId, bool useColors) {
  std::vector<float> matrix(12, 0.0);
  auto [nbr_st, nbr_ct] =
      AlignMolecule(ref, probe, matrix, refConfId, probeConfId, useColors);
  return python::make_tuple(nbr_st, nbr_ct);
}
} // namespace helpers

void wrap_pubchemshape() {
  python::def("AlignMol", &helpers::alignMol,
              (python::arg("ref"), python::arg("probe"),
               python::arg("refConfId") = -1, python::arg("probeConfId") = -1,
               python::arg("useColors") = true),
              "aligns probe to ref, probe is modified");
}

BOOST_PYTHON_MODULE(cshapealign) { wrap_pubchemshape(); }