/*! \file minimol_io.h
  Header file for atomic model io types */

// C Copyright (C) 2000-2006 Kevin Cowtan and University of York
// L
// L  This library is free software and is distributed under the terms
// L  and conditions of version 2.1 of the GNU Lesser General Public
// L  Licence (LGPL) with the following additional clause:
// L
// L     `You may also combine or link a "work that uses the Library" to
// L     produce a work containing portions of the Library, and distribute
// L     that work under terms of your choice, provided that you give
// L     prominent notice with each copy of the work that the specified
// L     version of the Library is used in it, and that you include or
// L     provide public access to the complete corresponding
// L     machine-readable source code for the Library including whatever
// L     changes were used in the work. (i.e. If you make changes to the
// L     Library you must distribute those, but you do not need to
// L     distribute source or object code to those portions of the work
// L     not covered by this licence.)'
// L
// L  Note that this clause grants an additional right and does not impose
// L  any additional restriction, and so does not affect compatibility
// L  with the GNU General Public Licence (GPL). If you wish to negotiate
// L  other terms, please contact the maintainer.
// L
// L  You can redistribute it and/or modify the library under the terms of
// L  the GNU Lesser General Public License as published by the Free Software
// L  Foundation; either version 2.1 of the License, or (at your option) any
// L  later version.
// L
// L  This library is distributed in the hope that it will be useful, but
// L  WITHOUT ANY WARRANTY; without even the implied warranty of
// L  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// L  Lesser General Public License for more details.
// L
// L  You should have received a copy of the CCP4 licence and/or GNU
// L  Lesser General Public License along with this library; if not, write
// L  to the CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
// L  The GNU Lesser General Public can also be obtained by writing to the
// L  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// L  MA 02111-1307 USA

#ifndef CLIPPER_MINIMOL_IO
#define CLIPPER_MINIMOL_IO

#include "../mmdb/clipper_mmdb.h"
#include "minimol_seq.h"
#include "../gemmi/clipper_gemmi_model.h"
#include <gemmi/mmread.hpp>
#include <gemmi/select.hpp>
#include <gemmi/fstream.hpp>
#include <gemmi/to_json.hpp>
#include <gemmi/align.hpp>
#include <gemmi/assembly.hpp>

namespace clipper
{

  //! MMDB file object for MiniMol i/o
  /*! This object is an i/o object for MiniMol, representing an
    interface between MiniMol and PDB or CIF files. It is implemented
    as an MMDB Manager, so MMDB methods may be additionally used on
    the MMDBfile object. */
  class MMDBfile : public MMDBManager
  {
  public:
    //! load MMDB hierarchy from file
    void read_file(const String &file);
    //! save MMDB hierarchy to file
    void write_file(const String &file, TYPE type = Default);
    //! import MiniMol from MMDB hierarchy
    void import_minimol(MiniMol &minimol, const int hnd = -1);
    //! export MiniMol to MMDB hierarchy
    void export_minimol(const MiniMol &minimol, TYPE type = Default);
  };

  //! SEQ file object for MiniMol sequence i/o
  class SEQfile
  {
  public:
    //! load SEQ data from file
    void read_file(const String &file);
    //! read a single sequence from the SEQ file
    void import_polymer_sequence(MPolymerSequence &target);
    //! read a molecule from the SEQ file
    void import_molecule_sequence(MMoleculeSequence &target);

  private:
    String contents;
  };

  //! GEMMI file object for MiniMol i/o
  /*! This object is an i/o object for MiniMol, representing an
    interface between MiniMol and PDB or CIF files. The model
    is read using GEMMI implemented methods and then imported to
    MiniMol object and vice versa. */
  class GEMMIFile
  {
  public:
    enum TYPE
    {
      Default = -1,
      PDB,
      CIF,
      Mmjson
    };
    enum CifStyle
    {
      Aligned,     // ::gemmi::cif::Style::Aligned
      Pdbx,        // ::gemmi::cif::Style::Pdbx
      PreferPairs, // ::gemmi::cif::Style::PreferPairs
    };

    //! null constructor
    GEMMIFile();
    //! read structure from file
    void read_file(const String &file, GemmiPdbReadOptions pdb_read_options = GemmiPdbReadOptions());
    //! write structure to file, minimal
    void write_file(const String &file, TYPE type = Default, GemmiWriteOptions write_options = GemmiWriteOptions());
    //! import MiniMol from Gemmi Structure
    void import_minimol(MiniMol &minimol, const int hnd = -1);
    //! export MiniMol to Gemmi Structure
    void export_minimol(MiniMol &minimol);
    //! set structure if not read from file
    void set_gemmi_structure(const GemmiStructure &st) { structure_ = st; };
    //! get structure held in GEMMIFile
    GemmiStructure get_gemmi_structure() const { return structure_; };
    //! set mmcif output groups, default: GemmiMmcifOutputGroups groups(true)
    void set_mmcif_output_groups(GemmiMmcifOutputGroups groups) { groups_ = groups; };
    //! set pdb write options
    void set_pdb_write_options(GemmiPdbWriteOptions opts) { pdb_write_opts_ = opts; };
    //! set spacegroup to structure
    void set_spacegroup(const Spacegroup &sg);
    //! set cell to structure
    void set_cell(const Cell &cell_in);
    //! get spacegroup in clipper::Spacegroup format
    Spacegroup spacegroup() const { return structure_.spacegroup(); };
    //! get cell in clipper::Cell format
    Cell cell() const { return structure_.get_cell(); };
    //! set cif output style
    void set_cif_style(CifStyle cifstyle);

  private:
    GemmiStructure structure_;
    gemmi::CGCifDocument st_doc_; 
    gemmi::CGCifStyle cifstyle_; // = gemmi::CGCifStyle::PreferPairs;
    GemmiMmcifOutputGroups groups_ = GemmiMmcifOutputGroups(true);
    GemmiPdbWriteOptions pdb_write_opts_ = GemmiPdbWriteOptions();
  };

} // namespace clipper

#endif
