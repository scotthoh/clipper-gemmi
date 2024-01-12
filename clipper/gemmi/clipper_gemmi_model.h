/* clipper_gemmi.cpp: Gemmi wrapper */
// Author: Soon Wen Hoh, 2023, University of York

#ifndef CLIPPER_GEMMI_MODEL
#define CLIPPER_GEMMI_MODEL

#include "../core/coords.h"
#include "clipper_gemmi.h"

#include <gemmi/model.hpp> // for Gemmi's model data structures
#include <gemmi/to_pdb.hpp>
#include <gemmi/to_mmcif.hpp>
#include <gemmi/to_cif.hpp>

// #include <clipper/clipper.h>

namespace clipper
{

  namespace gemmi
  {
    typedef ::gemmi::Structure CGStructure;
    typedef ::gemmi::CoorFormat CGCoorFormat;

    typedef ::gemmi::Model CGModel;
    typedef ::gemmi::Chain CGChain;
    typedef ::gemmi::Residue CGResidue;
    typedef ::gemmi::Atom CGAtom;
    typedef ::gemmi::CRA CGCRA;
    typedef ::gemmi::CraProxy CGCraProxy;
    typedef ::gemmi::ResidueId CGRId;

    typedef ::gemmi::Model *PGModel;
    typedef ::gemmi::Chain *PGChain;
    typedef ::gemmi::Residue *PGResidue;
    typedef ::gemmi::Atom *PGAtom;

    typedef ::gemmi::PdbReadOptions CGPdbReadOptions;
    typedef ::gemmi::PdbWriteOptions CGPdbWriteOptions;
    typedef ::gemmi::MmcifOutputGroups CGMmcifOutputGroups;
    typedef ::gemmi::cif::Style CGCifStyle;
    typedef ::gemmi::cif::Document CGCifDocument;
  }

  //! GEMMI Atom object wrapper
  /*! This class is a trivial derivation of the corresponding GEMMI,
    providing access in terms of Clipper types. Thus, when you need
    such access, simply cast the GEMMI object reference to this type
    to access the additional functions. */
  class GemmiAtom : public gemmi::CGAtom
  {
  public:
    //! null constructor
    GemmiAtom() {}
    //! constructor: from GEMMI atom
    GemmiAtom(const gemmi::CGAtom &a) : gemmi::CGAtom(a) {}

    // standard atom properties
    String id() const;                            //!< Atom id, e.g. CA, CB, CH3
    String element() const;                       //!< Atom element, e.g. C, H, Zn2+
    Coord_orth coord_orth() const;                //!< Atom coordinate (orthogonal Angstroms)
    ftype occupancy() const;                      //!< Atom occupancy (0...1)
    ftype u_iso() const;                          //!< Atom isotropic U
    U_aniso_orth u_aniso_orth() const;            //!< Atom anisotropic U (orthogonal As)
    void set_id(const String &n);                 //!< set id
    void set_element(const String &n);            //!< set element
    void set_coord_orth(const Coord_orth &v);     //!< set coordinate
    void set_occupancy(const ftype &v);           //!< set occupancy
    void set_u_iso(const ftype &v);               //!< set iso U
    void set_u_aniso_orth(const U_aniso_orth &v); //!< set aniso U
    // other atom properties
    String altconf() const;            //!< get atom alternate conformation code
    void set_altconf(const String &n); //!< set altconf
    int serial_num() const;            //!< get atom serial number
    String charge() const;             //!< get atom charge
  };

  //! GEMMI residue object wrapper
  /*! This class is a trivial derivation of the corresponding GEMMI,
      providing access in terms of Clipper types. Thus, when you need
      such access, simply cast the GEMMI object reference to this type
      to access the additional functions. */
  class GemmiResidue : public gemmi::CGResidue
  {
  public:
    //! null constructor
    GemmiResidue() {}
    //! constructor: from GEMMI residue
    GemmiResidue(const gemmi::CGResidue &r) : gemmi::CGResidue(r) {}

    // residue properties
    String type() const;               //!< Residue type
    int seqnum() const;                //!< Sequence number
    String inscode() const;            //!< Insertion code
    void set_type(const String &n);    //!< Set residue type
    void set_seqnum(const int &n);     //!< Set sequence number
    void set_inscode(const String &n); //!< Set insertion code
  };

  //! GEMMI chain object wrapper
  /*! This class is a trivial derivation of the corresponding GEMMI,
    providing access in terms of Clipper types. Thus, when you need
    such access, simply cast the GEMMI object reference to this type
    to access the additional functions. */
  class GemmiChain : public gemmi::CGChain
  {
  public:
    //! null constructor
    // GemmiChain() {}
    //! constructor: from GEMMI chain
    GemmiChain(const gemmi::CGChain &c) : gemmi::CGChain(c) {}

    // chain properties
    String id() const;            //!< Chain id
    void set_id(const String &n); //!< Set chain id
  };

  //! GEMMI model object wrapper
  /*! This class is a trivial derivation of the corresponding GEMMI,
    providing access in terms of Clipper types. Thus, when you need
    such access, simply cast the GEMMI object reference to this type
    to access the additional functions. */
  class GemmiModel : public gemmi::CGModel
  {
  public:
    //! null constructor
    // GemmiModel() {}
    //! constructor: from GEMMI model
    GemmiModel(const gemmi::CGModel &m) : gemmi::CGModel(m) {}

    // model properties
    String id() const;
    void set_id(const String &n);
  };

  //! Gemmi Structure object wrapper
  /*! This class is a trivial derivation of the corresponding Gemmi,
    providing access in terms of Clipper types. Thus, when you need
    such access, simply cast the Gemmi object reference to this type
    to access additional functions. */
  class GemmiStructure : public gemmi::CGStructure
  {
  public:
    //! null constructor
    GemmiStructure() {}
    //! copy constructor: from GEMMI Structure
    GemmiStructure(const gemmi::CGStructure &s) : gemmi::CGStructure(s) {}

    //! get IDs for model, chain, residue or atom
    String GetID_str(const String &model_name, const gemmi::CGCRA &cra, const String &hnd);
    //! get spacegroup and return in clipper::Spacegroup format
    Spacegroup spacegroup() const;
    //! get cell and return in clipper::Cell format
    Cell get_cell() const;
    //! set spacegroup to Gemmi Structure
    void set_spacegroup(const Spacegroup &spacegroup);
    //! set cell to Gemmi Structure
    void set_cell(const Cell &cell);
    //! print experimental method
    String get_exptlmethod() const;
    //! get transformation/origx in clipper::RTop_orth format: ORIGXn; _database_PDB_matrix.origx*
    RTop_orth get_transform() const;
    //! set transformation/origx to Gemmi Structure
    void set_transform(const RTop_orth &rtop);
    //! get resolution in clipper::Resolution format
    Resolution get_resolution() const;
    //! set resolution to Gemmi Structure
    void set_resolution(const Resolution &resol);
  };

  //! GEMMI atom list class
  /*! This class is used to convert GEMMI CRA to Clipper Atom_list
   It is a trivial derivation of clipper::Atom_list, and may be used
   wherever an Atom_list is required. */
  class GemmiAtom_list : public Atom_list
  {
  public:
    //! constructor: from CRA
    GemmiAtom_list(gemmi::CGCraProxy cra, const int natom);
  };

  struct GemmiPdbReadOptions : public gemmi::CGPdbReadOptions
  {
  public:
    //! null constructor
    GemmiPdbReadOptions() {}
    //! copy constructor: from Gemmi PdbReadOptions
    GemmiPdbReadOptions(gemmi::CGPdbReadOptions o) : gemmi::CGPdbReadOptions(o) {}
    //! options
    bool force_label = false; //!< assign label_seq even if full sequence is not known (assumes no gaps)
    bool copy_remarks = true; //!< to keep or discard resolution remark and REMARK 350
  };

  struct GemmiPdbWriteOptions : public gemmi::CGPdbWriteOptions
  {
  public:
    //! null constructor
    GemmiPdbWriteOptions() {}
    //! copy constructor: from Gemmi PdbWriteOptions
    GemmiPdbWriteOptions(gemmi::CGPdbWriteOptions o) : gemmi::CGPdbWriteOptions(o) {}
  };

  struct GemmiMmcifOutputGroups : public gemmi::CGMmcifOutputGroups
  {
  public:
    //! null constructor
    // GemmiMmcifOutputGroups() {}
    //! copy constructor: from Gemmi MmcifOutputGroups
    GemmiMmcifOutputGroups(gemmi::CGMmcifOutputGroups o) : gemmi::CGMmcifOutputGroups(o) {}
    GemmiMmcifOutputGroups(bool o) : gemmi::CGMmcifOutputGroups(o) {}
  };

  //! Write options for writing PDB/CIF file
  /*! A few simplified quick options for writing PDB/CIF files.
      Options such as Minimal and AllAuth will always overwrite
      what is set in GEMMIFile::set_mmcif_output_groups
      if they contradicts. */
  struct GemmiWriteOptions
  {
    bool Minimal = false;         //!< disable many records (HEADER, TITLE, ...) for PDB and mmcif
    bool ShortTer = false;        //!< write PDB TER records without numbers
    bool LinkR = false;           //!< use non-standard Refmac LINKR record instead of LINK
    bool AllAuth = false;         //!< write _atom_site.auth_atom_id (same as label_atom_id) and auth_comp_id (same as label_comp_id) for mmcif
    bool ShortChainNames = false; //!< shorten chain names, useful for PDB output if chain name size > 2
    bool UpdateCifDoc = false;    //!< update existing mmcif blocks from read in document
  };
}

#endif