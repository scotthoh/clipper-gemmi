/* minimol_io.cpp: atomic model types */
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

#include "minimol_io.h"

#include <fstream>
#include <sstream>
extern "C"
{
#include <string.h>
}

namespace clipper
{

  // MMDBfile

  // local types for referencing between MMDB and minimol
  struct SAtom
  {
    mmdb::PCAtom db;
    const MAtom *mm;
  };
  struct SMono
  {
    mmdb::PCResidue db;
    const MMonomer *mm;
    std::vector<SAtom> data;
  };
  struct SPoly
  {
    mmdb::PCChain db;
    const MPolymer *mm;
    std::vector<SMono> data;
  };
  struct SModl
  {
    mmdb::PCModel db;
    const MModel *mm;
    std::vector<SPoly> data;
  };

  /*! The file may be either a PDB or mmCIF file.
    If the spacegroup or cell are not set, they will be taken from the
    file, otherwise the existing values override the ones in the file.
    \param file The filename (or pathname) of the file to read. */
  void MMDBfile::read_file(const String &file)
  {
    int err = ReadCoorFile((char *)file.c_str());
    if (err)
      Message::message(Message_fatal("MMDBfile: read_file error: " + file + " : " + String(err)));
  }

  /*! The output file type will be the same as the file read, otherwise PDB.
    \param file The filename (or pathname) of the file to write.
    \param type 0=PDB, 1=CIF, 2=binary, default=same as input file, or PDB. */
  void MMDBfile::write_file(const String &file, TYPE type)
  {
    const TYPE types[3] = {PDB, CIF, Binary};
    int rtype = GetFileType();
    if (type == Default && rtype >= 0 && rtype <= 2)
      type = types[rtype];
    int err;
    switch (type)
    {
    case Binary:
      err = WriteMMDBF((char *)file.c_str());
      break;
    case CIF:
      err = WriteCIFASCII((char *)file.c_str());
      break;
    case PDB:
    default:
      err = WritePDBASCII((char *)file.c_str());
      break;
    }
    if (err)
      Message::message(Message_fatal("MMDBfile: write_file error: " + file + " : " + String(err)));
  }

  /*! Import data from the MMDB hierarchy into a MiniMol object. All the
    atoms, residues and chains are imported unless a selection handle is
    supplied, in which case only the selected atoms and their containers
    are imported.

    Each data is tagged with an optional Property<String> (see
    clipper::PropertyManager), called "CID", which describes the
    position in the MMDB hierarchy from which it was taken. This may
    (optionally) be used later to put a modified element back into the
    hierarchy.

    Any atoms which have alternate conformation codes will also be given
    an "AltConf" Property<String>.

    For more details, see clipper::PropertyManager::exists_property(),
    clipper::PropertyManager::get_property().

    \param minimol The minimol to import.
    \param hnd (optional) MMDB selection handle. */
  void MMDBfile::import_minimol(MiniMol &minimol, const int hnd)
  {
    // clear the model
    minimol = MiniMol(spacegroup(), cell());

    // make atom, residue, chain selections
    int h_atm = NewSelection();
    int h_res = NewSelection();
    int h_chn = NewSelection();
    int h_mod = NewSelection();
    if (hnd < 0)
    {
      SelectAtoms(h_atm, 0, 0, ::mmdb::SKEY_NEW);
    }
    else
    {
      Select(h_atm, ::mmdb::STYPE_ATOM, hnd, ::mmdb::SKEY_NEW);
    }
    Select(h_res, ::mmdb::STYPE_RESIDUE, h_atm, ::mmdb::SKEY_NEW);
    Select(h_chn, ::mmdb::STYPE_CHAIN, h_atm, ::mmdb::SKEY_NEW);
    Select(h_mod, ::mmdb::STYPE_MODEL, h_atm, ::mmdb::SKEY_NEW);

    // now import objects
    char txt[256];
    MModel &mol = minimol.model();
    mmdb::PCModel p_mod = GetModel(1);
    if (p_mod != NULL)
    {
      for (int c = 0; c < p_mod->GetNumberOfChains(); c++)
      {
        mmdb::PCChain p_chn = p_mod->GetChain(c);
        if (p_chn != NULL)
          if (p_chn->isInSelection(h_chn))
          {
            // import the chain
            MPolymer pol;
            for (int r = 0; r < p_chn->GetNumberOfResidues(); r++)
            {
              mmdb::PCResidue p_res = p_chn->GetResidue(r);
              if (p_res != NULL)
                if (p_res->isInSelection(h_res))
                {
                  // import the residue
                  MMonomer mon;
                  for (int a = 0; a < p_res->GetNumberOfAtoms(); a++)
                  {
                    mmdb::PCAtom p_atm = p_res->GetAtom(a);
                    if (p_atm != NULL)
                      if (p_atm->isInSelection(h_atm))
                        if (!p_atm->Ter)
                        {
                          // import the atom
                          MAtom atm(Atom::null());
                          atm.set_name(p_atm->GetAtomName(), p_atm->altLoc);
                          atm.set_element(p_atm->element);
                          if (p_atm->WhatIsSet & ::mmdb::ASET_Coordinates)
                            atm.set_coord_orth(
                                Coord_orth(p_atm->x, p_atm->y, p_atm->z));
                          if (p_atm->WhatIsSet & ::mmdb::ASET_Occupancy)
                            atm.set_occupancy(p_atm->occupancy);
                          if (p_atm->WhatIsSet & ::mmdb::ASET_tempFactor)
                            atm.set_u_iso(Util::b2u(p_atm->tempFactor));
                          if (p_atm->WhatIsSet & ::mmdb::ASET_Anis_tFac)
                            atm.set_u_aniso_orth(
                                U_aniso_orth(p_atm->u11, p_atm->u22, p_atm->u33,
                                             p_atm->u12, p_atm->u13, p_atm->u23));
                          p_atm->GetAtomID(txt);
                          atm.set_property("CID", Property<String>(String(txt)));
                          if (p_atm->altLoc[0] != '\0')
                            atm.set_property("AltConf",
                                             Property<String>(String(p_atm->altLoc)));
                          mon.insert(atm); // store the atom
                        }
                  }
                  mon.set_seqnum(p_res->GetSeqNum(), String(p_res->GetInsCode()));
                  mon.set_type(p_res->GetResName());
                  p_res->GetResidueID(txt);
                  mon.set_property("CID", Property<String>(String(txt)));
                  pol.insert(mon); // store the residue
                }
            }
            pol.set_id(p_chn->GetChainID());
            p_chn->GetChainID(txt);
            pol.set_property("CID", Property<String>(String(txt)));
            mol.insert(pol); // store the chain
          }
      }
      p_mod->GetModelID(txt);
      mol.set_property("CID", Property<String>(String(txt)));
    }

    // clean up
    DeleteSelection(h_atm);
    DeleteSelection(h_res);
    DeleteSelection(h_chn);
    DeleteSelection(h_mod);
  }

  /*! Export data to the MMDB hierarchy from a MiniMol object. All the
    atoms, residues and chains are exported.

    If any MiniMol object has a "CID" Property<String> (see
    clipper::PropertyManager), then the information from that object
    will be used to update the corresponding object in the MMDB
    hierarchy, if it exists. If there is no such entry in the MMDB
    hierarchy, or if no "CID" Property<String> exists, then a new object
    will be created in the MMDB hierarchy.
    \param type 0=PDB, 1=CIF, 2=binary, default=same as input file, or PDB. */
  void MMDBfile::export_minimol(const MiniMol &minimol, TYPE type)
  {
    // remove blank spaces for nucleic acid residue name
    // when cif format is required (e.g. "  U" to U ) SWH-May22
    const TYPE types[3] = {PDB, CIF, Binary};
    int rtype = GetFileType();
    if (type == Default && rtype >= 0 && rtype <= 2)
      type = types[rtype];
    // export spacegroup/cell
    if (!minimol.spacegroup().is_null())
      set_spacegroup(minimol.spacegroup());
    if (!minimol.cell().is_null())
      set_cell(minimol.cell());

    // create structure for relationships between Minimol and MMDB
    SModl smod;
    clipper::String cid;

    // fill structure
    smod.mm = &(minimol.model());
    smod.db = NULL;
    smod.data.resize(smod.mm->size());
    for (int p = 0; p < smod.data.size(); p++)
    { // loop over chains
      SPoly &spol = smod.data[p];
      spol.mm = &((*smod.mm)[p]);
      spol.db = NULL;
      spol.data.resize(spol.mm->size());
      for (int r = 0; r < spol.data.size(); r++)
      { // loop over residues
        SMono &smon = spol.data[r];
        smon.mm = &((*spol.mm)[r]);
        smon.db = NULL;
        smon.data.resize(smon.mm->size());
        for (int a = 0; a < smon.data.size(); a++)
        { // loop over atoms
          SAtom &satm = smon.data[a];
          satm.mm = &((*smon.mm)[a]);
          satm.db = NULL;
        }
      }
    }

    // make the MMDB references by CID if present
    if (smod.mm->exists_property("CID"))
    {
      cid = dynamic_cast<const Property<String> &>(smod.mm->get_property("CID")).value() + "/A/0/A";
      smod.db = GetModel((char *)cid.c_str());
    }
    if (smod.db != NULL)
      for (int p = 0; p < smod.data.size(); p++)
      { // loop over chains
        SPoly &spol = smod.data[p];
        if (spol.mm->exists_property("CID"))
        {
          cid = dynamic_cast<const Property<String> &>(spol.mm->get_property("CID")).value() + "/0/A";
          spol.db = GetChain((char *)cid.c_str());
        }
        if (spol.db != NULL)
          for (int r = 0; r < spol.data.size(); r++)
          { // loop over residues
            SMono &smon = spol.data[r];
            if (smon.mm->exists_property("CID"))
            {
              cid = dynamic_cast<const Property<String> &>(smon.mm->get_property("CID")).value() + "/A";
              smon.db = GetResidue((char *)cid.c_str());
            }
            if (smon.db != NULL)
              for (int a = 0; a < smon.data.size(); a++)
              { // loop over atoms
                SAtom &satm = smon.data[a];
                if (satm.mm->exists_property("CID"))
                {
                  cid = dynamic_cast<const Property<String> &>(satm.mm->get_property("CID")).value();
                  satm.db = GetAtom((char *)cid.c_str());
                }
              }
          }
      }

    // Now create MMDB objects for anything which is missing
    if (smod.db == NULL)
    {
      smod.db = new mmdb::CModel();
      AddModel(smod.db);
    }
    for (int p = 0; p < smod.data.size(); p++)
    { // loop over chains
      SPoly &spol = smod.data[p];
      if (spol.db == NULL)
      {
        spol.db = new mmdb::CChain();
        smod.db->AddChain(spol.db);
      }
      for (int r = 0; r < spol.data.size(); r++)
      { // loop over residues
        SMono &smon = spol.data[r];
        if (smon.db == NULL)
        {
          smon.db = new mmdb::CResidue();
          spol.db->AddResidue(smon.db);
        }
        for (int a = 0; a < smon.data.size(); a++)
        { // loop over atoms
          SAtom &satm = smon.data[a];
          if (satm.db == NULL)
          {
            satm.db = new mmdb::CAtom();
            smon.db->AddAtom(satm.db);
          }
        }
      }
    }

    // now fill in information in mmdb from MiniMol
    for (int p = 0; p < smod.data.size(); p++)
    { // loop over chains
      SPoly &spol = smod.data[p];
      spol.db->SetChainID((char *)spol.mm->id().substr(0, 9).c_str()); // set id
      for (int r = 0; r < spol.data.size(); r++)
      { // loop over residues
        SMono &smon = spol.data[r];
        smon.db->seqNum = smon.mm->seqnum(); // set residue info
        // remove blank spaces for nucleic acid residue name
        // when cif format is required (e.g. "  U" to U ) SWH-May22
        if (type == CIF)
        {
          smon.db->SetResName((char *)smon.mm->type().trim().substr(0, 19).c_str());
        }
        else
          smon.db->SetResName((char *)smon.mm->type().substr(0, 19).c_str());
        int pos = smon.mm->id().find(":");
        if (pos != String::npos)
          strcpy(smon.db->insCode, smon.mm->id().substr(pos + 1, 9).c_str());
        for (int a = 0; a < smon.data.size(); a++)
        { // loop over atoms
          SAtom &satm = smon.data[a];
          if (!satm.mm->coord_orth().is_null())
          { // set atom coord
            satm.db->x = satm.mm->coord_orth().x();
            satm.db->y = satm.mm->coord_orth().y();
            satm.db->z = satm.mm->coord_orth().z();
            satm.db->WhatIsSet |= ::mmdb::ASET_Coordinates;
          }
          if (!Util::is_nan(satm.mm->occupancy()))
          { // set atom occ
            satm.db->occupancy = satm.mm->occupancy();
            satm.db->WhatIsSet |= ::mmdb::ASET_Occupancy;
          }
          if (!Util::is_nan(satm.mm->u_iso()))
          { // set atom u_iso
            satm.db->tempFactor = Util::u2b(satm.mm->u_iso());
            satm.db->WhatIsSet |= ::mmdb::ASET_tempFactor;
          }
          if (!satm.mm->u_aniso_orth().is_null())
          { // set atom u_aniso
            satm.db->u11 = satm.mm->u_aniso_orth()(0, 0);
            satm.db->u22 = satm.mm->u_aniso_orth()(1, 1);
            satm.db->u33 = satm.mm->u_aniso_orth()(2, 2);
            satm.db->u12 = satm.mm->u_aniso_orth()(0, 1);
            satm.db->u13 = satm.mm->u_aniso_orth()(0, 2);
            satm.db->u23 = satm.mm->u_aniso_orth()(1, 2);
            satm.db->WhatIsSet |= ::mmdb::ASET_Anis_tFac;
          }
          if (satm.mm->id() != "") // atom id
            satm.db->SetAtomName((char *)satm.mm->id().substr(0, 19).c_str());
          if (satm.mm->element() != "") // atom element
            satm.db->SetElementName((char *)satm.mm->element().substr(0, 19).c_str());
          if (satm.mm->exists_property("AltConf"))
          { // alt conf code
            String a = dynamic_cast<const Property<String> &>(satm.mm->get_property("AltConf")).value();
            if (a != "")
              strcpy(satm.db->altLoc, a.substr(0, 19).c_str());
          }
        }
      }
    }
    FinishStructEdit();
  }

  void SEQfile::read_file(const String &file)
  {
    std::ifstream seqfile(file.c_str());
    std::ostringstream s;
    s << seqfile.rdbuf();
    contents = s.str();
    // non-portable to old Sun
    // contents = std::string(std::istreambuf_iterator<char>(seqfile),
    //            std::istreambuf_iterator<char>());
  }

  void SEQfile::import_polymer_sequence(MPolymerSequence &target)
  {
    MMoleculeSequence mms;
    import_molecule_sequence(mms);
    target = mms[0];
  }

  void SEQfile::import_molecule_sequence(MMoleculeSequence &target)
  {
    std::vector<clipper::String> lines = contents.split("\n");
    clipper::String id, seq = "";
    for (int l = 0; l < lines.size(); l++)
    {
      clipper::String line = lines[l].trim();
      if (line[0] == '>')
      {
        if (seq != "")
        {
          MPolymerSequence s;
          s.set_id(id);
          s.set_sequence(seq);
          target.insert(s);
        }
        id = line.substr(1);
        id = id.trim();
        seq = "";
      }
      else if (isalpha(line[0]))
      {
        for (int i = 0; i < line.length(); i++)
          if (isalpha(line[i]))
            seq += toupper(line[i]);
      }
    }
    if (seq != "")
    {
      MPolymerSequence s;
      s.set_id(id);
      s.set_sequence(seq);
      target.insert(s);
    }
  }

  // local types for referencing between GEMMI and MiniMol
  // work for GEMMI >= 0.6.3
  struct GAtom
  {
    gemmi::PGAtom db;
    const MAtom *mm;
  };
  struct GMono
  {
    gemmi::PGResidue db;
    const MMonomer *mm;
    std::vector<GAtom> data;
  };
  struct GPoly
  {
    gemmi::PGChain db;
    const MPolymer *mm;
    std::vector<GMono> data;
  };
  struct GModl
  {
    gemmi::PGModel db;
    const MModel *mm;
    std::vector<GPoly> data;
  };

  /*! null constructor, initialise cifstyle_*/
  GEMMIFile::GEMMIFile()
  {
    cifstyle_ = gemmi::CGCifStyle::PreferPairs;
  }

  /*! This reads file of either PDB, mmCIF or mmjson format.
    Input file will be determined from extension. i.e. ".pdb", ".ent", ".cif", ".mmcif", ".json"
    \param file The filename (or pathname) of the file to read.
    \param pdb_read_opts PdbReadOptions struct
    \parblock
    max_line_length=0, split_chain_on_ter=false, skip_remarks=false, force_label=false, copy_remarks=true
    \endparblock */
  void GEMMIFile::read_file(const String &file, GemmiPdbReadOptions pdb_read_opts)
  {
    gemmi::CGCoorFormat format = ::gemmi::coor_format_from_ext(file.trim().tail());
    switch (format)
    {
    case gemmi::CGCoorFormat::Pdb:
    {
      structure_ = ::gemmi::read_pdb(::gemmi::BasicInput(file), pdb_read_opts);
      // keep or clear remarks read, i.e. resolution and REMARK 350
      if (!pdb_read_opts.copy_remarks)
        structure_.raw_remarks.clear();
      // setup entities for mmcif format
      setup_entities(structure_);
      // Uses sequence alignment (model to SEQRES) to assign label_seq used in mmcif format
      // force_label=true will assign label_seq even if full sequence is not known (assumes no gaps)
      // refer GEMMI documentation
      assign_label_seq_id(structure_, pdb_read_opts.force_label);
      // pdb_input_routines(pdb_read_opts);
      break;
    }
    case gemmi::CGCoorFormat::Mmcif:
    {
      st_doc_.clear();
      structure_ = ::gemmi::make_structure(::gemmi::cif::read(::gemmi::BasicInput(file)), &st_doc_);
      break;
    }
    case gemmi::CGCoorFormat::Mmjson:
    {
      st_doc_.clear();
      structure_ = ::gemmi::make_structure(::gemmi::cif::read_mmjson(::gemmi::BasicInput(file)), &st_doc_);
      break;
    }
    case gemmi::CGCoorFormat::Unknown:
    case gemmi::CGCoorFormat::Detect:
    {
      String msg = "GEMMIFile: Error reading file, unknown format: " + file + "\n";
      Message::message(Message_fatal(msg));
    }
    }
  }

  /*! The default output file type will be the same as the file read, otherwise PDB.
    write_options.Minimal = true will always take precedence even when
    GemmiMmcifOutputGroups is set through GEMMIFile::set_mmcif_output_groups or
    PdbWriteOptions is set through GEMMIFile::set_pdb_write_options.
    \param file The filename (or pathname) of the file to be written.
    \param type 0=PDB, 1=CIF, default=same as input file or PDB
    \param write_options
    \parblock
    GemmiWriteOptions Minimal=false, ShortTER=false, LinkR=false, AllAuth=false,
    ShortChainNames = false, UpdateCifDoc = false;
    \endparblock */
  void GEMMIFile::write_file(const String &file, TYPE type, GemmiWriteOptions write_options)
  {
    // GEMMI file types
    // { Unknown, Detect, Pdb, Mmcif, Mmjson, ChemComp }
    const TYPE types[4] = {PDB, CIF, Mmjson};
    auto intype = structure_.input_format;
    if (type == Default)
    {
      switch (intype)
      {
      case gemmi::CGCoorFormat::Mmcif:
        type = types[1];
        break;
      case gemmi::CGCoorFormat::Mmjson:
        type = types[2];
        break;
      case gemmi::CGCoorFormat::Pdb:
      default:
        type = types[0];
      }
    }
    gemmi::CGCifDocument doc = st_doc_;
    ::gemmi::Ofstream os(file, &std::cout);
    switch (type)
    {
    case CIF: // write to mmcif/mmjson refer Gemmi's convert program
    case Mmjson:
    {
      if (doc.blocks.size() == 0)
      {
        doc.clear();
        doc = ::gemmi::make_mmcif_document(structure_);
      }
      if (write_options.Minimal)
      {
        doc.clear();
        doc.blocks.resize(1);
        ::gemmi::add_minimal_mmcif_data(structure_, doc.blocks[0]);
        GemmiMmcifOutputGroups groups(false);
        groups.atoms = true;
        groups.group_pdb = true;
        ::gemmi::update_mmcif_block(structure_, doc.blocks[0], groups);
      }
      else
      {
        if (!write_options.UpdateCifDoc)
        {
          // UpdateCifDoc = false so update only atoms
          // this will overwrite set_mmcif_output_groups
          GemmiMmcifOutputGroups groups(false);
          groups.atoms = true;
          groups.group_pdb = true;
          groups.auth_all = write_options.AllAuth;
          ::gemmi::update_mmcif_block(structure_, doc.blocks[0], groups);
        }
        else
        {
          // cif groups in document can be updated individually if set to true
          // i.e. groups.symmetry = true
          // use GEMMIFile::set_mmcif_output_groups to set groups
          // refer ::gemmi::MmcifOutputGroups or gemmi's to_mmcif.hpp
          ::gemmi::update_mmcif_block(structure_, doc.blocks[0], groups_);
        }
      }
      if (type == CIF)
      {
        ::gemmi::cif::write_cif_to_stream(os.ref(), doc, cifstyle_);
      }
      else
      {
        ::gemmi::cif::JsonWriter jswriter(os.ref());
        jswriter.set_mmjson();
        jswriter.write_json(doc);
      }
    }
    break;
    case PDB:
    default:
    {
      ::gemmi::shorten_ccd_codes(structure_);
      if (write_options.Minimal)
        pdb_write_opts_ = GemmiPdbWriteOptions::minimal();
      if (write_options.ShortTer)
        pdb_write_opts_.numbered_ter = false;
      if (write_options.LinkR)
        pdb_write_opts_.use_linkr = true;
      if (write_options.ShortChainNames)
        shorten_chain_names(structure_);
      try
      {
        ::gemmi::write_pdb(structure_, os.ref(), pdb_write_opts_);
      }
      catch (const std::exception &err)
      {
        String msg = "GEMMIFile: write file error: " + file + "\n>> ";
        msg += err.what();
        Message::message(Message_fatal(msg));
      }
    }
    break;
    }
  }

  /*! Import data from the GEMMI hierarchy into a MiniMol object. All the
    atoms, residues and chains are imported unless a selection handle is
    supplied, in which case only the selected atoms and their containers
    are imported.

    Each data is tagged with an optional Property<String> (see
    clipper::PropertyManager), called "CID", which describes the
    position in the GEMMI hierarchy from which it was taken. This may
    (optionally) be used later to put a modified element back into the
    hierarchy.

    Any atoms which have alternate conformation codes all also be given
    an "AltConf" Property<String>.

    For more details, see clipper::PropertyManager::exists_property(),
    clipper::PropertyManager::get_property(). (e.g. label = "CID")
    \param minimol The minimol to import.
    \param sel_hnd (optional) Model selection handle when number of models > 1. */
  void GEMMIFile::import_minimol(MiniMol &minimol, const int sel_hnd)
  {
    // clear the model
    if (!structure_.spacegroup().is_null())
      minimol = MiniMol(structure_.spacegroup(), structure_.get_cell());
    else
      minimol = MiniMol(Spacegroup::null(), structure_.get_cell());

    // import objects
    String txt;
    MModel &mol = minimol.model();

    // get name of first model or model number from sel_hnd
    int model_num = 0;
    if (sel_hnd >= 1)
      model_num = sel_hnd - 1;

    std::string model_name = structure_.models[model_num].name;
    gemmi::PGModel p_mod = structure_.find_model(model_name);

    mol.set_property("StrucName", Property<String>(structure_.name));
    gemmi::CGCRA cra;

    if (p_mod != NULL)
    {
      for (int c = 0; c < p_mod->chains.size(); c++)
      {
        gemmi::PGChain p_chn = &p_mod->chains.at(c);
        if (p_chn != NULL)
        {
          // import the chain
          cra.chain = p_chn;
          MPolymer pol;
          for (gemmi::CGResidue res : p_chn->whole())
          {
            // import residue
            cra.residue = &res;
            MMonomer mon;
            for (gemmi::CGAtom a : res.children())
            {
              // import atoms
              cra.atom = &a;
              MAtom atm(Atom::null());
              String altLoc = a.has_altloc() ? String(&a.altloc, 1) : "";
              atm.set_name(a.padded_name(), altLoc, true);
              atm.set_element(a.element.uname());
              atm.set_coord_orth(Coord_orth(a.pos.at(0), a.pos.at(1), a.pos.at(2)));
              atm.set_occupancy(a.occ);
              atm.set_u_iso(Util::b2u(a.b_iso));
              if (a.aniso.u11 != 0)
                atm.set_u_aniso_orth(U_aniso_orth(a.aniso.u11, a.aniso.u22, a.aniso.u33,
                                                  a.aniso.u12, a.aniso.u13, a.aniso.u23));
              txt = structure_.GetID_str(p_mod->name, cra, a.what()); // Get atom ID
              atm.set_property("CID", Property<String>(txt));
              if (a.has_altloc())
                atm.set_property("AltConf", Property<String>(altLoc));
              if (a.charge)
                atm.set_property("Charge", Property<signed char>(a.charge));
              mon.insert(atm); // store the atom
            }
            cra.atom = nullptr;
            String inscode = res.seqid.has_icode() ? String(&res.seqid.icode, 1) : "";
            mon.set_seqnum(res.seqid.num.value, inscode); // String(&res.seqid.icode, 1));
            mon.set_type(res.name);
            txt = structure_.GetID_str(p_mod->name, cra, res.what()); // Get residue ID
            mon.set_property("CID", Property<String>(String(txt)));
            pol.insert(mon); // store the residue
          }

          pol.set_id(p_chn->name);
          cra.residue = nullptr;
          txt = structure_.GetID_str(p_mod->name, cra, p_chn->what()); // Get chain ID
          pol.set_property("CID", Property<String>(String(txt)));
          mol.insert(pol); // store the chain
        }
      }
      cra.chain = nullptr;
      txt = structure_.GetID_str(p_mod->name, cra, p_mod->what()); // Get model ID
      mol.set_property("CID", Property<String>(String(txt)));
    }
  }

  /*! Export data to the GEMMI hierachy from a MiniMol object*/
  void GEMMIFile::export_minimol(MiniMol &minimol)
  {
    // export spacegroup/cell
    if (!minimol.spacegroup().is_null())
      structure_.set_spacegroup(minimol.spacegroup());
    if (!minimol.cell().is_null())
      structure_.set_cell(minimol.cell());

    // create structure for relationships between MiniMol and GEMMI
    GModl gmod;
    clipper::String cid;

    // fill structure
    gmod.mm = &(minimol.model());
    gmod.db = NULL;
    gmod.data.resize(gmod.mm->size());
    for (int p = 0; p < gmod.data.size(); p++)
    { // loop over chains
      GPoly &gpol = gmod.data[p];
      gpol.mm = &((*gmod.mm)[p]);
      gpol.db = NULL;
      gpol.data.resize(gpol.mm->size());
      for (int r = 0; r < gpol.data.size(); r++)
      { // loop over residues
        GMono &gmon = gpol.data[r];
        gmon.mm = &((*gpol.mm)[r]);
        gmon.data.resize(gmon.mm->size());
        for (int a = 0; a < gmon.data.size(); a++)
        { // loop over atoms
          GAtom &gatm = gmon.data[a];
          gatm.mm = &((*gmon.mm)[a]);
          gatm.db = NULL;
        }
      }
    }

    // make the GEMMI reference by CID if present
    if (gmod.mm->exists_property("CID"))
    {
      cid = dynamic_cast<const Property<String> &>(gmod.mm->get_property("CID")).value(); // + "/A/0/A";
      std::pair<gemmi::PGModel, gemmi::CGCRA> PM_CRA = ::gemmi::Selection(cid.trim()).first(structure_);
      if (PM_CRA.first != nullptr)
        gmod.db = PM_CRA.first;

      if (gmod.db != NULL)
      {
        for (int p = 0; p < gmod.data.size(); p++)
        {
          GPoly &gpol = gmod.data[p];
          if (gpol.mm->exists_property("CID"))
          {
            cid = dynamic_cast<const Property<String> &>(gpol.mm->get_property("CID")).value(); // + "/0/A";
            PM_CRA = ::gemmi::Selection(cid.trim()).first(structure_);
            gpol.db = PM_CRA.second.chain;
          }
          if (gpol.db != NULL)
          {
            for (int r = 0; r < gpol.data.size(); r++)
            {
              GMono &gmon = gpol.data[r];
              if (gmon.mm->exists_property("CID"))
              {
                cid = dynamic_cast<const Property<String> &>(gmon.mm->get_property("CID")).value(); // + "/A";
                PM_CRA = ::gemmi::Selection(cid.trim()).first(structure_);
                gmon.db = PM_CRA.second.residue;
              }
              if (gmon.db != NULL)
              {
                for (int a = 0; a < gmon.data.size(); a++)
                {
                  GAtom &gatm = gmon.data[a];
                  if (gatm.mm->exists_property("CID"))
                  {
                    cid = dynamic_cast<const Property<String> &>(gatm.mm->get_property("CID")).value();
                    PM_CRA = ::gemmi::Selection(cid.trim()).first(structure_);
                    gatm.db = PM_CRA.second.atom;
                  }
                }
              }
            }
          }
        }
      }
    }

    // create GEMMI object for anything that is missing
    // fill information in Gemmi from MiniMol
    if (gmod.db == NULL)
    {
      structure_.models.emplace_back("1");
      gmod.db = &structure_.models.back();
    }
    for (int p = 0; p < gmod.data.size(); p++)
    { // chains
      GPoly &gpol = gmod.data[p];
      if (gpol.db == NULL)
      { // nullptr
        String chn_name = gpol.mm->id().trim();
        gmod.db->chains.emplace_back(chn_name);
        gpol.db = &gmod.db->chains.back();
      }
      else
        gpol.db->name = gpol.mm->id().trim();
      for (int r = 0; r < gpol.data.size(); r++)
      { // residues
        GMono &gmon = gpol.data[r];
        if (gmon.db == NULL) // nullptr
        {
          gemmi::CGRId rid;
          gpol.db->residues.emplace_back(rid);
          gmon.db = &gpol.db->residues.back();
        }
        gmon.db->seqid.num = gmon.mm->seqnum();
        gmon.db->name = gmon.mm->type(); //.substr(0, 19);
        int pos = gmon.mm->id().find(":");
        if (pos != String::npos)
        {
          char inscode = gmon.mm->id()[pos + 1];
          if (inscode != '\r' && inscode != '\n')
            gmon.db->seqid.icode = inscode;
        }
        if (!gmon.db->het_flag)
        {
          bool std_res = ::gemmi::find_tabulated_residue(gmon.db->name).is_standard();
          gmon.db->het_flag = std_res ? 'A' : 'H';
        }
        for (int a = 0; a < gmon.data.size(); a++)
        { // loop through atoms
          GAtom &gatm = gmon.data[a];
          if (gatm.db == NULL) // nullptr
          {
            gemmi::CGAtom atm1;
            gmon.db->atoms.emplace_back(atm1);
            gatm.db = &gmon.db->atoms.back();
          }

          if (!gatm.mm->coord_orth().is_null())
          { // set atom coord
            gatm.db->pos.x = gatm.mm->coord_orth().x();
            gatm.db->pos.y = gatm.mm->coord_orth().y();
            gatm.db->pos.z = gatm.mm->coord_orth().z();
          }
          if (!Util::is_nan(gatm.mm->occupancy())) // set atom occupancy
            gatm.db->occ = gatm.mm->occupancy();
          if (!Util::is_nan(gatm.mm->u_iso())) // set atom u_iso
            gatm.db->b_iso = Util::u2b(gatm.mm->u_iso());
          if (!gatm.mm->u_aniso_orth().is_null())
          { // set atom u_aniso
            gatm.db->aniso.u11 = gatm.mm->u_aniso_orth().mat00();
            gatm.db->aniso.u22 = gatm.mm->u_aniso_orth().mat11();
            gatm.db->aniso.u33 = gatm.mm->u_aniso_orth().mat22();
            gatm.db->aniso.u12 = gatm.mm->u_aniso_orth().mat01();
            gatm.db->aniso.u13 = gatm.mm->u_aniso_orth().mat02();
            gatm.db->aniso.u23 = gatm.mm->u_aniso_orth().mat12();
          }
          if (gatm.mm->name() != "") // set atom id/name for Gemmi
            gatm.db->name = gatm.mm->name().trim();
          if (gatm.mm->element() != "") // set atom element
            gatm.db->element = ::gemmi::Element(gatm.mm->element().c_str());
          if (gatm.mm->exists_property("Charge"))
          { // set atom charge if any
            signed char a = dynamic_cast<const Property<signed char> &>(gatm.mm->get_property("Charge")).value();
            gatm.db->charge = a;
          }
          if (gatm.mm->exists_property("AltConf"))
          { // set alt conf code
            String a = dynamic_cast<const Property<String> &>(gatm.mm->get_property("AltConf")).value();
            if (a != "")
              gatm.db->altloc = (a[0] == ' ') ? '\0' : a[0];
          }
        }
      }
    }
  }

  /*! Set spacegroup to structure.
    To set the spacegroup to the structure held in
    GEMMIFile object.
    \param sg The spacegroup to set */
  void GEMMIFile::set_spacegroup(const Spacegroup &sg)
  {
    structure_.set_spacegroup(sg);
  }

  /*! Set cell to structure.
    To set the cell to the structure held in GEMMIFile object.
    \param cell The cell to set */
  void GEMMIFile::set_cell(const Cell &cell_in)
  {
    structure_.set_cell(cell_in);
  }

  /*! Set cif output style.
    To set the cif output style, default is PreferPairs.
    Refer to GEMMI's cif::Style
    \param cifstyle 0=PreferPairs, Default, 1=Aligned, 2=Pdbx */
  void GEMMIFile::set_cif_style(CifStyle cifstyle)
  {
    switch (cifstyle)
    {
    case Aligned:
      cifstyle_ = gemmi::CGCifStyle::Aligned;
      break;
    case Pdbx:
      cifstyle_ = gemmi::CGCifStyle::Pdbx;
      break;
    case PreferPairs:
    default:
      cifstyle_ = gemmi::CGCifStyle::PreferPairs;
    }
  }

} // namespace clipper
