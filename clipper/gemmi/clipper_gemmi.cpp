/* clipper_gemmi.cpp: Gemmi wrapper */

#include "clipper_gemmi.h"
#include <cstddef>
#include <gemmi/grid.hpp>
#include <map>
#include <ostream>
#include <set>

namespace clipper
{

/* Get dataset name and a list of column names from clipper path name. */
const std::vector<String> extract_column_names(const String &assign, const int &f_size) {
  std::vector<String> col_names(f_size, "MNF");
  // interpret list name in terms of columns
  if (assign.find("[") == String::npos) {
    // name is a single path
    col_names[0] = assign.split("/")[2];
  } else {
    // name is a list of mtz columns: extract column names from list
    String pref = assign.substr(0, assign.find("["));
    String post = assign.substr(assign.find("[") + 1, assign.find("]") - assign.find("[") - 1);
    std::vector<String> list = post.split(", ");
    for (int i = 0; i < list.size(); i++)
      col_names[i] = list[i];
  }
  return col_names;
}

const ::gemmi::SpaceGroup * GEMMI::spacegroup(const Spacegroup &spgr) {
  gemmi::GroupOps symops = gemmi::symops_from_hall(spgr.symbol_hall().c_str());
  return gemmi::find_spacegroup_by_ops(symops);
}

Spacegroup GEMMI::spacegroup(const gemmi::SpaceGroup &spgr) {
  String ops;
  for (auto op : spgr.operations())
    ops += op.triplet() + ";";
  return Spacegroup(Spgr_descr(ops, Spgr_descr::Symops));
}

::gemmi::Transform GEMMI::transform(const RTop_orth &rtop) {
  gemmi::Transform grtop; //identity mat33, zero vec3
  if (!rtop.rot().is_null()) {
    for (int i = 0; i < 3; i++) {
      grtop.mat[i][0] = rtop.rot()(i, 0);
      grtop.mat[i][1] = rtop.rot()(i, 1);
      grtop.mat[i][2] = rtop.rot()(i, 2);
    }
  }
  if (!rtop.trn().is_null()) {
    for (int i = 0; i < 3; i++) {
      grtop.vec.at(i) = rtop.trn()[i];
    }
  }
  return grtop;
}

RTop_orth GEMMI::transform(const gemmi::Transform &grtop) {
  Mat33<> m(grtop.mat[0][0], grtop.mat[0][1], grtop.mat[0][2], grtop.mat[1][0], grtop.mat[1][1], grtop.mat[1][2],
            grtop.mat[2][0], grtop.mat[2][1], grtop.mat[2][2]);
  Vec3<> v(grtop.vec.at(0), grtop.vec.at(1), grtop.vec.at(2));
  return RTop_orth(m, v);
}

Miller GEMMI::Hkl(const HKL &hkl) { return Miller{hkl.h(), hkl.k(), hkl.l()}; }

HKL GEMMI::Hkl(const Miller &hkl) { return HKL(hkl[0], hkl[1], hkl[2]); }

gemmi::UnitCell GEMMI::cell(const Cell &cell) {
  return gemmi::UnitCell(cell.a(), cell.b(), cell.c(), cell.alpha_deg(), cell.beta_deg(), cell.gamma_deg());
}

Cell GEMMI::cell(const gemmi::UnitCell &cell) {
  return Cell(Cell_descr(cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma));
}

void GEMMI::import_mtz_crystal(std::vector<CCP4MTZfile::crystalinf> &crystals, const gemmi::Mtz &mtzobj) {
  crystals.clear();
  CCP4MTZfile::crystalinf newxtl;
  CCP4MTZfile::datasetinf newset;
  CCP4MTZfile::datacolinf newcol;
  // std::set<String> crystname;
  std::map<String, int> crystname;
  int count = 0;
  for (int s = 0; s < mtzobj.datasets.size(); s++) {
    const gemmi::Mtz::Dataset *dataset = &mtzobj.datasets[s];
    auto result = crystname.insert({String(dataset->crystal_name), count});
    if (result.second) {
      newxtl.crystal = MTZcrystal(String(dataset->crystal_name), String(dataset->project_name), cell(dataset->cell));
      // if (newxtl.crystal.crystal_name != dataset->crystal_name)
      crystals.push_back(newxtl);
      count++;
    }

    newset.dataset = MTZdataset(dataset->dataset_name, dataset->wavelength);
    int cryst_idx = crystname.find(dataset->crystal_name)->second;
    crystals[cryst_idx].datasets.push_back(newset);
    for (int c = 0; c < mtzobj.columns.size(); c++) {
      const gemmi::Mtz::Column *col = &mtzobj.columns[c];
      newcol.label = col->label;
      newcol.type = String(&col->type, 1);
      newcol.source = col->source;
      newcol.grpname = "";
      newcol.grptype = "";
      newcol.grpposn = -1;
      // gemmi not using Column group
      crystals[cryst_idx].datasets.back().columns.push_back(newcol);
    }
  }
}

/*! Import data from gemmi's Mtz object into HKL_data object,
  given the column paths.

  An MTZ column type must be present in the MTZ_type_registry for the
  HKL_data type element name concerned.

  \param cdata The HKL_data object into which data is to be imported.
  \param mtzobj Gemmi Mtz class where the data is store.
  \param mtzpath The MTZ column names, as a path. See \ref MTZpaths for details. */
//template <class T>
void GEMMI::import_hkl_data(HKL_data_base &cdata, const gemmi::Mtz &mtzobj, const String mtzpath) {
  String colpath = mtzpath;
  // add exported data columns to local list
  int ncols = cdata.data_size();
  size_t cols[ncols];
  std::vector<String> col_names = extract_column_names(colpath, ncols);
  std::vector<String> dat_names = cdata.data_names().split(" ");
  std::vector<CCP4MTZfile::hkldatacol> newcols(ncols);
  String dataset_name = colpath.split("/")[1];
  // get the column indices

  for (int c = 0; c < ncols; c++) {
    if (col_names[c] != "MNF" && col_names[c] != "NAN" && col_names[c] != "mnf" && col_names[c] != "nan") {
      if (dataset_name != "*")
        cols[c] = mtzobj.column_with_label(col_names[c].c_str(), mtzobj.dataset_with_name(dataset_name.c_str()))->idx;
      else
        cols[c] = mtzobj.column_with_label(col_names[c].c_str())->idx;
    }
  }
  // loop through nreflections,data
  xtype values[ncols];
  for (size_t i = 0; i < mtzobj.data.size(); i += mtzobj.columns.size()) {
    // import data
    for (int c = 0; c < ncols; c++) {
      // set to Nan unless readable and present
      Util::set_null(values[c]);
      if (!Util::is_nan(xtype(mtzobj.data[i + cols[c]])))
        values[c] = xtype(mtzobj.data[i + cols[c]]);
    }
    cdata.data_import(Hkl(mtzobj.get_hkl(i)), values);
  }
}

/*! Return initialised HKL_info using data from Gemmi::Mtz which
  then can be used to construct HKL_data. 
  \param mtzobj Gemmi Mtz class where the data is store.
  \param tolerance Tolerance limit. */
HKL_info GEMMI::as_HKL_info(const gemmi::Mtz &mtzobj, double tolerance)
{
  if(mtzobj.max_1_d2 == 0)
    Message::message(Message_fatal("GEMMI: max_d_star_sq is zero."));
  std::vector<HKL> hkllist;
  hkllist.reserve(mtzobj.nreflections);
  
  for(size_t i = 0; i < mtzobj.data.size(); i += mtzobj.columns.size())
    hkllist.push_back(Hkl(mtzobj.get_hkl(i)));
  
  HKL_info result = as_HKL_info(mtzobj.cell, *mtzobj.spacegroup, mtzobj.resolution_high(), tolerance);
  result.add_hkl_list(hkllist);
  return result;
}

/*! Return initialised HKL_info using Gemmi unit cell and spacegroup which
  then can be used to construct HKL_data. 
  \param unit_cell The unit cell
  \param sg The spacegroup
  \param d_min The resolution limit
  \param tolerance Tolerance limit. */
HKL_info GEMMI::as_HKL_info(const gemmi::UnitCell &unit_cell, const gemmi::SpaceGroup &sg, double d_min,
                            double tolerance) {
  return HKL_info(spacegroup(sg), cell(unit_cell), Resolution(d_min - tolerance));
}

/*! Return initialised HKL_info using Gemmi unit cell and spacegroup which
  then can be used to construct HKL_data. 
  \param unit_cell The unit cell
  \param sg The spacegroup
  \param d_min The resolution limit
  \param miller Vector of Miller indices
  \param tolerance Tolerance limit. */
HKL_info GEMMI::as_HKL_info(const gemmi::UnitCell &unit_cell, const gemmi::SpaceGroup &space_group,
                            const std::vector<gemmi::Miller> &miller_indices, double tolerance) {
  std::vector<HKL> hkl_list;
  hkl_list.reserve(miller_indices.size());
  double max_d_star_sq = 0;
  // hkl_list.reserve();
  for (std::size_t i = 0; i < miller_indices.size(); i++) {
    max_d_star_sq = std::max(max_d_star_sq, unit_cell.calculate_1_d2(miller_indices[i]));
    hkl_list.push_back(Hkl(miller_indices[i]));
  }
  if (max_d_star_sq == 0)
    Message::message(Message_fatal("GEMMI: max_d_star_sq is zero."));
  HKL_info result = as_HKL_info(unit_cell, space_group, 1 / std::sqrt(max_d_star_sq), tolerance);
  result.add_hkl_list(hkl_list);
  return result;
}

HKL_data<data64::F_sigF> GEMMI::as_HKL_data(HKL_info &hkl_info, const std::vector<gemmi::Miller> &miller_indices,
                                            const std::vector<double> &data, const std::vector<double> &sigmas) {
  if (data.size() != miller_indices.size())
    Message::message(Message_fatal("GEMMI: Vectors for data and miller indices have different lengths."));
  if (sigmas.size() != miller_indices.size())
    Message::message(Message_fatal("GEMMI: Vectors for sigmas and miller indices have different lengths."));
  HKL_data<data64::F_sigF> hkl_data(hkl_info);
  for (std::size_t i = 0; i < miller_indices.size(); i++) {
    data64::F_sigF datum;
    datum.f() = data[i];
    datum.sigf() = sigmas[i];
    if (!hkl_data.set_data(Hkl(miller_indices[i]), datum))
      Message::message(Message_fatal("GEMMI: Unable to set data for " + Hkl(miller_indices[i]).format()));
  }
  return hkl_data;
}

HKL_data<data64::F_phi> GEMMI::as_HKL_data(HKL_info &hkl_info, const gemmi::AsuData<std::complex<double>> &data) {
  HKL_data<data64::F_phi> hkl_data(hkl_info);
  for (std::size_t i = 0; i < data.size(); i++) {
    if (!hkl_data.set_data(Hkl(data.get_hkl(i)), data64::F_phi(data.v[i].value)))
      Message::message(Message_fatal("GEMMI: Unable to set data for " + Hkl(data.get_hkl(i)).format()));
  }
  return hkl_data;
}

HKL_data<data64::F_phi> GEMMI::as_HKL_data_fphi(HKL_info &hkl_info, const std::vector<gemmi::Miller> &miller_indices,
                                                const std::vector<double> &data_f,
                                                const std::vector<double> &data_phi) {
  if ((data_f.size() != miller_indices.size()) || (data_phi.size() != miller_indices.size()))
    Message::message(Message_fatal("GEMMI: Vectors for data and miller indices have different lengths."));
  HKL_data<data64::F_phi> hkl_data(hkl_info);
  for (std::size_t i = 0; i < data_f.size(); i++) {
    if (!hkl_data.set_data(Hkl(miller_indices[i]), data64::F_phi(data_f[i], data_phi[i])))
      Message::message(Message_fatal("GEMMI: Unable to set data for " + Hkl(miller_indices[i]).format()));
  }
  return hkl_data;
}

std::vector<std::complex<double>> GEMMI::extract_complex(const HKL_data<data64::F_phi> &hkl_data,
                                                         const std::vector<Miller> &miller_indices) {
  std::vector<std::complex<double>> result(miller_indices.size());
  data64::F_phi datum;
  for (std::size_t i = 0; i < miller_indices.size(); i++) {
    if (!hkl_data.get_data(Hkl(miller_indices[i]), datum))
      Message::message(Message_fatal("CLIPPER: Unable to get data for " + Hkl(miller_indices[i]).format()));
    if (datum.missing())
      Message::message(Message_fatal("CLIPPER: Missing F_phi data for " + Hkl(miller_indices[i]).format()));
    result.push_back(datum);
  }
  return result;
}

/*! Converts a std::vector<const gemmi::Mtz::Column*> returned
  from gemmi method Mtz::columns_with_type('A') to
  clipper's HKL_data<data64::ABCD>.
  \param hkl_info HKL_info object
  \param miller_indices a vector of gemmi::Miller
  \param data a vector of  Mtz::Column* */
HKL_data<data64::ABCD> GEMMI::as_HKL_data(HKL_info &hkl_info, const std::vector<gemmi::Miller> &miller_indices,
                                          const std::vector<const gemmi::Mtz::Column *> &data) {
  if (data.size() != 4)
    Message::message(Message_fatal("GEMMI: Size for vector<const gemmi::Mtz::Column*> data not 4."));
  if (!std::all_of(data.cbegin(), data.cend(),
                   [&](const gemmi::Mtz::Column *c) { return c->size() == miller_indices.size(); }))
    Message::message(Message_fatal("GEMMI: Vectors for data and miller indices have different lengths."));

  HKL_data<data64::ABCD> hkl_data(hkl_info);
  data64::ABCD abcd;
  for (std::size_t i = 0; i < miller_indices.size(); i++) {
    std::array<double, 4> datum = {data[0]->at(i), data[1]->at(i), data[2]->at(i), data[3]->at(i)};
    abcd.data_import(datum.begin());
    if (!hkl_data.set_data(Hkl(miller_indices[i]), abcd))
      Message::message(Message_fatal("GEMMI: Unable to set data for " + Hkl(miller_indices[i]).format()));
  }
  return hkl_data;
}

std::vector<std::vector<double>> GEMMI::extract_hendrickson_lattman(const HKL_data<data64::ABCD> &hkl_data,
                                                                    const std::vector<gemmi::Miller> &miller_indices) {
  std::vector<std::vector<double>> result(miller_indices.size());
  data64::ABCD datum;
  for (std::size_t i = 0; i < miller_indices.size(); i++) {
    if (!hkl_data.get_data(Hkl(miller_indices[i]), datum))
      Message::message(Message_fatal("CLIPPER: Failed to retrieve Hendrickson-Lattman coefficients for " +
                                     Hkl(miller_indices[i]).format()));
    if (datum.missing())
      Message::message(
          Message_fatal("CLIPPER: Missing Hendrickson-Lattman coefficients for " + Hkl(miller_indices[i]).format()));
    std::vector<double> hl_coef = {datum.a(), datum.b(), datum.c(), datum.d()};
    result.push_back(hl_coef);
  }
  return result;
}

std::vector<double> GEMMI::extract_centroid_phases(const HKL_data<data64::Phi_fom> &hkl_data,
                                                   const std::vector<gemmi::Miller> &miller_indices) {
  std::vector<double> result(miller_indices.size());
  data64::Phi_fom datum;
  for (std::size_t i = 0; i < miller_indices.size(); i++) {
    if (!hkl_data.get_data(Hkl(miller_indices[i]), datum))
      Message::message(Message_fatal("CLIPPER: Failed to retrieve Phi for " + Hkl(miller_indices[i]).format()));
    if (datum.missing())
      Message::message(Message_fatal("CLIPPER: Missing Phi for " + Hkl(miller_indices[i]).format()));
    result.push_back(datum.phi());
  }
  return result;
}

std::vector<double> GEMMI::extract_figure_of_merit(const HKL_data<data64::Phi_fom> &hkl_data,
                                                   const std::vector<gemmi::Miller> &miller_indices) {
  std::vector<double> result(miller_indices.size());
  data64::Phi_fom datum;
  for (std::size_t i = 0; i < miller_indices.size(); i++) {
    if (!hkl_data.get_data(Hkl(miller_indices[i]), datum))
      Message::message(Message_fatal("CLIPPER: Failed to retrieve FOM for " + Hkl(miller_indices[i]).format()));
    if (datum.missing())
      Message::message(Message_fatal("CLIPPER: Missing FOM for " + Hkl(miller_indices[i]).format()));
    result.push_back(datum.fom());
  }
  return result;
}

} // namespace clipper

