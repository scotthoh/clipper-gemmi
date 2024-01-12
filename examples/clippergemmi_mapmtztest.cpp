
#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-gemmi.h>
#include <clipper/clipper-minimol.h>
#include <clipper/ccp4/ccp4_map_io.h>

#include <algorithm>
#include <chrono>
#include <clipper/core/coords.h>
#include <clipper/core/hkl_datatypes.h>
#include <clipper/core/map_utils.h>
#include <clipper/core/spacegroup.h>
#include <clipper/core/symop.h>
#include <fstream>
#include <gemmi/ccp4.hpp>
#include <gemmi/grid.hpp>
#include <gemmi/gz.hpp>
#include <gemmi/math.hpp>
#include <gemmi/mtz.hpp>      // mtz
#include <gemmi/symmetry.hpp>
#include <gemmi/unitcell.hpp> //Miller
#include <gemmi/asumask.hpp>
#include <gemmi/mmread_gz.hpp>
#include <map>
#include <gemmi/it92.hpp>
#include <gemmi/c4322.hpp>
#include <gemmi/dencalc.hpp>

int main(int argc, char **argv) {
  // CCP4MTZfile::crystalinf crystal;
using namespace gemmi;
using namespace clipper;
using Miller = ::gemmi::Miller;
using namespace std::chrono;

  String pdbin = argv[1];
  auto suffixout = pdbin.split("/.").end()[-2];
  auto st = ::gemmi::read_structure_gz(pdbin);
  std::cout << "Test converting spacegroup and unit cell" << std::endl;
  std::cout << "Clipper Spacegroup : \n";
  auto clipper_spg = GEMMI::spacegroup(*st.find_spacegroup());
  clipper_spg.debug();
  std::cout << "\nClipper-Gemmi Spacegroup : \n";
  const ::gemmi::SpaceGroup * clipper_gemmi_spg = GEMMI::spacegroup(clipper_spg);
  std::cout << clipper_gemmi_spg->number << " " << clipper_gemmi_spg->operations().order() << " " << clipper_gemmi_spg->operations().sym_ops.size() << " " << clipper_gemmi_spg->hall << std::endl;
  for (auto op : clipper_gemmi_spg->operations())
    std::cout << op.triplet() << "\n";
  std::cout << "\nInitial Gemmi Spacegroup : \n";
  auto g_spg = st.find_spacegroup();
  std::cout << g_spg->number << " " << g_spg->operations().order() << " " << g_spg->operations().sym_ops.size() << " " <<  g_spg->hall << std::endl;
  for (auto op : g_spg->operations())
    std::cout << op.triplet() << "\n";
  Cell clipper_uc = GEMMI::cell(st.cell);
  std::cout <<"\nClipper UnitCell : " << clipper_uc.format() << std::endl;
  ::gemmi::UnitCell uc = GEMMI::cell(clipper_uc);
  std::cout <<"Clipper-Gemmi Unitcell : " << uc.a << "," << uc.b << "," << uc.c << "," << uc.alpha << "," << uc.beta << "," << uc.gamma << std::endl;
  std::cout <<"Initial Gemmi Unitcell : " << st.cell.a << "," << st.cell.b << "," << st.cell.c << "," << st.cell.alpha << "," << st.cell.beta << "," << st.cell.gamma << std::endl;
  

  
  Xmap<float> xmap;
  NXmap<float> nxm;
  Ccp4<float> ccp4, ccp4_p1;
  String etype = "electron";
  std::cout << "Test converting Gemmi map to Xmap/NXmap vice versa." << std::endl;
  if (etype.compare(String(argv[2])) == 0){
    std::cout << "using electon SF" << std::endl;
    clipper::ScatteringFactorsType El = clipper::ScatteringFactorsType::SF_ELECTRON;
    clipper::ScatteringFactors::selectScattteringFactorsType(El);
    DensityCalculator<C4322<double>, float> dencalc;
    //calctype dencalc;
    dencalc.d_min = 2.5;
    //dencalc.rate = 1.5; default
    //dencalc.addends = calc.addends;
    dencalc.set_refmac_compatible_blur(st.models[0]);
    std::cout << "B_min = " << get_minimum_b(st.models[0]) << ", B_add : " << dencalc.blur << std::endl;
    dencalc.set_grid_cell_and_spacegroup(st);
    dencalc.put_model_density_on_grid(st.models[0]);
    ccp4.grid = dencalc.grid;
    ccp4.update_ccp4_header(2);
    ccp4.write_ccp4_map("gemmi_dencalc_map_" + suffixout + ".ccp4");
    // clipper-gemmi import xmap
    auto start1 = high_resolution_clock::now();
    GEMMI::import_xmap(xmap, ccp4);
    auto stop1 = high_resolution_clock::now();
    std::cout << "Clipper-gemmi import Xmap time : "
              << float((duration_cast<microseconds>(stop1 - start1)).count())/1000000.0
              << " seconds" << std::endl;
    //auto st2 = st;
    std::cout << "Gemmi dencalc map stats : \n";
    std::cout << "  Dmin = " << ccp4.hstats.dmin << ", Dmax = " << ccp4.hstats.dmax << ", Dmean = " << ccp4.hstats.dmean << ", RMS = " << ccp4.hstats.rms << std::endl;
    std::cout << "  Grid shape = " << dencalc.grid.nu << "," <<  dencalc.grid.nv << ","  << dencalc.grid.nw << std::endl;
    auto st2 = st;
    st2.spacegroup_hm = "P 1";
    std::cout << "Spacegroup : " << st.spacegroup_hm << " change to : " << st2.spacegroup_hm << std::endl;
    //dencalc.set_refmac_compatible_blur(st2.models[0]);
    std::cout << "B_min = " << get_minimum_b(st2.models[0]) << ", B_add : " << dencalc.blur << std::endl;
    dencalc.set_grid_cell_and_spacegroup(st2);
    dencalc.put_model_density_on_grid(st2.models[0]);
    Map_stats stats(xmap);
    std::cout << "Clipper imported Xmap stats : \n";
    std::cout << "  Dmin = " << stats.min() << ", Dmax = "<< stats.max() << ", Dmean = " << stats.mean() << ", RMS = " << stats.std_dev() << std::endl;
    std::cout << "  Grid shape = " << xmap.grid_sampling().nu() << ", " << xmap.grid_sampling().nv() << ", " << xmap.grid_sampling().nw() << std::endl;
    ccp4_p1.grid = dencalc.grid;
    ccp4_p1.update_ccp4_header(2);
    ccp4_p1.write_ccp4_map("gemmi_dencalc_map_p1_" + suffixout + ".ccp4");
    // clipper-gemmi import nxmap
    start1 = high_resolution_clock::now();
    GEMMI::import_nxmap(nxm, ccp4_p1);
    stop1 = high_resolution_clock::now();
    std::cout << "Clipper-gemmi import NXmap time : "
              << float((duration_cast<microseconds>(stop1 - start1)).count())/1000000.0
              << " seconds" << std::endl;
    std::cout << "Gemmi P1 map stats : \n";
    std::cout << "  Dmin = " << ccp4_p1.hstats.dmin << ", Dmax = "<< ccp4_p1.hstats.dmax << ", Dmean = " << ccp4_p1.hstats.dmean << ", RMS = " << ccp4_p1.hstats.rms << std::endl;
    std::cout << "  Grid shape = " << ccp4_p1.grid.nu << "," <<  ccp4_p1.grid.nv << ","  << ccp4_p1.grid.nw << std::endl;
    Map_stats stats_p1(nxm);
    std::cout << "Clipper imported NXmap map stats : \n";
    std::cout << "  Dmin = " << stats_p1.min() << ", Dmax = "<< stats_p1.max() << ", Dmean = " << stats_p1.mean() << ", RMS = " << stats_p1.std_dev() << std::endl;
    std::cout << "  Grid shape = " << nxm.grid().nu() << "," <<  nxm.grid().nv() << ","  << nxm.grid().nw() << std::endl;
  
  }
  else {
  //calctype dencalc;
  std::cout << "using X-ray SF" << std::endl;
  DensityCalculator<IT92<double>, float> dencalc;
  dencalc.d_min = 2.5;
  //dencalc.rate = 1.5; default
  //dencalc.addends = calc.addends;
  dencalc.set_refmac_compatible_blur(st.models[0]);
  std::cout << "B_min = " << get_minimum_b(st.models[0]) << ", B_add : " << dencalc.blur << std::endl;
  dencalc.set_grid_cell_and_spacegroup(st);
  dencalc.put_model_density_on_grid(st.models[0]);
  ccp4.grid = dencalc.grid;
  ccp4.update_ccp4_header(2);
  ccp4.write_ccp4_map("gemmi_dencalc_map_" + suffixout + ".ccp4");
  //auto st2 = st;
  // clipper-gemmi import xmap
  auto start1 = high_resolution_clock::now();
  GEMMI::import_xmap(xmap, ccp4);
  auto stop1 = high_resolution_clock::now();
  std::cout << "Clipper-gemmi import Xmap time : "
            << float((duration_cast<microseconds>(stop1 - start1)).count())/1000000.0
            << " seconds" << std::endl;
  std::cout << "Gemmi dencalc map stats : \n";
  std::cout << "  Dmin = " << ccp4.hstats.dmin << ", Dmax = " << ccp4.hstats.dmax << ", Dmean = " << ccp4.hstats.dmean << ", RMS = " << ccp4.hstats.rms << std::endl;
  std::cout << "  Grid shape = " << ccp4.grid.nu << "," <<  ccp4.grid.nv << ","  << ccp4.grid.nw << std::endl;
  Map_stats stats(xmap);
  std::cout << "Clipper imported Xmap stats : \n";
  std::cout << "  Dmin = " << stats.min() << ", Dmax = "<< stats.max() << ", Dmean = " << stats.mean() << ", RMS = " << stats.std_dev() << std::endl;
  std::cout << "  Grid shape = " << xmap.grid_sampling().nu() << "," <<  xmap.grid_sampling().nv() << ","  << xmap.grid_sampling().nw() << std::endl;
  auto st2 = st;
  st2.spacegroup_hm = "P 1";
  std::cout << "Spacegroup : " << st.spacegroup_hm << " ; now: " << st2.spacegroup_hm << std::endl;
  //dencalc.set_refmac_compatible_blur(st2.models[0]);
  std::cout << "B_min = " << get_minimum_b(st2.models[0]) << ", B_add : " << dencalc.blur << std::endl;
  dencalc.set_grid_cell_and_spacegroup(st2);
  dencalc.put_model_density_on_grid(st2.models[0]);
  ccp4_p1.grid = dencalc.grid;
  ccp4_p1.update_ccp4_header(2);
  ccp4_p1.write_ccp4_map("gemmi_dencalc_map_p1_" + suffixout + ".ccp4");
  // clipper-gemmi import nxmap
  start1 = high_resolution_clock::now();
  GEMMI::import_nxmap(nxm, ccp4_p1);
  stop1 = high_resolution_clock::now();
  std::cout << "Clipper-gemmi import NXmap time : "
            << float((duration_cast<microseconds>(stop1 - start1)).count())/1000000.0
            << " seconds" << std::endl;
  std::cout << "Gemmi P1 map stats : \n";
  std::cout << "  Dmin = " << ccp4_p1.hstats.dmin << ", Dmax = "<< ccp4_p1.hstats.dmax << ", Dmean = " << ccp4_p1.hstats.dmean << ", RMS = " << ccp4_p1.hstats.rms << std::endl;
  std::cout << "  Grid shape = " << ccp4_p1.grid.nu << "," <<  ccp4_p1.grid.nv << ","  << ccp4_p1.grid.nw << std::endl;
  Map_stats stats_p1(nxm);
  std::cout << "Clipper imported NXmap map stats : \n";
  std::cout << "  Dmin = " << stats_p1.min() << ", Dmax = "<< stats_p1.max() << ", Dmean = " << stats_p1.mean() << ", RMS = " << stats_p1.std_dev() << std::endl;
  std::cout << "  Grid shape = " << nxm.grid().nu() << "," <<  nxm.grid().nv() << ","  << nxm.grid().nw() << std::endl;
  }
  std::cout << "Read Gemmi calculated map into clipper, and export back to Gemmi map object" << std::endl;
  // Read in gemmi calculated map with clipper
  CCP4MAPfile mapio;
  mapio.open_read("gemmi_dencalc_map_" + suffixout + ".ccp4");
  clipper::Xmap<float> xmap1;
  auto start2 = high_resolution_clock::now();
  mapio.import_xmap(xmap1);
  mapio.close_read();
  auto stop2 = high_resolution_clock::now();
  std::cout << "  Gemmi map grid nu,nv,nw : " << ccp4.grid.nu << ", " << ccp4.grid.nv << ", " << ccp4.grid.nw << std::endl;
  std::cout << "  Clipper Xmap grid sampling : " << xmap1.grid_sampling().format() << std::endl;
  std::cout << "  Clipper Xmap asu min : " << xmap1.grid_asu().min().format() << std::endl;
  std::cout << "  Clipper Xmap asu max : " << xmap1.grid_asu().max().format() << std::endl;
  std::cout << "Import clipper Xmap time : "
            << float((duration_cast<microseconds>(stop2 - start2)).count())/1000000.0
            << " seconds" << std::endl;
  Map_stats stats(xmap1);
  std::cout << "Clipper Xmap stats : \n";
  std::cout << " Dmin = " << stats.min() << ", Dmax = "<< stats.max() << ", Dmean = " << stats.mean() << ", RMS = " << stats.std_dev() << std::endl;
  Ccp4<float> gmap;
  // export xmap to clipper map object
  start2 = high_resolution_clock::now();
  GEMMI::export_xmap(xmap1, gmap);
  stop2 = high_resolution_clock::now();
  std::cout << "Clipper-gemmi export Xmap time : "
              << float((duration_cast<microseconds>(stop2 - start2)).count())/1000000.0
              << " seconds" << std::endl;
  std::cout << "Clipper exported Xmap-Gemmi map stats :-\n";
  std::cout << "  Dmin = " << gmap.hstats.dmin << ", Dmax = "<< gmap.hstats.dmax << ", Dmean = " << gmap.hstats.dmean << ", RMS = " << gmap.hstats.rms << std::endl;
  std::cout << "  Grid shape = " << gmap.grid.nu << "," <<  gmap.grid.nv << ","  << gmap.grid.nw << std::endl;
  gmap.write_ccp4_map("export_clippergemmi_" + suffixout + ".map");
  // to check if map exported correctly
  // calculate correlation between clipper-gemmi exported map with previously gemmi calculated map from model
  auto corr = calculate_correlation(ccp4.grid, gmap.grid);
  std::cout << "  Correlation between Gemmi calculated map and clipper-gemmi exported map"  << std::endl;
  std::cout << "  Expecting correlation coefficient of 1.0" << std::endl;
  std::cout << "  Correlation Coefficient : " << corr.coefficient();
  std::cout << " ; Mean ratio : " << corr.mean_ratio();
  std::cout << " ; Correlation n : " << corr.n  << std::endl;
  
  // read in gemmi calculated p1 map
  mapio.open_read("gemmi_dencalc_map_p1_" + suffixout + ".ccp4");
  clipper::NXmap<float> nxmap1;
  start2 = high_resolution_clock::now();
  mapio.import_nxmap(nxmap1);
  mapio.close_read();
  stop2 = high_resolution_clock::now();
  std::cout << "  Gemmi map grid nu,nv,nw : " << ccp4_p1.grid.nu << ", " << ccp4_p1.grid.nv << ", " << ccp4_p1.grid.nw << std::endl;
  std::cout << "  Clipper NXmap grid sampling : " << nxmap1.grid().format() << std::endl;
  std::cout << "Import clipper NXmap time : "
            << float((duration_cast<microseconds>(stop2 - start2)).count())/1000000.0
            << " seconds" << std::endl;
  Map_stats stats_p1(nxmap1);
  std::cout << "  Clipper NXmap stats : \n";
  std::cout << "  Dmin = " << stats_p1.min() << ", Dmax = "<< stats_p1.max() << ", Dmean = " << stats_p1.mean() << ", RMS = " << stats_p1.std_dev() << std::endl;
  Ccp4<float> gmap_p1;
  // export clipper nxmap to gemmi map object
  start2 = high_resolution_clock::now();
  GEMMI::export_nxmap(nxmap1, gmap_p1, mapio.cell());
  stop2 = high_resolution_clock::now();
  std::cout << "Clipper-gemmi export NXmap time : "
              << float((duration_cast<microseconds>(stop2 - start2)).count())/1000000.0
              << " seconds" << std::endl;
  std::cout << "Clipper exported NXmap-Gemmi map stats :-\n";
  std::cout << "  Dmin = " << gmap_p1.hstats.dmin << ", Dmax = "<< gmap_p1.hstats.dmax << ", Dmean = " << gmap_p1.hstats.dmean << ", RMS = " << gmap_p1.hstats.rms << std::endl;
  std::cout << "  Grid shape = " << gmap_p1.grid.nu << "," <<  gmap_p1.grid.nv << ","  << gmap_p1.grid.nw << std::endl;
  gmap_p1.write_ccp4_map("export_clippergemmi_p1" + suffixout + ".map");
  // to check if map exported correctly
  // calculate correlation between clipper-gemmi exported map with previously gemmi calculated map from model
  auto corr_p1 = calculate_correlation(ccp4_p1.grid, gmap_p1.grid);
  std::cout << "  Correlation between Gemmi calculated map and clipper-gemmi exported map"  << std::endl;
  std::cout << "  Expecting correlation coefficient of 1.0" << std::endl;
  std::cout << "  Correlation Coefficient : " << corr_p1.coefficient();
  std::cout << " ; Mean ratio : " << corr_p1.mean_ratio();
  std::cout << " ; Correlation n : " << corr_p1.n  << std::endl;
  
  std::cout << "Convert Gemmi Transform to Clipper RTop_orth\n" << std::endl;
  RTop_orth rtop = GEMMI::transform(st.origx);
  std::cout << " "  << rtop.format() << std::endl;
  std::cout << "Convert Clipper RTop_orth back to Gemmi Transform\n" << std::endl;
  Transform tr = GEMMI::transform(rtop);
  std::cout << "Rot Mat33 : \n";
  for (int i = 0 ; i < 3; i++){
    std::cout << " ";
    for (int j = 0; j < 3; j++){
      std::cout << tr.mat[i][j] << ", ";
    }
    std::cout << std::endl;
  }
  std::cout << "Trans Vec3 : \n";
  std::cout << tr.vec.x << ", " << tr.vec.y << ", " << tr.vec.z << std::endl;
  std::cout << "\nTest converting reflection data from Gemmi to Clipper and vice versa." << std::endl;
  
  HKL hkl(1,2,3);
  Miller miller{2,3,4};

  Miller hkl2m = GEMMI::Hkl(hkl);
  std::cout << " HKL to Miller : " << hkl2m[0] << "," << hkl2m[1] << "," << hkl2m[2] << std::endl;
  std::cout << " Expecting : " << hkl.format() << std::endl;
  HKL m2hkl = GEMMI::Hkl(miller);
  std::cout << " Miller to HKL : " << m2hkl.format() << std::endl;
  std::cout << " Expecting : " << miller[0] << "," << miller[1] << "," << miller[2] <<std::endl;

  Mtz mtzobj;
  mtzobj.read_file(argv[3]);
  HKL_info hklinfo=GEMMI::as_HKL_info(mtzobj.get_cell(), *mtzobj.spacegroup, mtzobj.resolution_high());
  
  std::vector<HKL> hkllist;
  std::vector<Miller> millervec;
  hkllist.reserve(mtzobj.nreflections);
  double max_d_star_sq=0.0;
  for (std::size_t i = 0; i < mtzobj.data.size(); i += mtzobj.columns.size())
  {
    Miller miller = mtzobj.get_hkl(i);
    millervec.push_back(miller);
    max_d_star_sq = std::max(max_d_star_sq, mtzobj.cell.calculate_1_d2(miller));
    hkllist.push_back(HKL(miller[0], miller[1],miller[2]));
  }
  std::cout << max_d_star_sq << ", " << mtzobj.max_1_d2 << std::endl;
  hklinfo.add_hkl_list(hkllist);
  HKL_data<data64::F_sigF> fsigf(hklinfo);
  std::cout << "Before import hkl data : " << std::endl;
  //fsigf.debug();
  data64::F_sigF datum;
  bool exist_data = fsigf.get_data(hkllist[0], datum);
  if (exist_data)
    std::cout << " " << hkllist[0].format() << "=" << datum.f() << ", " << datum.sigf() << std::endl;
  else
    std::cout << "no data" << std::endl;
  std::cout << "After import hkl data : " << std::endl;
  // general method can be used for other HKL datatypes i.e. I_sigI, F_phi, Phi_fom etc
  GEMMI::import_hkl_data(fsigf, mtzobj, String("/*/*/[FP,SIGFP]"));
  exist_data = fsigf.get_data(hkllist[0], datum);
  if (exist_data)
    std::cout << " " << hkllist[0].format() << "=" << datum.f() << ", " << datum.sigf() << std::endl;
  else
   std::cout << "no data" << std::endl;

  HKL_data<data64::F_phi> fphi(hklinfo);
  data64::F_phi fphidatum;
  exist_data = fphi.get_data(hkllist[0], fphidatum);
  if(exist_data)
    std::cout << " " << hkllist[0].format() << "=" << fphidatum.f() << ", " << fphidatum.phi() << std::endl;
  else std::cout << "no data" << std::endl;

  GEMMI::import_hkl_data(fphi, mtzobj, String("/*/*/[FWT,PHWT]"));
  exist_data = fphi.get_data(hkllist[0], fphidatum);
  if(exist_data)
    std::cout << " " << hkllist[0].format() << "=" << fphidatum.f() << ", " << fphidatum.phi() << std::endl;
  else std::cout << "no data" << std::endl;
  
  // specific
  //HKL_data<data64::F_phi> f_phidata;
  //for (size_t i = 0; i < mtzobj.data.size(); i+=mtzobj.columns.size())
  //{
//
  //}
  //f_phidata = GEMMI::as_HKL_data_fphi(hklinfo, millervec, const std::vector<double> &data_f, const std::vector<double> &data_phi)

  auto hklinfo2 = GEMMI::as_HKL_info(mtzobj);
  std::cout << "Mtz nreflections : " << mtzobj.nreflections << " ; ";
  hklinfo2.debug();
  HKL_data<data64::F_phi> fphi2(hklinfo2);
  GEMMI::import_hkl_data(fphi2, mtzobj, String("/*/*/[FWT,PHWT]"));
  exist_data = fphi2.get_data(hkllist[0], fphidatum);
  if(exist_data)
    std::cout << " " << hkllist[0].format() << "=" << fphidatum.f() << ", " << fphidatum.phi() << std::endl;
  else std::cout << "no data" << std::endl;
  
  }