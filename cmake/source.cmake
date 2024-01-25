set( clipper-core_sources 
${CMAKE_SOURCE_DIR}/clipper/core/atomsf.cpp
${CMAKE_SOURCE_DIR}/clipper/core/coords.cpp
${CMAKE_SOURCE_DIR}/clipper/core/nxmap_operator.cpp
${CMAKE_SOURCE_DIR}/clipper/core/cell.cpp
${CMAKE_SOURCE_DIR}/clipper/core/derivs.cpp
${CMAKE_SOURCE_DIR}/clipper/core/ramachandran.cpp
${CMAKE_SOURCE_DIR}/clipper/core/clipper_instance.cpp
${CMAKE_SOURCE_DIR}/clipper/core/fftmap.cpp
${CMAKE_SOURCE_DIR}/clipper/core/resol_basisfn.cpp
${CMAKE_SOURCE_DIR}/clipper/core/clipper_memory.cpp
${CMAKE_SOURCE_DIR}/clipper/core/fftmap_sparse.cpp
${CMAKE_SOURCE_DIR}/clipper/core/resol_fn.cpp
${CMAKE_SOURCE_DIR}/clipper/core/clipper_message.cpp
${CMAKE_SOURCE_DIR}/clipper/core/hkl_compute.cpp
${CMAKE_SOURCE_DIR}/clipper/core/resol_targetfn.cpp
${CMAKE_SOURCE_DIR}/clipper/core/clipper_stats.cpp
${CMAKE_SOURCE_DIR}/clipper/core/hkl_data.cpp
${CMAKE_SOURCE_DIR}/clipper/core/rotation.cpp
${CMAKE_SOURCE_DIR}/clipper/core/clipper_test.cpp
${CMAKE_SOURCE_DIR}/clipper/core/hkl_datatypes.cpp
${CMAKE_SOURCE_DIR}/clipper/core/spacegroup.cpp
${CMAKE_SOURCE_DIR}/clipper/core/clipper_types.cpp
${CMAKE_SOURCE_DIR}/clipper/core/hkl_info.cpp
${CMAKE_SOURCE_DIR}/clipper/core/spacegroup_data.cpp
${CMAKE_SOURCE_DIR}/clipper/core/clipper_util.cpp
${CMAKE_SOURCE_DIR}/clipper/core/hkl_lookup.cpp
${CMAKE_SOURCE_DIR}/clipper/core/symop.cpp
${CMAKE_SOURCE_DIR}/clipper/core/container.cpp
${CMAKE_SOURCE_DIR}/clipper/core/hkl_operators.cpp
${CMAKE_SOURCE_DIR}/clipper/core/container_hkl.cpp
${CMAKE_SOURCE_DIR}/clipper/core/map_interp.cpp
${CMAKE_SOURCE_DIR}/clipper/core/container_map.cpp
${CMAKE_SOURCE_DIR}/clipper/core/map_utils.cpp
${CMAKE_SOURCE_DIR}/clipper/core/xmap.cpp
${CMAKE_SOURCE_DIR}/clipper/core/container_types.cpp
${CMAKE_SOURCE_DIR}/clipper/core/nxmap.cpp
${CMAKE_SOURCE_DIR}/clipper/core/clipper_thread.cpp
${CMAKE_SOURCE_DIR}/clipper/core/test_core.cpp
${CMAKE_SOURCE_DIR}/clipper/core/test_data.cpp
)
set(clipper-core_headers 
${CMAKE_SOURCE_DIR}/clipper/core/atomsf.h
${CMAKE_SOURCE_DIR}/clipper/core/coords.h
${CMAKE_SOURCE_DIR}/clipper/core/nxmap_operator.h
${CMAKE_SOURCE_DIR}/clipper/core/cell.h
${CMAKE_SOURCE_DIR}/clipper/core/derivs.h
${CMAKE_SOURCE_DIR}/clipper/core/ramachandran.h
${CMAKE_SOURCE_DIR}/clipper/core/clipper_instance.h
${CMAKE_SOURCE_DIR}/clipper/core/fftmap.h
${CMAKE_SOURCE_DIR}/clipper/core/resol_basisfn.h
${CMAKE_SOURCE_DIR}/clipper/core/clipper_memory.h
${CMAKE_SOURCE_DIR}/clipper/core/fftmap_sparse.h
${CMAKE_SOURCE_DIR}/clipper/core/clipper_precision.h
${CMAKE_SOURCE_DIR}/clipper/core/resol_fn.h
${CMAKE_SOURCE_DIR}/clipper/core/clipper_message.h
${CMAKE_SOURCE_DIR}/clipper/core/hkl_compute.h
${CMAKE_SOURCE_DIR}/clipper/core/resol_targetfn.h
${CMAKE_SOURCE_DIR}/clipper/core/clipper_stats.h
${CMAKE_SOURCE_DIR}/clipper/core/hkl_data.h
${CMAKE_SOURCE_DIR}/clipper/core/rotation.h
${CMAKE_SOURCE_DIR}/clipper/core/clipper_test.h
${CMAKE_SOURCE_DIR}/clipper/core/hkl_datatypes.h
${CMAKE_SOURCE_DIR}/clipper/core/spacegroup.h
${CMAKE_SOURCE_DIR}/clipper/core/clipper_sysdep.h
${CMAKE_SOURCE_DIR}/clipper/core/clipper_types.h
${CMAKE_SOURCE_DIR}/clipper/core/hkl_info.h
${CMAKE_SOURCE_DIR}/clipper/core/spacegroup_data.h
${CMAKE_SOURCE_DIR}/clipper/core/clipper_util.h
${CMAKE_SOURCE_DIR}/clipper/core/hkl_lookup.h
${CMAKE_SOURCE_DIR}/clipper/core/symop.h
${CMAKE_SOURCE_DIR}/clipper/core/container.h
${CMAKE_SOURCE_DIR}/clipper/core/hkl_operators.h
${CMAKE_SOURCE_DIR}/clipper/core/container_hkl.h
${CMAKE_SOURCE_DIR}/clipper/core/map_interp.h
${CMAKE_SOURCE_DIR}/clipper/core/container_map.h
${CMAKE_SOURCE_DIR}/clipper/core/map_utils.h
${CMAKE_SOURCE_DIR}/clipper/core/xmap.h
${CMAKE_SOURCE_DIR}/clipper/core/container_types.h
${CMAKE_SOURCE_DIR}/clipper/core/nxmap.h
${CMAKE_SOURCE_DIR}/clipper/core/clipper_thread.h
${CMAKE_SOURCE_DIR}/clipper/core/test_core.h
${CMAKE_SOURCE_DIR}/clipper/core/test_data.h
)

set(clipper-contrib_sources
${CMAKE_SOURCE_DIR}/clipper/contrib/convolution_search.cpp
${CMAKE_SOURCE_DIR}/clipper/contrib/sfcalc.cpp
${CMAKE_SOURCE_DIR}/clipper/contrib/edcalc.cpp
${CMAKE_SOURCE_DIR}/clipper/contrib/sfcalc_obs.cpp
${CMAKE_SOURCE_DIR}/clipper/contrib/fffear.cpp
${CMAKE_SOURCE_DIR}/clipper/contrib/sfscale.cpp
${CMAKE_SOURCE_DIR}/clipper/contrib/function_object_bases.cpp
${CMAKE_SOURCE_DIR}/clipper/contrib/sfweight.cpp
${CMAKE_SOURCE_DIR}/clipper/contrib/mapfilter.cpp
${CMAKE_SOURCE_DIR}/clipper/contrib/skeleton.cpp
${CMAKE_SOURCE_DIR}/clipper/contrib/originmatch.cpp
${CMAKE_SOURCE_DIR}/clipper/contrib/test_contrib.cpp
)

set(clipper-contrib_headers
${CMAKE_SOURCE_DIR}/clipper/contrib/convolution_search.h
${CMAKE_SOURCE_DIR}/clipper/contrib/sfcalc.h
${CMAKE_SOURCE_DIR}/clipper/contrib/edcalc.h
${CMAKE_SOURCE_DIR}/clipper/contrib/sfcalc_obs.h
${CMAKE_SOURCE_DIR}/clipper/contrib/fffear.h
${CMAKE_SOURCE_DIR}/clipper/contrib/sfscale.h
${CMAKE_SOURCE_DIR}/clipper/contrib/function_object_bases.h
${CMAKE_SOURCE_DIR}/clipper/contrib/sfweight.h
${CMAKE_SOURCE_DIR}/clipper/contrib/mapfilter.h
${CMAKE_SOURCE_DIR}/clipper/contrib/skeleton.h
${CMAKE_SOURCE_DIR}/clipper/contrib/originmatch.h
${CMAKE_SOURCE_DIR}/clipper/contrib/test_contrib.h
)

set(clipper-phs_sources
${CMAKE_SOURCE_DIR}/clipper/phs/phs_io.cpp
)

set(clipper-phs_headers
${CMAKE_SOURCE_DIR}/clipper/phs/phs_io.h
)

set(clipper-cns_sources
${CMAKE_SOURCE_DIR}/clipper/cns/cns_hkl_io.cpp
${CMAKE_SOURCE_DIR}/clipper/cns/cns_map_io.cpp
)

set(clipper-cns_headers
${CMAKE_SOURCE_DIR}/clipper/cns/cns_hkl_io.h
${CMAKE_SOURCE_DIR}/clipper/cns/cns_map_io.h
)

set(clipper-mmdb_sources
${CMAKE_SOURCE_DIR}/clipper/mmdb/clipper_mmdb.cpp
)

set(clipper-mmdb_headers
${CMAKE_SOURCE_DIR}/clipper/mmdb/clipper_mmdb.h
)

set(clipper-cif_sources
${CMAKE_SOURCE_DIR}/clipper/cif/cif_data_io.cpp
)

set(clipper-cif_headers
${CMAKE_SOURCE_DIR}/clipper/cif/cif_data_io.h
)

set(clipper-ccp4_sources
${CMAKE_SOURCE_DIR}/clipper/ccp4/ccp4_map_io.cpp
${CMAKE_SOURCE_DIR}/clipper/ccp4/ccp4_mtz_io.cpp
${CMAKE_SOURCE_DIR}/clipper/ccp4/ccp4_mtz_types.cpp
${CMAKE_SOURCE_DIR}/clipper/ccp4/ccp4_utils.cpp
)

set(clipper-ccp4_headers
${CMAKE_SOURCE_DIR}/clipper/ccp4/ccp4_map_io.h
${CMAKE_SOURCE_DIR}/clipper/ccp4/ccp4_mtz_io.h
${CMAKE_SOURCE_DIR}/clipper/ccp4/ccp4_mtz_types.h
${CMAKE_SOURCE_DIR}/clipper/ccp4/ccp4_utils.h
)

# added clipper-gemmi files
set(clipper-gemmi_sources
${CMAKE_SOURCE_DIR}/clipper/gemmi/clipper_gemmi.cpp
${CMAKE_SOURCE_DIR}/clipper/gemmi/clipper_gemmi_model.cpp
)

set(clipper-gemmi_headers
${CMAKE_SOURCE_DIR}/clipper/gemmi/clipper_gemmi.h
${CMAKE_SOURCE_DIR}/clipper/gemmi/clipper_gemmi_model.h
)

# changed minimol_io to minimol_io_seq, minimol_io_mmdb, minimol_io_gemmi
set(clipper-minimol_sources
${CMAKE_SOURCE_DIR}/clipper/minimol/container_minimol.cpp
${CMAKE_SOURCE_DIR}/clipper/minimol/minimol.cpp	
${CMAKE_SOURCE_DIR}/clipper/minimol/minimol_data.cpp
${CMAKE_SOURCE_DIR}/clipper/minimol/minimol_io_seq.cpp
${CMAKE_SOURCE_DIR}/clipper/minimol/minimol_io_gemmi.cpp
${CMAKE_SOURCE_DIR}/clipper/minimol/minimol_io_mmdb.cpp
${CMAKE_SOURCE_DIR}/clipper/minimol/minimol_seq.cpp
${CMAKE_SOURCE_DIR}/clipper/minimol/minimol_utils.cpp
)

set(clipper-minimol_headers
${CMAKE_SOURCE_DIR}/clipper/minimol/container_minimol.h
${CMAKE_SOURCE_DIR}/clipper/minimol/minimol.h
${CMAKE_SOURCE_DIR}/clipper/minimol/minimol_data.h
${CMAKE_SOURCE_DIR}/clipper/minimol/minimol_io_gemmi.h
${CMAKE_SOURCE_DIR}/clipper/minimol/minimol_io_mmdb.h
${CMAKE_SOURCE_DIR}/clipper/minimol/minimol_io_seq.h
${CMAKE_SOURCE_DIR}/clipper/minimol/minimol_seq.h
${CMAKE_SOURCE_DIR}/clipper/minimol/minimol_utils.h
)

