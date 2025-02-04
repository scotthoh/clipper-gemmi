cmake_minimum_required(VERSION 3.4)
project(clipper CXX)

include(GNUInstallDirs)

option(BUILD_SHARED_LIBS "build shared libs (true)" TRUE)

find_package(CCP4 REQUIRED ccp4c)
find_package(gemmi 0.6.7 CONFIG REQUIRED)

# FFTW2 - first search in the install path, then in the system
find_library(FFTW3_LIBRARY NAMES fftw3f
             HINTS ${LIB_INSTALL_DIR} ${CMAKE_INSTALL_PREFIX}/lib
             NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)
find_path(FFTW3_INCLUDE_DIRS fftw3.h
          HINTS ${CMAKE_INSTALL_PREFIX}/include)
message(STATUS "FFTW3 libraries - ${FFTW3_LIBRARY}")
message(STATUS "FFTW3 header directory - ${FFTW3_INCLUDE_DIRS}")

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads)

include(${CMAKE_SOURCE_DIR}/cmake/source.cmake)

set(PACKAGE "clipper")
set(VERSION "3.0.3-dev") # might need to bump major version to 3?
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmakein
               ${CMAKE_CURRENT_BINARY_DIR}/config.h @ONLY)

set(prefix      ${CMAKE_INSTALL_PREFIX})
set(exec_prefix ${CMAKE_INSTALL_PREFIX})
set(libdir      ${CMAKE_INSTALL_PREFIX}/lib)
set(includedir  ${CMAKE_INSTALL_PREFIX}/include)
set(CLIPPER_REQUIRES "ccp4c") #mmdb2")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/clipper.pc.cmakein
   ${CMAKE_CURRENT_BINARY_DIR}/clipper.pc.in @ONLY)
 file(GENERATE OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/clipper.pc
   INPUT ${CMAKE_CURRENT_BINARY_DIR}/clipper.pc.in)

add_library(clipper-core ${clipper-core_sources} ${clipper-core_headers})
add_library(clipper::core ALIAS clipper-core)

target_link_libraries(clipper-core PRIVATE 
                             ${FFTW3_LIBRARY}
)
if(Threads_FOUND)
target_link_libraries(clipper-core PRIVATE Threads::Threads)
endif()

target_include_directories(clipper-core PRIVATE .
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/clipper/core
    ${FFTW3_INCLUDE_DIRS}
    ${CMAKE_BINARY_DIR}
)

set_target_properties(clipper-core PROPERTIES PUBLIC_HEADER 
                      ${CMAKE_SOURCE_DIR}/clipper/clipper.h 
                      ${CMAKE_SOURCE_DIR}/clipper/minimal-clipper-hkl.h
                      ${CMAKE_SOURCE_DIR}/clipper/minimal-clipper-map.h
)
set_target_properties(clipper-core PROPERTIES PRIVATE_HEADER "${clipper-core_headers}")
set_target_properties(clipper-core PROPERTIES SOVERSION ${VERSION})
if(MSVC)
    set_target_properties(clipper-core PROPERTIES DEFINE_SYMBOL "CLIPPER_DLL_EXPORT")
endif()

install(TARGETS clipper-core
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/clipper/
      PRIVATE_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/clipper/core
)

add_library(clipper-contrib ${clipper-contrib_sources} ${clipper-contrib_headers})
add_library(clipper::contrib ALIAS clipper-contrib)

target_link_libraries(clipper-contrib PRIVATE 
                             clipper::core 
                             ${FFTW3_LIBRARY}
)
if(Threads_FOUND)
target_link_libraries(clipper-contrib PRIVATE Threads::Threads)
endif()

target_include_directories(clipper-contrib PRIVATE .
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/clipper/core
    ${CMAKE_SOURCE_DIR}/clipper/contrib
    ${FFTW3_INCLUDE_DIRS}
    ${CMAKE_BINARY_DIR}
)

set_target_properties(clipper-contrib PROPERTIES PUBLIC_HEADER "${CMAKE_SOURCE_DIR}/clipper/clipper-contrib.h")
set_target_properties(clipper-contrib PROPERTIES PRIVATE_HEADER "${clipper-contrib_headers}")
set_target_properties(clipper-contrib PROPERTIES SOVERSION ${VERSION})

install(TARGETS clipper-contrib
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/clipper/
      PRIVATE_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/clipper/contrib
)

add_library(clipper-cns ${clipper-cns_sources} ${clipper-cns_headers})
add_library(clipper::cns ALIAS clipper-cns)

target_link_libraries(clipper-cns PRIVATE 
                             clipper::core 
)
if(Threads_FOUND)
target_link_libraries(clipper-cns PRIVATE Threads::Threads)
endif()

target_include_directories(clipper-cns PRIVATE .
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/clipper/core
    ${CMAKE_SOURCE_DIR}/clipper/cns
    ${CMAKE_BINARY_DIR}
)

set_target_properties(clipper-cns PROPERTIES PUBLIC_HEADER "${CMAKE_SOURCE_DIR}/clipper/clipper-cns.h")
set_target_properties(clipper-cns PROPERTIES PRIVATE_HEADER "${clipper-cns_headers}")
set_target_properties(clipper-cns PROPERTIES SOVERSION ${VERSION})

install(TARGETS clipper-cns
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/clipper/
      PRIVATE_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/clipper/cns
)

add_library(clipper-phs ${clipper-phs_sources} ${clipper-phs_headers})
add_library(clipper::phs ALIAS clipper-phs)

target_link_libraries(clipper-phs PRIVATE 
                             clipper::core 
)
if(Threads_FOUND)
target_link_libraries(clipper-phs PRIVATE Threads::Threads)
endif()

target_include_directories(clipper-phs PRIVATE .
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/clipper/core
    ${CMAKE_SOURCE_DIR}/clipper/phs
    ${CMAKE_BINARY_DIR}
)

set_target_properties(clipper-phs PROPERTIES PUBLIC_HEADER "${CMAKE_SOURCE_DIR}/clipper/clipper-phs.h")
set_target_properties(clipper-phs PROPERTIES PRIVATE_HEADER "${clipper-phs_headers}")
set_target_properties(clipper-phs PROPERTIES SOVERSION ${VERSION})

install(TARGETS clipper-phs
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/clipper/
      PRIVATE_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/clipper/phs
)

#add_library(clipper-mmdb ${clipper-mmdb_sources} ${clipper-mmdb_headers})
#add_library(clipper::mmdb ALIAS clipper-mmdb)

#target_link_libraries(clipper-mmdb PRIVATE 
#                             clipper::core 
#                             ${MMDB2_LIBRARY}
#)
#if(Threads_FOUND)
#target_link_libraries(clipper-mmdb PRIVATE Threads::Threads)
#endif()

#target_include_directories(clipper-mmdb PRIVATE .
#    ${CMAKE_SOURCE_DIR}
#    ${CMAKE_SOURCE_DIR}/clipper/core
#    ${CMAKE_SOURCE_DIR}/clipper/mmdb
#    ${MMDB2_INCLUDE_DIR}
#    ${CCP4C_INCLUDE_DIR}
#    ${FFTW3_INCLUDE_DIRS}
#    ${CMAKE_BINARY_DIR}
#)

#set_target_properties(clipper-mmdb PROPERTIES PUBLIC_HEADER "${CMAKE_SOURCE_DIR}/clipper/clipper-mmdb.h")
#set_target_properties(clipper-mmdb PROPERTIES PRIVATE_HEADER "${clipper-mmdb_headers}")
#set_target_properties(clipper-mmdb PROPERTIES SOVERSION ${VERSION})

#install(TARGETS clipper-mmdb
#      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
#      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
#      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
#      PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/clipper/
#      PRIVATE_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/clipper/mmdb
#)

add_library(clipper-ccp4 ${clipper-ccp4_sources} ${clipper-ccp4_headers})
add_library(clipper::ccp4 ALIAS clipper-ccp4)

target_link_libraries(clipper-ccp4 PRIVATE 
                             clipper::core 
                             ${CCP4C_LIBRARY}
                             #${MMDB2_LIBRARY} # mmdb library might not be needed here?
)
if(Threads_FOUND)
target_link_libraries(clipper-ccp4 PRIVATE Threads::Threads)
endif()

target_include_directories(clipper-ccp4 PRIVATE .
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/clipper/core
    ${CMAKE_SOURCE_DIR}/clipper/ccp4
    #${MMDB2_INCLUDE_DIR}  # mmdb library might not be needed here?
    ${CCP4C_INCLUDE_DIR}
    ${CMAKE_BINARY_DIR}
)

set_target_properties(clipper-ccp4 PROPERTIES PUBLIC_HEADER "${CMAKE_SOURCE_DIR}/clipper/clipper-ccp4.h")
set_target_properties(clipper-ccp4 PROPERTIES PRIVATE_HEADER "${clipper-ccp4_headers}")
set_target_properties(clipper-ccp4 PROPERTIES SOVERSION ${VERSION})

install(TARGETS clipper-ccp4
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/clipper/
      PRIVATE_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/clipper/ccp4
)

# gemmi
add_library(clipper-gemmi ${clipper-gemmi_sources} ${clipper-gemmi_headers})
add_library(clipper::gemmi ALIAS clipper-gemmi)

target_link_libraries(clipper-gemmi PRIVATE 
                             clipper::core gemmi::gemmi_cpp)

if(Threads_FOUND)
target_link_libraries(clipper-gemmi PRIVATE Threads::Threads)
endif()

target_include_directories(clipper-gemmi PRIVATE .
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/clipper/core
    ${CMAKE_SOURCE_DIR}/clipper/gemmi
    ${CMAKE_BINARY_DIR}
)

set_target_properties(clipper-gemmi PROPERTIES PUBLIC_HEADER "${CMAKE_SOURCE_DIR}/clipper/clipper-gemmi.h")
set_target_properties(clipper-gemmi PROPERTIES PRIVATE_HEADER "${clipper-gemmi_headers}")
set_target_properties(clipper-gemmi PROPERTIES SOVERSION ${VERSION})

install(TARGETS clipper-gemmi
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/clipper/
      PRIVATE_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/clipper/gemmi
)
#gemmi

add_library(clipper-minimol ${clipper-minimol_sources} ${clipper-minimol_headers})
add_library(clipper::minimol ALIAS clipper-minimol)

target_link_libraries(clipper-minimol PRIVATE 
                             clipper::core 
			     #clipper::mmdb
                             clipper::gemmi
                             ${CCP4C_LIBRARY}
			     #${MMDB2_LIBRARY}
			                 gemmi::gemmi_cpp
)
if(Threads_FOUND)
target_link_libraries(clipper-minimol PRIVATE Threads::Threads)
endif()

target_include_directories(clipper-minimol PRIVATE .
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/clipper/core
    #${CMAKE_SOURCE_DIR}/clipper/mmdb
    ${CMAKE_SOURCE_DIR}/clipper/gemmi
    ${CMAKE_SOURCE_DIR}/clipper/minimol
    ${MMDB2_INCLUDE_DIR}
    ${CCP4C_INCLUDE_DIR}
    ${CMAKE_BINARY_DIR}
)

set_target_properties(clipper-minimol PROPERTIES PUBLIC_HEADER "${CMAKE_SOURCE_DIR}/clipper/clipper-minimol.h")
set_target_properties(clipper-minimol PROPERTIES PRIVATE_HEADER "${clipper-minimol_headers}")
set_target_properties(clipper-minimol PROPERTIES SOVERSION ${VERSION})
if(MSVC)
    set_target_properties(clipper-minimol PROPERTIES DEFINE_SYMBOL "CLIPPER_MINIMOL_DLL_EXPORT")
endif()

install(TARGETS clipper-minimol
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/clipper/
      PRIVATE_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/clipper/minimol
)

add_library(clipper-cif ${clipper-cif_sources} ${clipper-cif_headers})
add_library(clipper::cif ALIAS clipper-cif)

target_link_libraries(clipper-cif PRIVATE 
                             clipper::core 
                             gemmi::gemmi_cpp
)
if(Threads_FOUND)
target_link_libraries(clipper-cif PRIVATE Threads::Threads)
endif()

target_include_directories(clipper-cif PRIVATE .
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/clipper/core
    ${CMAKE_SOURCE_DIR}/clipper/cif
    ${CMAKE_BINARY_DIR}
)

set_target_properties(clipper-cif PROPERTIES PUBLIC_HEADER "${CMAKE_SOURCE_DIR}/clipper/clipper-cif.h")
set_target_properties(clipper-cif PROPERTIES PRIVATE_HEADER "${clipper-cif_headers}")
set_target_properties(clipper-cif PROPERTIES SOVERSION ${VERSION})

install(TARGETS clipper-cif
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/clipper/
      PRIVATE_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/clipper/cif
)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/clipper.pc DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

