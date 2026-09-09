cmake_minimum_required(VERSION 3.18)

# flake.nix is the canonical version source (bumped there only), so parse the
# version line out of it rather than keeping a second copy that drifts.
# Expected line format:  version = "<YY.M.DD>";
file(STRINGS "${CMAKE_SOURCE_DIR}/flake.nix" _au_version_lines
     REGEX "^ *version = \"")
if(NOT _au_version_lines)
    message(FATAL_ERROR "Could not find a version line (version = \"...\") in flake.nix")
endif()
list(GET _au_version_lines 0 _au_version_line)
string(REGEX REPLACE "^ *version = \"([^\"]*)\"; *$" "\\1" AU_VERSION "${_au_version_line}")
if(AU_VERSION STREQUAL _au_version_line)
    message(FATAL_ERROR "Could not parse version from flake.nix: ${_au_version_line}")
endif()

project(analysis-utilities VERSION ${AU_VERSION} LANGUAGES CXX)

option(AU_USE_CUDA "Build with CUDA acceleration for RooFit PDFs" OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Match the doEval hot-loop flags the Makefile was using
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

find_package(OpenMP REQUIRED)
find_package(ROOT 6.38 REQUIRED COMPONENTS RooFit RooFitCore)
include(${ROOT_USE_FILE})

if(AU_USE_CUDA)
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 17)
    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        set(CMAKE_CUDA_ARCHITECTURES 89)
    endif()
    find_package(CUDAToolkit REQUIRED)
    add_compile_definitions(AU_ROOFIT_BACKEND_CUDA=1)
endif()

file(GLOB CPP_SOURCES CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/src/*.cpp)
set(CU_SOURCES "")
if(AU_USE_CUDA)
    file(GLOB CU_SOURCES CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/gpu/*.cu)
endif()

set(LIB_NAME analysis-utils)

# ROOT_GENERATE_DICTIONARY locates headers via INCLUDE_DIRECTORIES on the
# current directory, not via target properties — set it before the call.
include_directories(${CMAKE_SOURCE_DIR}/include)

# Generate ROOT dictionary for custom RooAbsPdf classes
ROOT_GENERATE_DICTIONARY(${LIB_NAME}Dict
    RooFitPhotopeakPdfs.hpp
    LINKDEF include/LinkDef.h
    MODULE ${LIB_NAME}
)

# Shared library (with dictionary)
add_library(${LIB_NAME} SHARED
    ${CPP_SOURCES}
    ${CU_SOURCES}
    ${CMAKE_CURRENT_BINARY_DIR}/${LIB_NAME}Dict.cxx
)
target_include_directories(${LIB_NAME} PUBLIC
    ${CMAKE_SOURCE_DIR}/include
    ${ROOT_INCLUDE_DIRS}
)
target_link_libraries(${LIB_NAME} PUBLIC
    ROOT::RooFit
    ROOT::RooFitCore
    OpenMP::OpenMP_CXX
)
if(AU_USE_CUDA)
    target_link_libraries(${LIB_NAME} PUBLIC CUDA::cudart)
endif()

# Static library — no dictionary (bundling it here too would risk double
# TClassTable registration if a downstream links both)
add_library(${LIB_NAME}-static STATIC
    ${CPP_SOURCES}
    ${CU_SOURCES}
)
target_include_directories(${LIB_NAME}-static PUBLIC
    ${CMAKE_SOURCE_DIR}/include
    ${ROOT_INCLUDE_DIRS}
)
target_link_libraries(${LIB_NAME}-static PUBLIC
    ROOT::RooFit
    ROOT::RooFitCore
    OpenMP::OpenMP_CXX
)
if(AU_USE_CUDA)
    target_link_libraries(${LIB_NAME}-static PUBLIC CUDA::cudart)
endif()
set_target_properties(${LIB_NAME}-static PROPERTIES OUTPUT_NAME ${LIB_NAME})

# Public headers
file(GLOB HEADERS ${CMAKE_SOURCE_DIR}/include/*.hpp)

include(GNUInstallDirs)
install(TARGETS ${LIB_NAME} ${LIB_NAME}-static
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/lib${LIB_NAME}_rdict.pcm
    ${CMAKE_CURRENT_BINARY_DIR}/lib${LIB_NAME}.rootmap
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

# pkgconfig — substitute @PREFIX@ and @LIBNAME@ in the .pc.in template
set(PREFIX ${CMAKE_INSTALL_PREFIX})
set(LIBNAME ${LIB_NAME})
configure_file(
    ${CMAKE_SOURCE_DIR}/analysis-utilities.pc.in
    ${CMAKE_CURRENT_BINARY_DIR}/analysis-utilities.pc
    @ONLY
)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/analysis-utilities.pc
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
)

# ---------------------------------------------------------------------------
# API documentation
#
# Doxygen is not needed to build or install the library, so the target is only
# defined when doxygen is present — `nix build` never sees it. graphviz is
# likewise optional; without dot the inheritance and include graphs are simply
# left out of the generated HTML.
#
#   cmake --build build --target docs
# ---------------------------------------------------------------------------
option(AU_BUILD_DOCS "Define the 'docs' target when Doxygen is available" ON)
if(AU_BUILD_DOCS)
    find_package(Doxygen OPTIONAL_COMPONENTS dot)
endif()
if(AU_BUILD_DOCS AND DOXYGEN_FOUND)
    # Default to the in-tree docs/ that .gitignore already covers, so the
    # generated HTML is one relative path away from the source being read.
    set(AU_DOCS_OUTPUT_DIR "${CMAKE_SOURCE_DIR}/docs"
        CACHE PATH "Directory the 'docs' target writes generated API documentation to")

    # doxygen-awesome-css ships only data files, so it cannot be found on PATH.
    # The dev shells export DOXYGEN_AWESOME_CSS; packages.docs sets it directly.
    set(AU_DOXYGEN_AWESOME_DIR "$ENV{DOXYGEN_AWESOME_CSS}"
        CACHE PATH "Directory containing doxygen-awesome.css")

    set(AU_SOURCE_DIR ${CMAKE_SOURCE_DIR})
    set(AU_DOXYGEN_OUTPUT_DIR ${AU_DOCS_OUTPUT_DIR})
    if(TARGET Doxygen::dot)
        set(AU_DOXYGEN_HAVE_DOT YES)
    else()
        set(AU_DOXYGEN_HAVE_DOT NO)
        message(STATUS "graphviz (dot) not found — API docs will omit inheritance and include graphs")
    endif()

    if(NOT EXISTS "${AU_DOXYGEN_AWESOME_DIR}/doxygen-awesome.css")
        message(FATAL_ERROR
            "doxygen-awesome-css not found at '${AU_DOXYGEN_AWESOME_DIR}'. "
            "Enter a dev shell (nix develop), or pass "
            "-DAU_DOXYGEN_AWESOME_DIR=<dir containing doxygen-awesome.css>. "
            "Configure with -DAU_BUILD_DOCS=OFF to skip the docs target entirely.")
    endif()

    configure_file(
        ${CMAKE_SOURCE_DIR}/Doxyfile.in
        ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        @ONLY
    )

    add_custom_target(docs
        COMMAND ${CMAKE_COMMAND} -E make_directory ${AU_DOXYGEN_OUTPUT_DIR}
        COMMAND Doxygen::doxygen ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        COMMAND ${CMAKE_COMMAND} -E echo "API documentation: ${AU_DOXYGEN_OUTPUT_DIR}/html/index.html"
        COMMAND ${CMAKE_COMMAND} -E echo "Doxygen warnings:  ${AU_DOXYGEN_OUTPUT_DIR}/doxygen-warnings.log"
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM
    )
elseif(AU_BUILD_DOCS)
    message(STATUS "Doxygen not found — the 'docs' target is unavailable")
endif()
