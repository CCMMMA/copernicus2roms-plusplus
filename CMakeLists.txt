cmake_minimum_required(VERSION 3.20)
project(MyOcean2ROMS)

set(CMAKE_CXX_STANDARD 23)

find_package(OpenMP REQUIRED)
set(LIBOMP OpenMP::OpenMP_CXX)

include(ExternalProject)
set(EXTERNAL_INSTALL_LOCATION ${CMAKE_BINARY_DIR}/external)

set(NETCDFCXXLIBS_PARAMS "")
if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(NETCDFCXXLIBS_PARAMS "-lnetcdf -lhdf5_hl -lhdf5 -lm -lz -lcurl")
endif()

include_directories(SYSTEM ${EXTERNAL_INSTALL_LOCATION}/include)
link_directories(${EXTERNAL_INSTALL_LOCATION}/lib)

ExternalProject_Add(szlib
        URL https://support.hdfgroup.org/ftp/lib-external/szip/2.1.1/src/szip-2.1.1.tar.gz
        TIMEOUT 360
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND ./configure CFLAGS=-fPIC CPPFLAGS=-I${EXTERNAL_INSTALL_LOCATION}/include/ LDFLAGS=-L${EXTERNAL_INSTALL_LOCATION}/lib/ --prefix=${EXTERNAL_INSTALL_LOCATION}
        BUILD_COMMAND make
        INSTALL_COMMAND make install
        )
set(LIBSZLIB ${EXTERNAL_INSTALL_LOCATION}/lib/libsz.a)
if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(LIBSZLIB ${EXTERNAL_INSTALL_LOCATION}/lib/libsz.dylib)
endif()

ExternalProject_Add(zlib
        URL https://zlib.net/zlib-1.3.1.tar.gz
        TIMEOUT 360
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND CFLAGS=-fPIC ./configure --prefix=${EXTERNAL_INSTALL_LOCATION}
        BUILD_COMMAND make
        INSTALL_COMMAND make install
        )
set(LIBZLIB ${EXTERNAL_INSTALL_LOCATION}/lib/libz.a)
if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(LIBZLIB ${EXTERNAL_INSTALL_LOCATION}/lib/libz.dylib)
endif()

ExternalProject_Add(curl
        URL https://curl.haxx.se/download/curl-7.73.0.tar.gz
        TIMEOUT 360
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND ./configure --prefix=${EXTERNAL_INSTALL_LOCATION} CFLAGS=-fPIC CPPFLAGS=-I${EXTERNAL_INSTALL_LOCATION}/include/ LDFLAGS=-L${EXTERNAL_INSTALL_LOCATION}/lib/ --without-libidn --disable-ldap --without-brotli --without-openssl
        BUILD_COMMAND make
        INSTALL_COMMAND make install
        )
set(LIBCURL ${EXTERNAL_INSTALL_LOCATION}/lib/libcurl.a)
if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(LIBCURL ${EXTERNAL_INSTALL_LOCATION}/lib/libcurl.dylib)
endif()

ExternalProject_Add(hdf5
        URL https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.12/hdf5-1.12.0/src/hdf5-1.12.0.tar.gz
        TIMEOUT 360
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND ./configure --prefix=${EXTERNAL_INSTALL_LOCATION} CFLAGS=-fPIC CPPFLAGS=-I${EXTERNAL_INSTALL_LOCATION}/include/ LDFLAGS=-L${EXTERNAL_INSTALL_LOCATION}/lib/ --enable-hl --enable-shared --enable-build-mode=production --enable-unsupported --enable-cxx --with-zlib=${EXTERNAL_INSTALL_LOCATION} --with-szlib=${EXTERNAL_INSTALL_LOCATION} --enable-threadsafe --with-pthread
        BUILD_COMMAND make
        INSTALL_COMMAND make install
        )
set(LIBHDF5 ${EXTERNAL_INSTALL_LOCATION}/lib/libhdf5_hl.a;${EXTERNAL_INSTALL_LOCATION}/lib/libhdf5.a)
if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(LIBHDF5 ${EXTERNAL_INSTALL_LOCATION}/lib/libhdf5.dylib)
endif()

ExternalProject_Add(netcdf
        URL https://downloads.unidata.ucar.edu/netcdf-c/4.8.1/netcdf-c-4.8.1.tar.gz
        TIMEOUT 360
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND ./configure CFLAGS=-fPIC CPPFLAGS=-I${EXTERNAL_INSTALL_LOCATION}/include/ LDFLAGS=-L${EXTERNAL_INSTALL_LOCATION}/lib/ --prefix=${EXTERNAL_INSTALL_LOCATION}  --enable-shared --enable-netcdf-4 --enable-dap --enable-byterange --enable-erange-fill
        INSTALL_COMMAND make install
        )
set(LIBNETCDF ${EXTERNAL_INSTALL_LOCATION}/lib/libnetcdf.a)
if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(LIBNETCDF ${EXTERNAL_INSTALL_LOCATION}/lib/libnetcdf.dylib)
endif()

ExternalProject_Add(netcdfcxx
        URL https://downloads.unidata.ucar.edu/netcdf-cxx/4.3.1/netcdf-cxx4-4.3.1.tar.gz
        # URL /Users/ciro/Downloads/netcdf-cxx4-4.3.1.zip
        TIMEOUT 360
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND ./configure LIBS=${NETCDFCXXLIBS_PARAMS} CFLAGS=-fPIC CPPFLAGS=-I${EXTERNAL_INSTALL_LOCATION}/include/ LDFLAGS=-L${EXTERNAL_INSTALL_LOCATION}/lib/ --prefix=${EXTERNAL_INSTALL_LOCATION}   --enable-shared
        INSTALL_COMMAND make install
        )
set(LIBNETCDFCXX ${EXTERNAL_INSTALL_LOCATION}/lib/libnetcdf_c++4.a)
if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(LIBNETCDFCXX ${EXTERNAL_INSTALL_LOCATION}/lib/libnetcdf_c++4.dylib)
endif()

add_executable(MyOcean2ROMS MyOcean2ROMS.cpp jncregridder/roms/ROMSGrid.cpp jncregridder/roms/ROMSGrid.h jncregridder/data/copernicus/CopernicusBase.cpp jncregridder/data/copernicus/CopernicusBase.h jncregridder/data/copernicus/CopernicusTem.cpp jncregridder/data/copernicus/CopernicusTem.h jncregridder/data/copernicus/CopernicusCur.cpp jncregridder/data/copernicus/CopernicusCur.h jncregridder/data/copernicus/CopernicusSSH.cpp jncregridder/data/copernicus/CopernicusSSH.h jncregridder/data/copernicus/CopernicusSal.cpp jncregridder/data/copernicus/CopernicusSal.h jncregridder/util/BilinearInterpolator.cpp jncregridder/util/BilinearInterpolator.h jncregridder/util/InterpolatorBase.cpp jncregridder/util/InterpolatorBase.h jncregridder/roms/ROMSInit.cpp jncregridder/roms/ROMSInit.h jncregridder/util/BilinearInterpolator3D.cpp jncregridder/util/BilinearInterpolator3D.h jncregridder/roms/ROMSBoundary.cpp jncregridder/roms/ROMSBoundary.h)

add_dependencies(zlib szlib)
add_dependencies(curl zlib)
add_dependencies(hdf5 curl)
add_dependencies(netcdf hdf5)
add_dependencies(netcdfcxx netcdf)
add_dependencies(${PROJECT_NAME} netcdfcxx)

target_link_libraries(${PROJECT_NAME} ${LIBOMP} ${LIBNETCDFCXX} ${LIBNETCDF} ${LIBCURL} ${LIBHDF5} ${LIBZLIB} ${LIBSZLIB})
