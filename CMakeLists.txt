#find_package(MPI REQUIRED)
#add_executable(my_mpi_bin src1.cpp src2.cpp)
#target_include_directories(my_mpi_bin PRIVATE ${MPI_CXX_INCLUDE_PATH} src1.h src2.h)
#target_compile_options(my_mpi_bin PRIVATE ${MPI_CXX_COMPILE_FLAGS} my_compile_flags)
#target_link_libraries(my_mpi_bin ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS} my_link_flags)

cmake_minimum_required(VERSION 3.14)
project(conv_rewrite)

find_package(MPI REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_CXX_STANDARD 17)
message(${MPI_INCLUDE_PATH})

include_directories(${MPI_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_executable(conv_rewrite main.cpp)
target_link_libraries(conv_rewrite ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
    set_target_properties(conv_rewrite PROPERTIES
            COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set_target_properties(conv_rewrite PROPERTIES
            LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()