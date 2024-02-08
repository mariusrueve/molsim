cmake_minimum_required(VERSION 3.18)

project(simple_cxx_example)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(RDKit REQUIRED)
find_package(Boost COMPONENTS timer system program_options REQUIRED)
find_package(spdlog REQUIRED)

add_executable(molsim molsim.cpp)
target_link_libraries(molsim PUBLIC Boost::timer Boost::system Boost::program_options spdlog::spdlog 
RDKit::FileParsers)