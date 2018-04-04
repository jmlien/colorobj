#####################################################
#                COMMON Libraries                   #
#####################################################
set(APP_NAME colorobj)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

set(APP_SRC ${APP_SRC} main.cpp dump.cpp flattened_mesh.cpp angular-map.cpp draw.cpp)

find_package(CGAL REQUIRED)
find_package(CGAL REQUIRED COMPONENTS Core)
find_package(Eigen3 REQUIRED)
find_package(OpenGL REQUIRED)
find_package(glfw3 3.2 REQUIRED)

# search for pkg-config
include (FindPkgConfig)
if (NOT PKG_CONFIG_FOUND)
    message (FATAL_ERROR "pkg-config not found")
endif ()

# check for libpng
pkg_check_modules (LIBPNG libpng16 REQUIRED)
if (NOT LIBPNG_FOUND)
    message(FATAL_ERROR "You don't seem to have libpng16 development libraries installed")
else ()
    include_directories (${LIBPNG_INCLUDE_DIRS})
    link_directories (${LIBPNG_LIBRARY_DIRS})
    link_libraries (${LIBPNG_LIBRARIES})
endif ()


include( ${CGAL_USE_FILE} )
include( ${EIGEN3_USE_FILE} )


########## add subdirectories ############
add_subdirectory(mascgl-3)


set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
find_package(Boost REQUIRED COMPONENTS
             program_options
             system
             filesystem)

add_executable(${APP_NAME} ${APP_SRC})

include_directories(
      ${EIGEN3_INCLUDE_DIRS}
      ${CGAL_INCLUDE_DIRS}
      ${Boost_INCLUDE_DIR}
      ${GMP_INCLUDE_DIRS}
      ${PNG_INCLUDE_DIR}
      )

include_directories(
    "${CMAKE_SOURCE_DIR}/libs/FIST/include"
    #    "${CMAKE_SOURCE_DIR}/src"
       "${CMAKE_SOURCE_DIR}/src/mascgl-3"
       "${CMAKE_SOURCE_DIR}/src/mascgl-3/mathtool"
       "${CMAKE_SOURCE_DIR}/src/mascgl-3/modelgraph"
       "${CMAKE_SOURCE_DIR}/src/mascgl-3/soil")

link_directories(${Boost_LIBRARY_DIR})


target_link_libraries(${APP_NAME}
                    ${Boost_LIBRARIES}
                    ${GMP_LIBRARIES}  ${PNG_LIBRARY}
                    ${CGAL_LIBRARIES} glfw ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY})
