PROJECT (RBDLEXAMPLE CXX)

CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

# We need to add the project source path to the CMake module path so that
# the FindRBDL.cmake script can be found.
LIST( APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR} )


SET(CUSTOM_RBDL_PATH "" CACHE PATH "Path to specific RBDL Installation")

# Search for the RBDL include directory and library
FIND_PACKAGE (NLopt REQUIRED)
FIND_PACKAGE (RBDL COMPONENTS URDFREADER REQUIRED)
FIND_PACKAGE (Eigen3 REQUIRED)
FIND_PACKAGE (Boost REQUIRED)
FIND_PACKAGE (pinocchio REQUIRED)
FIND_PACKAGE (urdfdom REQUIRED)


# Add the include directory to the include paths
INCLUDE_DIRECTORIES (PRIVATE ${NLOPT_PRIVATE_INCLUDE_DIRS} ${RBDL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${PINOCCHIO_INCLUDE_DIR} ${Boost_INCLUDE_DIRS} ${pinocchio_INCLUDE_DIRS})

# Create an executable
ADD_EXECUTABLE (pinocchio pinocchio_urdf.cpp)

# And link the library against the executable
TARGET_LINK_LIBRARIES ( pinocchio
	${NLOPT_LIBRARIES}
	${RBDL_LIBRARY}
	${RBDL_URDFREADER_LIBRARY}
	${PINOCCHIO_LIBRARY}
	${pinocchio_LIBRARY_DIRS}
	${Boost_LIBRARIES}
	)

	
target_compile_options(pinocchio PUBLIC ${PINOCCHIO_CFLAGS} "-DPINOCCHIO_WITH_URDFDOM")

