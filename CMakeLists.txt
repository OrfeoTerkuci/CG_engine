cmake_minimum_required(VERSION 3.5)

project(engine)

############################################################
# Set compiler flags
############################################################
# If supported by your compiler, you can add the -Wall, -Wextra, –fstack-protector-all and -g3 flags here.
set(OWN_GXX_FLAGS "-std=c++14")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${OWN_GXX_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OWN_GXX_FLAGS}")
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${OWN_GXX_FLAGS}")

############################################################
# List all sources
############################################################
set(engine_sources
        easy_image.cc
        engine.cc
        ini_configuration.cc
        l_parser/l_parser.cc
        vector/vector3d.cc
        )
############################################################
# Create main executable
############################################################
set(exe_name "engine")
add_executable( ${exe_name} ${engine_sources} )
############################################################
# Create test specific executables
############################################################
add_executable( "l_systems" ${engine_sources} )
set_target_properties("l_systems" PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/ini/l_systems" )

add_executable( "line_drawings" ${engine_sources} )
set_target_properties("line_drawings" PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/ini/line_drawings" )

add_executable( "wireframes" ${engine_sources} )
set_target_properties("wireframes" PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/ini/wireframes" )

add_executable( "z_buffered_wireframes" ${engine_sources} )
set_target_properties("z_buffered_wireframes" PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/ini/z_buffered_wireframes" )

add_executable( "z_buffering" ${engine_sources} )
set_target_properties("z_buffering" PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/ini/z_buffering" )

add_executable( "3d_fractals" ${engine_sources} )
set_target_properties("3d_fractals" PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/ini/3d_fractals" )

add_executable( "clipping" ${engine_sources} )
set_target_properties("3d_fractals" PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/ini/clipping" )

add_executable( "ambient" ${engine_sources} )
set_target_properties("3d_fractals" PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/ini/ambient" )

add_executable( "ambient_clipping" ${engine_sources} )
set_target_properties("3d_fractals" PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/ini/ambient_clipping" )

add_executable( "diffuse" ${engine_sources} )
set_target_properties("3d_fractals" PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/ini/diffuse" )

add_executable( "diffuse_clipping" ${engine_sources} )
set_target_properties("3d_fractals" PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/ini/diffuse_clipping" )

add_executable( "specular" ${engine_sources} )
set_target_properties("3d_fractals" PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/ini/specular" )

add_executable( "specular_clipping" ${engine_sources} )
set_target_properties("3d_fractals" PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/ini/specular_clipping" )

install( TARGETS ${exe_name} DESTINATION ${PROJECT_SOURCE_DIR}/ )
