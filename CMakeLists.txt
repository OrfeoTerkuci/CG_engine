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
set(geometry_sources
        src/geometry/shapes.cpp
)

set(generator_sources
        src/generators/basic2d.cpp
        src/generators/draw2d.cpp
        src/generators/draw3d.cpp
        src/generators/fractals.cpp
        src/generators/lSystems3d.cpp
        src/generators/platonicBodies.cpp
        src/generators/zBuffering.cpp
        src/generators/wireframes.cpp
)

set(tooling_sources
        src/clipping.cpp
        src/transformations.cpp
        src/color/color.cpp
        src/lighting/light.cpp
)


set(engine_sources
        src/easy_image.cc
        src/engine.cc
        src/ini_configuration.cc
        src/l_parser/l_parser.cc
        src/vector/vector3d.cc
        ${geometry_sources}
        ${generator_sources}
        ${tooling_sources}
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
set_target_properties("clipping" PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/ini/clipping" )

add_executable( "ambient" ${engine_sources} )
set_target_properties("ambient" PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/ini/ambient" )

add_executable( "ambient_clipping" ${engine_sources} )
set_target_properties("ambient_clipping" PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/ini/ambient_clipping" )

add_executable( "diffuse" ${engine_sources} )
set_target_properties("diffuse" PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/ini/diffuse" )

add_executable( "diffuse_clipping" ${engine_sources} )
set_target_properties("diffuse_clipping" PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/ini/diffuse_clipping" )

add_executable( "specular" ${engine_sources} )
set_target_properties("specular" PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/ini/specular" )

add_executable( "specular_clipping" ${engine_sources} )
set_target_properties("specular_clipping" PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/ini/specular_clipping" )

add_executable( "shadowing" ${engine_sources} )
set_target_properties("shadowing" PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/ini/shadowing" )

add_executable( "shadowing_clipping" ${engine_sources} )
set_target_properties("shadowing_clipping" PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/ini/shadowing_clipping" )

add_executable( "obj_parser_en_voorbeeld" ${engine_sources} )
set_target_properties("obj_parser_en_voorbeeld" PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/ini/obj_parser_en_voorbeeld" )

add_executable( "spheres_and_cylinders" ${engine_sources} )
set_target_properties("spheres_and_cylinders" PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/ini/spheres_and_cylinders" )

install( TARGETS ${exe_name} DESTINATION ${PROJECT_SOURCE_DIR}/ )
