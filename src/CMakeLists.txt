add_library(${PROJECT_NAME} INTERFACE)

configure_file(SimpleFad_config.hpp.in SimpleFad_config.hpp)

target_include_directories(${PROJECT_NAME}
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
  
target_include_directories(${PROJECT_NAME}
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)
  
target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)

install(
  FILES 
    "SimpleFad.hpp"
    "SimpleFad_BinaryExpr.hpp"
    "SimpleFad_BinaryOps.hpp"
    "SimpleFad_ComparisonOps.hpp"
    "SimpleFad_DerivStorageTraits.hpp"
    "SimpleFad_Exception.hpp"
    "SimpleFad_ExprBase.hpp"
    "SimpleFad_ExprTraits.hpp"
    "SimpleFad_FadLiteral.hpp"
    "SimpleFad_FadVariable.hpp"
    "SimpleFad_UnaryExpr.hpp"
    "SimpleFad_UnaryOps.hpp"
    "${CMAKE_CURRENT_BINARY_DIR}/SimpleFad_config.hpp"
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
