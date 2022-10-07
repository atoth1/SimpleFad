set(UNIT unit_tests)
set(SOURCES
  FadVariableTests.cpp
  FadLiteralTests.cpp
  BinaryExprTests.cpp
  UnaryExprTests.cpp
  ComparisonOpsTests.cpp
  )
add_executable(${UNIT} ${SOURCES})
target_link_libraries(${UNIT} ${PROJECT_NAME} GTest::GTest GTest::Main)
gtest_discover_tests(${UNIT})

set(EVAL evaluation_tests)
add_executable(${EVAL} EvaluationTests.cpp)
target_link_libraries(${EVAL} ${PROJECT_NAME} GTest::GTest GTest::Main)
gtest_discover_tests(${EVAL})

if(LAPACK_FOUND AND ENABLE_NEWTON_TESTS)
  set(NEWTON newton_tests)
  add_executable(${NEWTON} NewtonsMethodTests.cpp)
  target_link_libraries(${NEWTON} ${PROJECT_NAME} GTest::GTest GTest::Main LAPACK::LAPACK)
  gtest_discover_tests(${NEWTON})
endif()
