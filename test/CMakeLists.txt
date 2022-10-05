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
