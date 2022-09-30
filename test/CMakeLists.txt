set(NAME unit_tests)
set(SOURCES
  FadVariableTests.cpp
  FadLiteralTests.cpp
  BinaryExprTests.cpp
  UnaryExprTests.cpp
  ComparisonOpsTests.cpp
  )
add_executable(${NAME} ${SOURCES})
target_link_libraries(${NAME} ${PROJECT_NAME} GTest::GTest GTest::Main)
gtest_discover_tests(${NAME})
