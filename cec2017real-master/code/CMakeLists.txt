cmake_minimum_required(VERSION 3.0)
PROJECT(cec17)
ADD_EXECUTABLE(test "test.cc")
ADD_EXECUTABLE(testrandom "testrandom.cc")
ADD_EXECUTABLE(testsolis "testsolis.cc")
ADD_EXECUTABLE(ARO "ARO.cpp")
ADD_EXECUTABLE(multiARO "multiARO.cpp")
ADD_EXECUTABLE(memeticARO "memeticARO.cpp")
ADD_LIBRARY("cec17_test_func" SHARED "cec17_test_func.c" "cec17.c")
include_directories(".")
TARGET_LINK_LIBRARIES(test "cec17_test_func")
TARGET_LINK_LIBRARIES(testrandom "cec17_test_func")
TARGET_LINK_LIBRARIES(testsolis "cec17_test_func")
TARGET_LINK_LIBRARIES(ARO "cec17_test_func")
TARGET_LINK_LIBRARIES(multiARO "cec17_test_func")
TARGET_LINK_LIBRARIES(memeticARO "cec17_test_func")

