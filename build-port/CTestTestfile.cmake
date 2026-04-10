# CMake generated Testfile for 
# Source directory: /Users/karel/git/modified-Fulgor
# Build directory: /Users/karel/git/modified-Fulgor/build-port
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test([=[phylign-port-regression]=] "/opt/homebrew/Frameworks/Python.framework/Versions/3.14/bin/python3.14" "/Users/karel/git/modified-Fulgor/tests/phylign_port_regression.py" "--new-binary" "/Users/karel/git/modified-Fulgor/build-port/fulgor")
set_tests_properties([=[phylign-port-regression]=] PROPERTIES  _BACKTRACE_TRIPLES "/Users/karel/git/modified-Fulgor/CMakeLists.txt;82;add_test;/Users/karel/git/modified-Fulgor/CMakeLists.txt;0;")
