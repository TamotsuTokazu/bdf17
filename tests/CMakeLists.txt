# tests/CMakeLists.txt

# Download and include GoogleTest using FetchContent
include(FetchContent)

FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        v1.15.2 # Specify the desired version
)

# Fetch the content
FetchContent_MakeAvailable(googletest)

# Now, we can use the gtest and gtest_main targets

# Collect all test source files
file(GLOB TEST_SOURCES "*.cpp")

# Add the test executable
add_executable("${PROJECT_NAME}_tests" ${TEST_SOURCES})

# Link the test executable against gtest and gtest_main
target_link_libraries("${PROJECT_NAME}_tests"
    gtest_main
    # If your project executable is needed for testing, link it here
    ${PROJECT_NAME}
)

# Include directories if needed
target_include_directories("${PROJECT_NAME}_tests" PRIVATE ${PROJECT_SOURCE_DIR}/include)

# Add test cases to CTest
include(GoogleTest)
gtest_discover_tests("${PROJECT_NAME}_tests")
