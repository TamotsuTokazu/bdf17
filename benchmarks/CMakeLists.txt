# benchmarks/CMakeLists.txt

# Include FetchContent to download dependencies
include(FetchContent)

# Declare Google Benchmark
FetchContent_Declare(
    googlebenchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG        v1.9.0 # Specify the desired version
)

# Optionally, include Google Test as a dependency for Google Benchmark's testing
set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)

# Fetch and make available Google Benchmark
FetchContent_MakeAvailable(googlebenchmark)

# Collect all benchmark source files
file(GLOB BENCHMARK_SOURCES "*.cpp")

# Add the benchmark executable
add_executable("${PROJECT_NAME}_benchmarks" ${BENCHMARK_SOURCES})

# Link the benchmark executable against Google Benchmark and your project executable or libraries
target_link_libraries("${PROJECT_NAME}_benchmarks"
    benchmark::benchmark_main
    ${PROJECT_NAME} # Link against your main project executable or its libraries
)

# Include directories if needed
target_include_directories("${PROJECT_NAME}_benchmarks" PRIVATE ${PROJECT_SOURCE_DIR}/include)

# Optionally, set compile options for benchmarks
target_compile_options("${PROJECT_NAME}_benchmarks" PRIVATE -Ofast -march=native -mavx2 -funroll-loops -ftree-vectorize -Wno-unused-parameter)

# Installation rules for benchmarks (optional)
# install(TARGETS gpv23_benchmarks DESTINATION bin)
