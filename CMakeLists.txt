cmake_minimum_required(VERSION 3.16)
project(ecc_pari_generator CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Optimize
set(CMAKE_BUILD_TYPE Release)
add_compile_options(-O3 -march=native -flto)

# ===== PARI (manual) =====
find_library(PARI_LIB pari REQUIRED)
find_path(PARI_INCLUDE_DIR pari/pari.h REQUIRED)

# ===== libsodium =====
find_package(PkgConfig REQUIRED)
pkg_check_modules(SODIUM REQUIRED libsodium)

# ===== executable =====
add_executable(app
    main.cpp
    helpers.cpp
    random.cpp
    prime.cpp
    curve_search.cpp
)

# ===== include dirs =====
target_include_directories(app PRIVATE
    ${PARI_INCLUDE_DIR}
    ${SODIUM_INCLUDE_DIRS}
)

# ===== link libs =====
target_link_libraries(app PRIVATE
    ${PARI_LIB}
    ${SODIUM_LIBRARIES}
)