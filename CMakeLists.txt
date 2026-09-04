cmake_minimum_required(VERSION 3.25)
project(DioDe LANGUAGES C CXX)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(FetchContent)
set(GEOGRAM_SUB_BUILD ON CACHE BOOL "" FORCE)
set(GEOGRAM_WITH_GRAPHICS OFF CACHE BOOL "" FORCE)
set(GEOGRAM_WITH_LEGACY_NUMERICS OFF CACHE BOOL "" FORCE)
set(GEOGRAM_WITH_HLBFGS OFF CACHE BOOL "" FORCE)
set(GEOGRAM_WITH_TETGEN OFF CACHE BOOL "" FORCE)
set(GEOGRAM_WITH_TRIANGLE OFF CACHE BOOL "" FORCE)
set(GEOGRAM_WITH_LUA OFF CACHE BOOL "" FORCE)
set(GEOGRAM_LIB_ONLY ON CACHE BOOL "" FORCE)
if(NOT VORPALINE_PLATFORM)
    string(TOLOWER "${CMAKE_SYSTEM_PROCESSOR}" DIODE_SYSTEM_PROCESSOR)
    if(WIN32)
        if(MSVC)
            set(VORPALINE_PLATFORM Win-vs-generic)
        elseif(MINGW AND DIODE_SYSTEM_PROCESSOR MATCHES "^(amd64|x86_64)$")
            set(VORPALINE_PLATFORM Win64-gcc)
        else()
            message(FATAL_ERROR
                "Unsupported Windows toolchain/architecture for Geogram: "
                "${CMAKE_CXX_COMPILER_ID}/${CMAKE_SYSTEM_PROCESSOR}")
        endif()
    elseif(APPLE)
        if(NOT CMAKE_CXX_COMPILER_ID MATCHES "^(AppleClang|Clang)$")
            message(FATAL_ERROR
                "Geogram's bundled macOS platform requires Clang")
        elseif(DIODE_SYSTEM_PROCESSOR MATCHES "^(arm64|aarch64)$")
            set(VORPALINE_PLATFORM Darwin-aarch64-clang)
        elseif(DIODE_SYSTEM_PROCESSOR MATCHES "^(amd64|x86_64)$")
            set(VORPALINE_PLATFORM Darwin-clang)
        else()
            message(FATAL_ERROR
                "Unsupported macOS architecture for Geogram: ${CMAKE_SYSTEM_PROCESSOR}")
        endif()
    elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        if(DIODE_SYSTEM_PROCESSOR MATCHES "^(amd64|x86_64)$")
            if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
                set(VORPALINE_PLATFORM Linux64-gcc)
            elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
                set(VORPALINE_PLATFORM Linux64-clang)
            else()
                message(FATAL_ERROR
                    "Unsupported Linux compiler for Geogram: ${CMAKE_CXX_COMPILER_ID}")
            endif()
        elseif(DIODE_SYSTEM_PROCESSOR MATCHES "^(arm64|aarch64)$" AND
               CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
            set(VORPALINE_PLATFORM Linux64-gcc-aarch64)
        else()
            message(FATAL_ERROR
                "Unsupported Linux toolchain/architecture for Geogram: "
                "${CMAKE_CXX_COMPILER_ID}/${CMAKE_SYSTEM_PROCESSOR}. "
                "Set VORPALINE_PLATFORM explicitly if Geogram supports it.")
        endif()
    else()
        message(FATAL_ERROR
            "Unsupported platform for Geogram: ${CMAKE_SYSTEM_NAME}")
    endif()
endif()
FetchContent_Declare(
    geogram
    URL https://github.com/BrunoLevy/geogram/releases/download/v1.9.9/geogram_1.9.9.tar.gz
    URL_HASH SHA256=65402f3ce4b40efab178874c4ec1a8852f0ca1ab4567d72000d53efcaff52214
    EXCLUDE_FROM_ALL
)
FetchContent_MakeAvailable(geogram)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(geogram PRIVATE -w)
endif()
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    find_package(OpenMP REQUIRED COMPONENTS C CXX)
    target_link_libraries(geogram PUBLIC OpenMP::OpenMP_C OpenMP::OpenMP_CXX)
endif()

add_library(diode_backend INTERFACE)
target_include_directories(diode_backend INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(diode_backend INTERFACE geogram)

add_subdirectory(examples)
add_subdirectory(bindings/python)
