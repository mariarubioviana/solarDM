cmake_minimum_required(VERSION 3.18)  # Versión mínima requerida de CMake

project(SolarDM)  # Nombre del proyecto

# Configurar el ejecutable
add_executable(solarDM solarDM.cpp)  # Nombre del ejecutable y lista de archivos fuente

# Opcional: Configurar las opciones de compilación
set(CMAKE_CXX_STANDARD 14)  # Establecer el estándar de C++ (por ejemplo, C++11)

add_compile_options(-std=c++14 -Wno-c++11-extensions)

# Opcional: Agregar directorios de inclusión
# include_directories(path/to/include/directory)

# Opcional: Agregar bibliotecas

#target_link_libraries(fftw fftw3)

# Opcional: Configurar opciones de compilación adicionales
# target_compile_options(MyExecutable PRIVATE option1 option2)

# Opcional: Configurar opciones de enlace adicionales
# target_link_options(MyExecutable PRIVATE option1 option2)

# Target de instalación
install(TARGETS solarDM
        DESTINATION bin)
