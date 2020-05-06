# Add all available projects

# Macro from https://stackoverflow.com/questions/7787823/cmake-how-to-get-the-name-of-all-subdirectories-of-a-directory
macro(subdirlist result curdir)
    file(GLOB children RELATIVE ${curdir} ${curdir}/*)
    set(dirlist "")
    foreach(child ${children})
        if(IS_DIRECTORY ${curdir}/${child})
            LIST(APPEND dirlist ${child})
        endif()
    endforeach()
    set(${result} ${dirlist})
endmacro()

subdirlist(LoS ${CMAKE_CURRENT_SOURCE_DIR})
foreach(subdir ${LoS})
    add_subdirectory(${subdir})
endforeach()
