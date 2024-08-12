if(EXISTS "/home/kqs/canoe/examples/install_manifest.txt")
    file(READ "/home/kqs/canoe/examples/install_manifest.txt" files)
    string(REGEX REPLACE "\n" ";" files "${files}")
    foreach(file ${files})
        message(STATUS "Uninstalling \"${file}\"")
        if(EXISTS "${file}")
            execute_process(
                COMMAND /usr/bin/cmake -E remove "${file}"
                OUTPUT_VARIABLE rm_out
                RESULT_VARIABLE rm_retval
            )
            if(NOT rm_retval EQUAL "0")
                message(FATAL_ERROR "Problem when removing \"${file}\"")
            endif()
        else()
            message(STATUS "File \"${file}\" does not exist.")
        endif()
    endforeach()
else()
    message(STATUS "Install manifest not found.")
endif()
