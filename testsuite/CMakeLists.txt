link_libraries(xc)
add_executable(xc-consistency xc-consistency.c)
add_executable(xc-error xc-error.c)
add_executable(xc-get_data xc-get_data.c)
add_executable(xc-regression xc-regression.c)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/workdir)

file(GLOB _polbz2s RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "regression/*.pol.bz2")
foreach(_polbz2 ${_polbz2s})
    # get test name
    get_filename_component(_polfile ${_polbz2} NAME)
    string(REPLACE "." ";" _pollist ${_polfile})
    separate_arguments(_pollist)
    list(GET _pollist 0 func)
    list(GET _pollist 1 system)
    set(refname "${func}.${system}")

    # extract reference values
    execute_process(COMMAND bunzip2 -c "${CMAKE_CURRENT_SOURCE_DIR}/regression/${refname}.pol.bz2"
                    OUTPUT_FILE "${CMAKE_CURRENT_BINARY_DIR}/workdir/${refname}.pol_ref")
    execute_process(COMMAND bunzip2 -c "${CMAKE_CURRENT_SOURCE_DIR}/regression/${refname}.unpol.bz2"
                    OUTPUT_FILE "${CMAKE_CURRENT_BINARY_DIR}/workdir/${refname}.unpol_ref")

    # this is properly the real test, but presumably it'll get rebuilt if libxc changes
    add_custom_command(TARGET xc-regression POST_BUILD
                       COMMAND xc-regression ARGS ${func} "${CMAKE_CURRENT_SOURCE_DIR}/input/${system}" "workdir/${func}.${system}.pol" "workdir/${func}.${system}.unpol" ">" "/dev/null")

    # the testing itself
        # Spin-unrestricted
    add_test(NAME "${func}-${system}-pol"
             COMMAND xc-error "workdir/${func}.${system}.pol" "workdir/${func}.${system}.pol_ref")
        # Spin-restricted
    add_test(NAME "${func}-${system}-unpol"
             COMMAND xc-error "workdir/${func}.${system}.unpol" "workdir/${func}.${system}.unpol_ref")
endforeach()

# Warning: I don't know how to make a test fail, so proper capture of failure not confirmed
