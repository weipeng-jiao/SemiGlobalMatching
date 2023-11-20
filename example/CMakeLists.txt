
# message(STATUS ${SRC_LIST})
# 生成可执行
if (NON_ADD_LIB)
    message(STATUS "Just add exe")
    add_executable(${PROJECT_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp ${SRC_LIST})
    target_link_libraries(${PROJECT_NAME} ${OpenCV_LIBS} OpenMP::OpenMP_CXX ${CMAKE_THREAD_LIBS_INIT})
else()
    message(STATUS "exe linked src lib")
    add_executable(${PROJECT_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp )
    target_link_libraries(shared_lib)

endif()

# 打包文件
install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION bin)

