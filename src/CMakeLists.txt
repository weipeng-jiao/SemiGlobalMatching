# 遍历文件
file(GLOB src_cxx ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/*.c)
file(GLOB src_hxx ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
set( SRC_LIST ${src_cxx} ${src_hxx} CACHE PATH "lib include path" FORCE)


# 生成库文件 打包动态库
if (ADD_SHARED_LIB)
    message(STATUS "Add shared lib")
    add_library (shared_lib SHARED ${SRC_LIST})
    # install(TARGETS shared_lib LIBRARY DESTINATION lib)
endif()   

if (ADD_STATIC_LIB)
    message(STATUS "Add static lib")
    add_library (static_lib STATIC ${SRC_LIST})
    # install(TARGETS static_lib LIBRARY DESTINATION lib)
endif() 
