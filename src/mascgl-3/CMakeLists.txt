set(APP_SRC
    ${APP_SRC}
    ${CMAKE_CURRENT_SOURCE_DIR}/model.cpp
    PARENT_SCOPE)

add_subdirectory(modelgraph)
add_subdirectory(mathtool)
add_subdirectory(soil)
