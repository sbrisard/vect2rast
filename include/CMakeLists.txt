# configure_file(
#   ${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}/${PROJECT_NAME}.hpp.in
#   ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}/${PROJECT_NAME}.hpp)

install(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}
  DESTINATION include)
#  PATTERN "*.in" EXCLUDE)

# install(
#   DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}
#   DESTINATION include)
