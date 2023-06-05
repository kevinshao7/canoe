include(FetchContent)

set(FETCHCONTENT_QUIET FALSE)

FetchContent_Declare(
  athenapp
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  URL https://github.com/chengcli/athenapp/archive/refs/tags/v0.3.1-alpha.tar.gz
)

FetchContent_MakeAvailable(athenapp)

include_directories(${athenapp_SOURCE_DIR})
