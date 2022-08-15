option(ENABLE_ASAN "Enable AddressSanitizer" FALSE)
option(ENABLE_UBSAN "Enable UndefinedBehaviorSanitizer" FALSE)
option(ENABLE_TSAN "Enable ThreadSanitizer" FALSE)

if (ENABLE_ASAN)
  string(APPEND CMAKE_CXX_FLAGS " -fsanitize=address")
  string(APPEND CMAKE_EXE_LINKER_FLAGS " -fsanitize=address")
endif ()

if (ENABLE_UBSAN)
  string(APPEND CMAKE_CXX_FLAGS " -fsanitize=undefined")
  string(APPEND CMAKE_EXE_LINKER_FLAGS " -fsanitize=undefined")
endif ()

if (ENABLE_TSAN)
  string(APPEND CMAKE_CXX_FLAGS " -fsanitize=thread")
  string(APPEND CMAKE_EXE_LINKER_FLAGS " -fsanitize=thread")
endif ()