#ifdef __APPLE__
#else
#include <NBody/Visualization/GLLoader.hpp>
#include <iostream>
#include <dlfcn.h>

namespace NBody::Visualization {

#define GLE(ret, name, ...) name##proc* gl##name;
PAPAYA_GL_LIST
#undef GLE

bool glLoaderInit() {
  auto libGL = dlopen("libGL.so", RTLD_LAZY);
  if (!libGL) {
    std::cerr << "libGL.so cannot be loaded!\n";
    return false;
  }

#define GLE(ret, name, ...) \
  gl##name = (name##proc*)dlsym(libGL, "gl"#name); \
  if (!gl##name) { \
  std::cerr << "Function gl" << #name << " cannot be loaded from libGL.so!\n"; \
  return false; \
  }
  PAPAYA_GL_LIST
#undef GLE

  return true;
}

} // NBody::Visualization
#endif