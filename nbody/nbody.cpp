#include <NBody/NBodyProgram.hpp>

int main(int argc, char* argv[]) try {
  // Marco FLOATING_POINT_TYPE is configured in CMakeLists.txt.
  const auto program = NBody::NBodyProgram<FLOATING_POINT_TYPE>{
      argc,
      argv
  };
  program.launch();
  return 0;
} catch (const std::exception& e) {
  std::cerr << e.what() << std::endl;
  return 1;
}